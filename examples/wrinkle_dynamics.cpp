#include "halfedge/trimesh_types.h"
#include <argus/argus.hpp>
#include <argus/util.hpp>
#include <cstdlib>
#include <glm/glm.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <set>
#include <stdexcept>
#include <sys/stat.h>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
#include <random>
#include <chrono>

// External Libraries
#include <LBFGS.h>
#include <LBFGSB.h>
#include <boost/filesystem.hpp>
#include <igl/readOBJ.h>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>
#include <Discregrid/All>
#include <Eigen/Geometry>
std::vector<int> clamppp;
bool saveCheckpoints = false;
struct RenderedMeshData
{
	std::mutex mutex;
	std::vector<std::string> names;
	std::vector<Eigen::Matrix<argus::Real, Eigen::Dynamic, 3>> faceNormals;
	std::vector<Eigen::Matrix<argus::Real, Eigen::Dynamic, 3>> faceGradPhis;
	std::vector<Eigen::Matrix<argus::Real, Eigen::Dynamic, 3>> vertexNormals;
	std::vector<Eigen::Matrix<argus::Real, Eigen::Dynamic, 3>> meshVs;
	std::vector<std::vector<std::vector<trimesh::index_t>>> meshFs;
	std::vector<argus::VectorN> meshStrains;
	std::vector<argus::VectorN> meshAmps;
	std::vector<argus::VectorN> meshSDfs;
};

struct StaticMeshData
{
	std::vector<std::string> names;
	std::vector<Eigen::Matrix<argus::Real, Eigen::Dynamic, 3>> meshVs;
	std::vector<std::vector<std::vector<trimesh::index_t>>> meshFs;
};

struct RenderedPCData
{
	std::mutex mutex;
	std::vector<std::string> names;
	std::vector<Eigen::Matrix<argus::Real, Eigen::Dynamic, 3>> vertices;
};

struct SimContext
{
	nlohmann::json simData;
	int totalFrames = 1000;
	std::atomic<bool> keepRunning { true };
	std::atomic<bool> isPaused { true }; // Start paused?

	std::unique_ptr<RenderedMeshData> clothData;
	std::vector<std::unique_ptr<RenderedMeshData>> dynamicObstacles;
	std::unique_ptr<StaticMeshData> staticMeshes;
	std::unique_ptr<StaticMeshData> staticObstacles;
	std::unique_ptr<RenderedPCData> renderSeams;
	std::unique_ptr<RenderedPCData> renderClamps;

	SimContext()
	{
		clothData = std::make_unique<RenderedMeshData>();
		staticMeshes = std::make_unique<StaticMeshData>();
		staticObstacles = std::make_unique<StaticMeshData>();
		renderSeams = std::make_unique<RenderedPCData>();
		renderClamps = std::make_unique<RenderedPCData>();
	}
};

std::vector<std::vector<trimesh::index_t>> convertTrianglesToPolyscope(const std::vector<trimesh::triangle_t>& list)
{
	std::vector<std::vector<trimesh::index_t>> result;
	result.reserve(list.size());
	for (const auto& l : list)
		result.push_back({ l.i(), l.j(), l.k() });
	return result;
}

Eigen::Matrix<argus::Real, Eigen::Dynamic, 3> getFaceGradPhis(argus::Sheet& sheet)
{
	Eigen::Matrix<argus::Real, Eigen::Dynamic, 3> gradphis;
	const std::vector<trimesh::triangle_t>& faces = sheet.getTriangles();
	gradphis.resize(sheet.getFaceCount(), 3);
	for (int i = 0; i < sheet.getFaceCount(); i++)
	{
		argus::Vector3 vi = sheet.getPositions().segment<3>(3 * faces[i].v[0]);
		argus::Vector3 vj = sheet.getPositions().segment<3>(3 * faces[i].v[1]);
		argus::Vector3 vk = sheet.getPositions().segment<3>(3 * faces[i].v[2]);

		argus::Matrix32 A;
		A.col(0) = vj - vi;
		A.col(1) = vk - vi;

		argus::Matrix22 ATA = A.transpose() * A;
		gradphis.row(i) = A * ATA.inverse() * sheet.getDphis1VectorPerFace().col(i);
	}
	return gradphis;
}

class Visualizer
{
public:
	explicit Visualizer(SimContext& ctx)
	    : ctx(ctx)
	{
	}

	void run(argus::Sheet& sheet)
	{
		polyscope::init();

		setupStaticScene(sheet);

		polyscope::state::userCallback = [this]()
		{ this->callback(); };

		polyscope::show();

		ctx.keepRunning = false;
	}

private:
	SimContext& ctx;
	int currentFrame = 0;

	// UI State for Exporter
	int ui_exportFrame = 0;
	int ui_targetMeshIdx = 0; // 0 for Sheet, 1+ for Obstacles
	std::string ui_statusMsg = "Ready";
	std::string m_exportBaseDir = "exports/";

	void ensureExportDir()
	{
		if (!fs::exists(m_exportBaseDir))
		{
			fs::create_directories(m_exportBaseDir);
		}
	}

	// Export Logic
	void exportToOBJ(int frame, int targetIdx)
	{
		ensureExportDir();
		std::string filename;
		std::string statusPrefix;

		std::lock_guard<std::mutex> lock(ctx.clothData->mutex);

		// Handle "Sheet" (Index 0)
		if (targetIdx == 0)
		{
			if (frame >= ctx.clothData->meshVs.size())
			{
				ui_statusMsg = "Error: Frame " + std::to_string(frame) + " not buffered.";
				return;
			}
			filename = m_exportBaseDir + "sheet_f" + std::to_string(frame) + ".obj";
			writeOBJFile(filename, ctx.clothData->meshVs[frame], ctx.clothData->meshFs[frame]);
		}
		// Handle Obstacles (Index 1+)
		else
		{
			int obsIdx = targetIdx - 1;
			if (obsIdx >= ctx.dynamicObstacles.size())
			{
				ui_statusMsg = "Error: Invalid Mesh Index.";
				return;
			}
			auto& obs = ctx.dynamicObstacles[obsIdx];
			std::lock_guard<std::mutex> obsLock(obs->mutex);
			if (frame >= obs->meshVs.size())
			{
				ui_statusMsg = "Error: Obstacle frame not buffered.";
				return;
			}
			filename = m_exportBaseDir + "obs_" + std::to_string(obsIdx) + "_f" + std::to_string(frame) + ".obj";
			writeOBJFile(filename, obs->meshVs[frame], obs->meshFs[frame]);
		}

		ui_statusMsg = "Exported: " + filename;
	}

	// Use a template for the face container to accept long, size_t, etc.
	template <typename T>
	void writeOBJFile(const std::string& path,
	    const Eigen::MatrixXd& V,
	    const std::vector<std::vector<T>>& F)
	{
		std::ofstream out(path);
		if (!out)
		{
			// Assuming this is inside your class with ui_statusMsg
			ui_statusMsg = "Error: Cannot write to " + path;
			return;
		}

		// Write Vertices
		for (int i = 0; i < V.rows(); i++)
		{
			out << "v " << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << "\n";
		}

		// Write Faces
		for (const auto& face : F)
		{
			out << "f";
			for (auto idx : face)
			{
				// OBJ is 1-indexed, so add 1
				out << " " << (idx + 1);
			}
			out << "\n";
		}
		out.close();
	}

	void drawExporterUI()
	{
		ImGui::Separator();
		ImGui::Text("--- Scene Exporter ---");

		// List Mappings
		if (ImGui::TreeNode("View Mesh Mappings"))
		{
			ImGui::BulletText("Index 0: Sheet (Cloth)");
			for (size_t i = 0; i < ctx.dynamicObstacles.size(); ++i)
			{
				ImGui::BulletText("Index %d: Obstacle %zu", (int)i + 1, i);
			}
			ImGui::TreePop();
		}
		ImGui::SetNextItemWidth(300.0f);
		ImGui::InputInt("Mesh Index", &ui_targetMeshIdx);
		ImGui::SetNextItemWidth(300.0f);

		ImGui::InputInt("Frame to Export", &ui_exportFrame);

		if (ImGui::Button("Export Single Frame"))
		{
			exportToOBJ(ui_exportFrame, ui_targetMeshIdx);
		}
		ImGui::SameLine();
		if (ImGui::Button("Export All Frames"))
		{
			size_t total = (ui_targetMeshIdx == 0) ? ctx.clothData->meshVs.size() : ctx.dynamicObstacles[ui_targetMeshIdx - 1]->meshVs.size();
			for (int f = 0; f < (int)total; ++f)
			{
				exportToOBJ(f, ui_targetMeshIdx);
			}
			ui_statusMsg = "Exported " + std::to_string(total) + " frames.";
		}

		ImGui::TextColored(ImVec4(1, 1, 0, 1), "Status: %s", ui_statusMsg.c_str());
	}

	void setupStaticScene(argus::Sheet& sheet)
	{
		int id = 0;
		for (const auto& name : ctx.staticMeshes->names)
		{
			auto* m = polyscope::registerSurfaceMesh(name, ctx.staticMeshes->meshVs[id], ctx.staticMeshes->meshFs[id]);
			m->setEnabled(false);
			id++;
		}

		id = 0;
		for (const auto& name : ctx.staticObstacles->names)
		{
			auto* m = polyscope::registerSurfaceMesh(name, ctx.staticObstacles->meshVs[id], ctx.staticObstacles->meshFs[id]);
			m->setEnabled(false);
			id++;
		}

		id = 0;
		for (const auto& name : ctx.renderClamps->names)
		{
			auto* p = polyscope::registerPointCloud(name, ctx.renderClamps->vertices[id]);
			p->setEnabled(false);
			id++;
		}

		// setupUVLayout(sheet);

		// setupSeams(sheet);

		polyscope::updateStructureExtents();
		auto [minB, maxB] = polyscope::state::boundingBox;
		glm::vec3 center = 0.5f * (minB + maxB);
		float padding = 1.2f * glm::length(maxB - minB);
		polyscope::state::boundingBox = std::make_tuple(center - glm::vec3(padding), center + glm::vec3(padding));
		polyscope::options::automaticallyComputeSceneExtents = false;
	}

	void setupUVLayout(argus::Sheet& sheet)
	{
		argus::MatrixNN V = sheet.getPositionsN3().transpose();
		Eigen::Vector3d min_v = V.colwise().minCoeff();
		Eigen::Vector3d max_v = V.colwise().maxCoeff();
		double mesh_height = std::max(max_v.y() - min_v.y(), 1.0);
		double mesh_width = std::max(max_v.x() - min_v.x(), 1.0);

		argus::MatrixNN TC = sheet.getUVs().transpose();
		auto uvtris = sheet.getUVTriangles();
		Eigen::MatrixXd TC_3D(TC.rows(), 3);
		TC_3D.col(0) = TC.col(0) * mesh_height;
		TC_3D.col(1) = TC.col(1) * mesh_height;
		TC_3D.col(2).setZero();

		auto* uvMesh = polyscope::registerSurfaceMesh("UV_Layout", TC_3D, convertTrianglesToPolyscope(uvtris));
		uvMesh->setEnabled(false);
		uvMesh->setEdgeWidth(1.0);
		uvMesh->setSurfaceColor({ 0.2, 0.8, 0.2 });
		double offset_x = max_v.x() + (mesh_width * 1.5);
		uvMesh->setPosition(glm::vec3(offset_x, min_v.y(), min_v.z()));
	}

	void setupSeams(argus::Sheet& sheet)
	{
		argus::MatrixNN V = sheet.getPositionsN3().transpose();
		auto uvtris = sheet.getUVTriangles();
		Eigen::MatrixXi F = sheet.getFaceMatrix();
		std::vector<std::set<int>> v_to_tc_map(V.rows());
		std::vector<Eigen::Vector3d> seamVertices;

		for (int i = 0; i < F.rows(); ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				v_to_tc_map[F(i, j)].insert(uvtris[i].v[j]);
			}
		}
		for (int i = 0; i < V.rows(); ++i)
		{
			if (v_to_tc_map[i].size() > 1)
				seamVertices.push_back(V.row(i));
		}

		if (!seamVertices.empty())
		{
			Eigen::MatrixXd seamCloud(seamVertices.size(), 3);
			for (size_t i = 0; i < seamVertices.size(); ++i)
				seamCloud.row(i) = seamVertices[i];

			auto* pc = polyscope::registerPointCloud("Seam_Vertices", seamCloud);
			pc->setPointColor({ 1.0, 0.0, 0.0 });
			pc->setPointRadius(0.02 * (polyscope::state::lengthScale));
		}
	}

	bool updateSceneToFrame(int frame)
	{
		std::lock_guard<std::mutex> lock(ctx.clothData->mutex);
		if (frame >= static_cast<int>(ctx.clothData->meshVs.size()))
			return false;

		// Sheet
		auto* mesh = polyscope::registerSurfaceMesh("Sheet", ctx.clothData->meshVs[frame], ctx.clothData->meshFs[frame]);
		if (frame < ctx.clothData->meshAmps.size())
			mesh->addVertexScalarQuantity("Amplitudes", ctx.clothData->meshAmps[frame]);
		if (frame < ctx.clothData->meshStrains.size())
			mesh->addFaceScalarQuantity("Strain", ctx.clothData->meshStrains[frame]);
		if (frame < ctx.clothData->faceGradPhis.size())
			mesh->addFaceVectorQuantity("GradPhis", ctx.clothData->faceGradPhis[frame]);

		// Dynamic Obstacles
		int iter = 0;
		for (auto& obsData : ctx.dynamicObstacles)
		{
			std::lock_guard<std::mutex> obsLock(obsData->mutex);
			if (frame < obsData->names.size())
			{
				polyscope::registerSurfaceMesh(obsData->names[frame], obsData->meshVs[frame], obsData->meshFs[frame]);
			}
			iter++;
		}

		// Seams
		std::lock_guard<std::mutex> seamLock(ctx.renderSeams->mutex);
		if (frame < ctx.renderSeams->vertices.size() && !ctx.renderSeams->vertices.empty())
		{
			polyscope::registerPointCloud(ctx.renderSeams->names[frame], ctx.renderSeams->vertices[frame]);
		}

		// Clamps
		std::lock_guard<std::mutex> clampLock(ctx.renderClamps->mutex);
		if (frame < ctx.renderClamps->vertices.size() && !ctx.renderClamps->vertices.empty())
		{
			polyscope::registerPointCloud(ctx.renderClamps->names[frame], ctx.renderClamps->vertices[frame]);
		}
		return true;
	}

	void callback()
	{
		ImGui::PushItemWidth(100);
		bool frameBuffered = true;

		if (ImGui::Button("Next"))
			frameBuffered = updateSceneToFrame(++currentFrame);
		ImGui::SameLine();
		if (ImGui::Button("Prev"))
			frameBuffered = updateSceneToFrame(currentFrame = std::max(0, currentFrame - 1));

		if (ctx.isPaused)
		{
			if (ImGui::Button("Play"))
				ctx.isPaused = false;
		}
		else
		{
			frameBuffered = updateSceneToFrame(currentFrame);
			if (frameBuffered)
				currentFrame = (currentFrame + 1) % ctx.totalFrames;
			if (ImGui::Button("Pause"))
				ctx.isPaused = true;
		}

		ImGui::SameLine();
		ImGui::Text("Frame: %d", currentFrame);

		// Slider
		int tempFrame = currentFrame;
		ImGui::SetNextItemWidth(400.0f);
		if (ImGui::SliderInt("##Frame", &tempFrame, 0, ctx.totalFrames - 1))
		{
			size_t bufferedFrames = 0;
			{
				std::lock_guard<std::mutex> lock(ctx.clothData->mutex);
				bufferedFrames = ctx.clothData->meshVs.size();
			}

			int maxFrame = (bufferedFrames > 0) ? static_cast<int>(bufferedFrames - 1) : 0;
			if (tempFrame > maxFrame)
			{
				currentFrame = maxFrame;
				frameBuffered = (bufferedFrames > 0) ? updateSceneToFrame(currentFrame) : false;
				// Keep play state when the requested frame isn't buffered yet.
			}
			else
			{
				currentFrame = tempFrame;
				frameBuffered = updateSceneToFrame(currentFrame);
				ctx.isPaused = true;
			}
		}
		if (!frameBuffered)
		{
			ImGui::SameLine();
			ImGui::Text("Frame not buffered yet");
		}
		ui_exportFrame = currentFrame; // Default to current frame
		drawExporterUI();

		ImGui::PopItemWidth();
	}
};

// --- SIMULATION LOGIC ---
void checkNans(argus::Sheet& sheet)
{
	//  if (!vec.allFinite()) {
	//     std::cout << "Vector contains NaN or Inf values." << std::endl;
	if (sheet.getPositions().allFinite() == false)
	{
		std::cerr << "Error nan positions!!";
		if (saveCheckpoints)
			exit(1);
	}
	if (sheet.getDphis1VectorPerFace().allFinite() == false)
	{
		std::cerr << "Error nan dphis!!";
		if (saveCheckpoints)

			exit(1);
	}
	if (sheet.getAmpVec().allFinite() == false)
	{
		std::cerr << "Error nan amps!!";
		if (saveCheckpoints)

			exit(1);
	}
	if (sheet.getVelocities().allFinite() == false)
	{
		std::cerr << "Error nan velocities!!";
		if (saveCheckpoints)

			exit(1);
	}
}
void runSimulationLoop(argus::World& world, argus::Sheet& sheet, argus::SceneAnimator& animator, SimContext& ctx)
{
	if (!ctx.simData.contains("timestepper_config"))
		throw std::runtime_error("No timestepper_config found!");

	nlohmann::json params = ctx.simData["timestepper_config"];
	world.prepare(params);

	auto t1 = argus::util::captureCurrentTime();

	if (ctx.simData["obstacles_present"])
	{
		for (auto& obs : world.getObstacles())
			obs->startAsyncLoad(0);
	}

	animator.update(0, 0.0, sheet, world.getObstacles());

	for (int i = 0; i < ctx.totalFrames; i++)
	{
		if (false && i == 300)
		{
			sheet.getClampedVertices().resize(clamppp.size());
			for (int cv = 0; cv < clamppp.size(); cv++)
				sheet.getClampedVertices()(cv) = clamppp[cv];
		}
		if (!ctx.keepRunning)
			break;

		std::cout << "--- Frame " << i << " / " << ctx.totalFrames << " ---" << std::endl;

		if (ctx.simData["obstacles_present"])
		{
			for (auto& obs : world.getObstacles())
				obs->updateToActive();
		}

		sheet.updateCurFrame(i);
		animator.update(i, ctx.simData["dt"].get<argus::Real>(), sheet, world.getObstacles());
		world.step(ctx.simData["dt"].get<argus::Real>(), i);
		checkNans(sheet);
		for (auto& obs : world.getObstacles())
			obs->startAsyncLoad(i + 1);

		// --- Data Sync to Context ---
		{
			std::lock_guard<std::mutex> lock(ctx.clothData->mutex);

			ctx.clothData->meshVs.push_back(sheet.getPositionsN3().transpose());
			ctx.clothData->meshFs.push_back(convertTrianglesToPolyscope(sheet.getTriangles()));
			ctx.clothData->meshAmps.push_back(sheet.getAmpVec());

			sheet.computeGreenStrainPerFace();
			ctx.clothData->meshStrains.push_back(sheet.getStrainValuesPerFace());
			ctx.clothData->faceGradPhis.push_back(getFaceGradPhis(sheet));
		}

		// Sync Dynamic Obstacles
		if (!world.getObstacles().empty())
		{
			int iter = 0;

			for (auto&& it : world.getObstacles())
			{
				if (it->meshAvailable())
				{
					auto [V, F] = it->getMesh();
					if (V.rows() != 0 && iter < ctx.dynamicObstacles.size())
					{
						std::lock_guard<std::mutex> obsLock(ctx.dynamicObstacles[iter]->mutex);
						ctx.dynamicObstacles[iter]->names.push_back(it->getName());
						ctx.dynamicObstacles[iter]->meshVs.push_back(V);
						ctx.dynamicObstacles[iter]->meshFs.push_back(F);
						iter++;
					}
				}
			}
		}

		// Sync Seams/Clamps
		std::map<int, int> seams = sheet.getSeamVertices();
		if (!seams.empty())
		{
			Eigen::Matrix<argus::Real, Eigen::Dynamic, 3> seamVerts(2 * seams.size(), 3);
			int counter = 0;
			for (auto it = seams.begin(); it != seams.end(); it++)
			{
				seamVerts.row(counter) = sheet.getPositionsN3().col(it->first);
				seamVerts.row(counter + 1) = sheet.getPositionsN3().col(it->second);
				counter += 2;
			}
			std::lock_guard<std::mutex> lock(ctx.renderSeams->mutex);
			ctx.renderSeams->names.push_back("Seam vertices");
			ctx.renderSeams->vertices.push_back(seamVerts);
		}

		Eigen::VectorXi clamp = sheet.getClampedVertices();
		if (clamp.size() > 0)
		{
			Eigen::Matrix<argus::Real, Eigen::Dynamic, 3> clampedVerts(clamp.size(), 3);
			for (int k = 0; k < clamp.size(); k++)
				clampedVerts.row(k) = sheet.getPositionsN3().col(clamp(k));

			std::lock_guard<std::mutex> lock(ctx.renderClamps->mutex);
			ctx.renderClamps->names.push_back("Clamped vertices");
			ctx.renderClamps->vertices.push_back(clampedVerts);
		}
	}

	auto t2 = argus::util::captureCurrentTime();
	std::cout << "Simulation finished in "
	          << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()
	          << " ms" << std::endl;
}

// --- HEADLESS SIMULATION LOOP ---

int runSimulationHeadless(argus::World& world, argus::Sheet& sheet, argus::SceneAnimator& animator, SimContext& ctx)
{
	if (!ctx.simData.contains("timestepper_config"))
		throw std::runtime_error("Error: timestepper_config missing!");

	nlohmann::json params = ctx.simData["timestepper_config"];
	world.prepare(params);

	if (saveCheckpoints)
	{
		world.saveCheckpoint({ "positions", "energy", "velocities", "amps",
		                         "dphisPerFace", "stvk_energy_per_face", "collision_vertices" },
		    argus::project_root + "checkpoints/", 0);
	}

	auto t1 = argus::util::captureCurrentTime();

	// --- Frame 0 Initialization ---
	std::cout << "Initializing obstacles..." << std::endl;
	if (ctx.simData["obstacles_present"])
	{
		for (auto& obs : world.getObstacles())
			obs->startAsyncLoad(0);
	}

	// Apply animation at Frame 0
	animator.update(0, 0.0, sheet, world.getObstacles());

	if (ctx.simData["obstacles_present"])
	{
		for (auto& obs : world.getObstacles())
		{
			obs->updateToActive();
		}
	}

	// --- Simulation Loop ---
	int totalFrames = ctx.simData["total_frames"].get<int>();

	// Run stabilization frame (Frame 1) if desired
	// std::cout << "Running frame 1 / " << totalFrames << std::endl;
	// sheet.updateCurFrame(1);
	// world.step(ctx.simData["dt"], 1);

	// if (saveCheckpoints)
	// {
	// 	world.saveCheckpoint({ "positions", "energy", "velocities", "amps",
	// 	                         "dphisPerFace", "stvk_energy_per_face", "collision_vertices" },
	// 	    argus::project_root + "checkpoints/", 1);
	// }

	// // Preload Frame 2 if obstacles exist
	// if (ctx.simData["obstacles_present"])
	// {
	// 	for (auto& obs : world.getObstacles())
	// 		obs->startAsyncLoad(2);
	// }

	// Main Loop (Starting from 2)
	for (unsigned int i = 0; i < static_cast<unsigned int>(totalFrames); i++)
	{
		if (false && i == 300)
		{
			sheet.getClampedVertices().resize(clamppp.size());
			for (int cv = 0; cv < clamppp.size(); cv++)
				sheet.getClampedVertices()(cv) = clamppp[cv];
		}
		std::cout << "--- Frame " << i << " / " << totalFrames << " ---" << std::endl;

		if (ctx.simData["obstacles_present"])
		{
			for (auto& obs : world.getObstacles())
				obs->updateToActive();
		}

		sheet.updateCurFrame(i);
		animator.update(i, ctx.simData["dt"].get<argus::Real>(), sheet, world.getObstacles());

		world.step(ctx.simData["dt"], i);

		// Preload Next
		for (auto& obs : world.getObstacles())
			obs->startAsyncLoad(i + 1);

		if (saveCheckpoints)
		{
			world.saveCheckpoint({ "positions", "energy", "velocities", "amps",
			                         "dphisPerFace", "stvk_energy_per_face", "collision_vertices" },
			    argus::project_root + "checkpoints/", i);
		}
		checkNans(sheet);
	}

	auto t2 = argus::util::captureCurrentTime();
	std::cout << "End simulation loop." << std::endl;

	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
	if (saveCheckpoints)
	{
		nlohmann::json stats;
		try
		{
			stats = world.getTimeStepStats();
		}
		catch (...)
		{
			std::cerr << "Error! Your timestepper does not provide getStats()!\n";
		}
		stats["total_runtime"] = duration.count();

		argus::simConf["optimization_stats"] = stats;
	}

	std::cout << "Simulation took " << duration.count() << " ms" << std::endl;
	std::cout << "##### Successful execution #####\n";

	return duration.count();
}

// --- SETUP HELPERS ---

void setConstraints(argus::Sheet& sheet)
{
	argus::Real sheet_thickness = sheet.getThickness();
	auto amp_constraints = [sheet_thickness](argus::VectorN& X, argus::VectorN& Xdot, std::vector<int> dofs, argus::Real dt)
	{
		if (saveCheckpoints)
		{
			if (argus::simConf.contains("dynamic_amp_constraints") == false)
				argus::simConf["dynamic_amp_constraints"] = "On";
		}
		argus::Real threshold = 0.1 * sheet_thickness;
		for (auto it = dofs.begin(); it != dofs.end(); it++)
		{
			if (X(*it) < threshold)
			{
				Xdot(*it) = (threshold - X(*it)) / dt;
				X(*it) = threshold;
			}
		}
	};

	std::vector<int> dofs(sheet.getVertexCount());
	std::iota(dofs.begin(), dofs.end(), 0);
	sheet.addConstraint(argus::DynamicConstraint { dofs, amp_constraints });
}

void initializeSheet(argus::Sheet& sheet, nlohmann::json& simData)
{
	sheet.getRestShape().requireDeformedCoordsToDeformationGradient();
	sheet.getRestShape().requireIntrinsicGeometry();
	sheet.buildDampingLaplacian();
	sheet.updateBaseMesh();

	std::vector<int> active_dofs(4 * sheet.getVertexCount() + 2 * sheet.getFaceCount());
	std::iota(active_dofs.begin(), active_dofs.end(), 0);
	sheet.setActiveDofs(active_dofs);

	if (simData["cloth"].contains("init_mesh_path"))
	{
		auto [tri, pos, uv] = argus::util::loadMeshFromObj(argus::project_root + simData["cloth"]["init_mesh_path"].get<std::string>());
		sheet.getPositionsN3() = pos;
	}
	// sheet.getPositions() = 1 * argus::Vector3::Random();
	// sheet.getPositionsN3().col(195) += argus::Vector3({ 0.02, 0, 0 });
}

std::tuple<std::vector<trimesh::triangle_t>, argus::Matrix3N, std::vector<trimesh::triangle_t>, argus::Matrix2N> preprocessMesh(nlohmann::json data)
{
	std::vector<trimesh::triangle_t> triangles, uvtriangles;
	argus::Matrix3N pos;
	argus::Matrix2N uv;

	if (data.contains("reference_mesh_path"))
		std::tie(triangles, pos, uvtriangles, uv) = argus::util::loadMeshWithUVFromObj(data["reference_mesh_path"].get<std::string>());
	else if (data.contains("mesh_path"))
		std::tie(triangles, pos, uvtriangles, uv) = argus::util::loadMeshWithUVFromObj(argus::project_root + data["mesh_path"].get<std::string>());
	else
		throw std::runtime_error("Mesh path information is not passed!");
	if (data.contains("motion"))
	{
		for (auto&& it : data["motion"])
		{
			std::string mode = it["mode"].get<std::string>();
			if (mode == "t")
			{
				argus::Vector3 update = it["updates"].get<argus::Vector3>();
				for (int i = 0; i < pos.cols(); i++)
					pos.col(i) += update;
			}
			else if (mode == "s")
			{
				argus::Vector3 update = it["updates"].get<argus::Vector3>();
				argus::Matrix33 scale = argus::Matrix33::Identity();
				scale(0, 0) = update(0);
				scale(1, 1) = update(1);
				scale(2, 2) = update(2);
				for (int i = 0; i < pos.cols(); i++)
					pos.col(i) = (scale * pos.col(i)).eval();
			}
			else if (mode == "r")
			{
				argus::Real theta = M_PI * it["theta"].get<argus::Real>() / 180.0;
				std::string axis = it["axis"].get<std::string>();
				argus::Matrix33 rotationMat;
				if (axis == "x")
				{
					rotationMat << 1, 0, 0,
					    0, std::cos(theta), -std::sin(theta),
					    0, std::sin(theta), std::cos(theta);
				}
				else if (axis == "y")
				{
					rotationMat << std::cos(theta), 0, -std::sin(theta),
					    0, 1, 0,
					    std::sin(theta), 0, std::cos(theta);
				}
				else if (axis == "z")
				{
					rotationMat << std::cos(theta), -std::sin(theta), 0,
					    std::sin(theta), std::cos(theta), 0,
					    0, 0, 1;
				}
				else
				{
					throw "Not implemented!";
				}
				pos = (rotationMat * pos).eval();
			}
		}
	}
	return std::make_tuple(triangles, pos, uvtriangles, uv);
}

void addObstaclesToContext(argus::World& world, SimContext& ctx)
{
	if (!ctx.simData.contains("obstacles_present") || !ctx.simData["obstacles_present"])
		return;

	auto& obstacleData = ctx.simData["obstacles"];
	for (auto it = obstacleData.begin(); it != obstacleData.end(); ++it)
	{
		std::string type = (*it)["type"];
		auto obstacleOpt = argus::getObstacle(type, *it);

		if (obstacleOpt.has_value())
		{
			argus::Obstacle* rawObs = obstacleOpt.value().get();
			world.addObstacle(std::move(obstacleOpt.value()));
			if (rawObs->meshAvailable())
				try
				{

					auto [V, F] = rawObs->getMesh();
					std::string name = (*it).value("name", type);
					if (V.rows() != 0)
					{
						// Create RenderData for this dynamic obstacle
						auto obd = std::make_unique<RenderedMeshData>();
						obd->names.push_back(name);
						obd->meshVs.push_back(V);
						obd->meshFs.push_back(F);
						ctx.dynamicObstacles.push_back(std::move(obd));
					}
				}
				catch (...)
				{
				}
		}
	}
}

int main(int argc, char** argv)
{
	SimContext ctx;
	bool guiMode = true;

	std::string config_path;

	for (int i = 1; i < argc; ++i)
	{
		std::string arg = argv[i];
		if (arg == "--nogui")
		{
			guiMode = false;
		}
		else if (arg == "--gui")
		{
			guiMode = true;
		}
		else if (arg == "--checkpoint")
		{
			saveCheckpoints = true;
		}
		else if (!arg.empty() && arg[0] == '-')
		{
			std::cerr << "Unknown flag: " << arg << std::endl;
			return 1;
		}
		else if (config_path.empty())
		{
			config_path = argus::project_root + arg;
		}
		else
		{
			std::cerr << "Unexpected extra argument: " << arg << std::endl;
			return 1;
		}
	}

	if (config_path.empty())
	{
		std::cerr << "Usage: ./wrinkle_dynamics <config_path> [--gui|--nogui] [--checkpoint]" << std::endl;
		return 1;
	}
	clamppp = { 4, 81, 89, 194, 218, 296, 298, 309, 644, 1745, 2249, 2318, 2331, 2332, 2341, 2342, 2950, 2951, 2963, 2964, 2969, 2970, 2987, 3001, 4084, 4085, 4086, 4139, 4156, 4169, 4744, 5326, 5327, 5328, 5373, 5916, 6546, 7050, 7133, 7143, 7751, 7752, 7764, 7765, 7770, 7771, 7802, 8885, 8886, 8887, 8940, 8957, 8970, 9545, 10127, 10128, 10129, 10174, 10717, 11347, 11851, 11854, 11861, 11920, 11921, 11933, 11934, 11943, 11944, 12507, 12508, 12552, 12553, 12565, 12566, 12571, 12572, 12589, 12590, 12591, 12603, 13686, 13687, 13688, 13741, 13758, 13771, 14346, 14378, 14527, 14539, 14929, 14930, 15518, 21742, 21743, 21747, 21748, 21749, 21750, 21751, 21752, 23607, 23608, 23609, 23637, 23638, 23767, 23787, 23788, 23789, 23818, 23819, 23826, 23827, 23828, 23833, 24027, 24028, 24029, 24114, 24115, 24116, 24144, 24145, 24146, 24174, 24175, 24176, 24211, 24291, 24292, 24293, 24330, 24331, 24332, 24384, 24385, 24386, 24396, 24397, 24398, 25089, 25090, 25091, 25107, 25108, 25109 };

	std::cout << "Loading config: " << config_path << std::endl;
	ctx.simData = argus::util::loadSimulationConfig(config_path);
	ctx.totalFrames = ctx.simData.value("total_frames", 1000);

	argus::World world;
	argus::SceneAnimator sceneAnimator;
	sceneAnimator.load(ctx.simData);
	std::vector<trimesh::triangle_t> triangles, uvtriangles;
	argus::Matrix3N positions;
	argus::Matrix2N uvs;

	try
	{
		std::tie(triangles, positions, uvtriangles, uvs) = preprocessMesh(ctx.simData["cloth"]);
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << std::endl;
		return 1;
	}

	argus::Material material(ctx.simData["rho"], ctx.simData["Y"], ctx.simData["mu"], ctx.simData["h"]);

	argus::SheetInfo info;
	if (ctx.simData["cloth"].contains("mesh_path"))
		info.sheetPath = argus::project_root + ctx.simData["cloth"]["mesh_path"].get<std::string>();
	if (ctx.simData["cloth"].contains("seam_vertices_path"))
		info.seam_path = argus::project_root + ctx.simData["cloth"]["seam_vertices_path"].get<std::string>();

	if (ctx.simData["cloth"].contains("clamped_vertices"))
	{
		auto& arr = ctx.simData["cloth"]["clamped_vertices"];
		if (arr.is_array())
		{
			for (auto& item : arr)
			{
				if (item.is_array() && item.size() == 4)
				{
					info.clampedVertices.insert(info.clampedVertices.end(), { item[0].get<int>(), item[1].get<int>(), item[2].get<int>(), item[3].get<int>() });
				}
			}
		}
	}
	else if (ctx.simData["cloth"].contains("clamped_vertices_path"))
	{
		std::string p = argus::project_root + ctx.simData["cloth"]["clamped_vertices_path"].get<std::string>();
		std::ifstream f(p);
		if (f.is_open())
		{
			int v, x, y, z;
			while (f >> v >> x >> y >> z)
				info.clampedVertices.insert(info.clampedVertices.end(), { v, x, y, z });
		}
	}
	if (ctx.simData["cloth"].contains("clamped_vertices_motion_path"))
	{

		std::string p = argus::project_root + ctx.simData["cloth"]["clamped_vertices_motion_path"].get<std::string>();
		argus::AnimationSequence anim = argus::parseKinematicTrack(p);

		info.clampedVerticesMotion = anim;
	}

	bool wrinkles_present = ctx.simData.value("sim_type", "wrinkles") == "wrinkles";
	info.wrinkle_mesh = wrinkles_present;
	std::unique_ptr<argus::Sheet> sheetPtr = std::make_unique<argus::Sheet>(
	    material, triangles, positions, uvtriangles, uvs, info);
	argus::Sheet& sheet = static_cast<argus::Sheet&>(world.addDynamic(std::move(sheetPtr)));

	if (ctx.simData.contains("forces"))
	{
		std::vector<nlohmann::json> forcesData = ctx.simData["forces"];
		std::optional<std::vector<std::shared_ptr<argus::SheetForces>>> sheetForces = argus::getSheetForcesVector(forcesData);
		sheet.clearSheetForces();
		if (sheetForces.has_value())
			sheet.addSheetForces(sheetForces.value());
	}

	setConstraints(sheet);
	initializeSheet(sheet, ctx.simData);

	addObstaclesToContext(world, ctx);

	std::string timestepper = ctx.simData["timestepper"].get<std::string>();

	world.setIntegrator(timestepper);

	int runtimeMs = -1;
	if (guiMode)
	{
		{
			std::lock_guard<std::mutex> lock(ctx.clothData->mutex);
			ctx.clothData->meshVs.push_back(sheet.getPositionsN3().transpose());
			ctx.clothData->meshFs.push_back(convertTrianglesToPolyscope(sheet.getTriangles()));
			ctx.clothData->meshAmps.push_back(sheet.getAmpVec());

			sheet.computeGreenStrainPerFace();
			ctx.clothData->meshStrains.push_back(sheet.getStrainValuesPerFace());
			ctx.clothData->faceGradPhis.push_back(getFaceGradPhis(sheet));

			ctx.staticMeshes->names.push_back("Rest Mesh");
			ctx.staticMeshes->meshVs.push_back(sheet.getRestPositionsN3().transpose());
			ctx.staticMeshes->meshFs.push_back(convertTrianglesToPolyscope(sheet.getTriangles()));
		}

		std::cout << "Starting simulation thread..." << std::endl;
		auto guiStart = argus::util::captureCurrentTime();
		std::thread simThread(runSimulationLoop, std::ref(world), std::ref(sheet), std::ref(sceneAnimator), std::ref(ctx));

		std::cout << "Starting visualizer..." << std::endl;
		Visualizer visualizer(ctx);
		visualizer.run(sheet);

		if (simThread.joinable())
			simThread.join();
		auto guiEnd = argus::util::captureCurrentTime();
		runtimeMs = static_cast<int>(std::chrono::duration_cast<std::chrono::milliseconds>(guiEnd - guiStart).count());
	}
	else
	{
		if (saveCheckpoints)
		{
			argus::initCheckpointing(sheet);
			std::cout << "Saving input configuration!\n";
			std::string cpDir = argus::project_root + "/checkpoints/" + sheet.getCheckpointDirName();
			boost::filesystem::create_directories(cpDir);
			std::ofstream o(cpDir + "/inputConfig.json");
			o << std::setw(4) << ctx.simData << std::endl;
			o.close();
		}
		runtimeMs = runSimulationHeadless(world, sheet, sceneAnimator, ctx);
	}

	if (!guiMode)
	{
		if (saveCheckpoints)
		{
			argus::simConf["total_time_wrinkle_optimization_ms"] = runtimeMs;
			argus::simConf["total_frames"] = ctx.simData["total_frames"];
			std::cout << "Saving simulation configuration!\n";
			std::string cpDirEnd = argus::project_root + "/checkpoints/" + sheet.getCheckpointDirName();
			std::ofstream op(cpDirEnd + "/simulationConfig.json");
			op << std::setw(4) << argus::simConf << std::endl;
			op.close();
		}
	}

	return 0;
}
