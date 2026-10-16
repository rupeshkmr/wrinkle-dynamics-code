#include <argus/sheet_state_io.hpp>
#include <argus/util.hpp>
namespace argus
{
	void saveSheetObj(std::string saveDir, int frame, Sheet& sheet, argus::Real& energy)
	{
		Matrix3N positions = sheet.getPositionsN3();
		const argus::Matrix2N uvs = sheet.getUVs();
		const std::vector<trimesh::triangle_t> triangles = sheet.getTriangles();
		std::ofstream file;
		std::string savePath = saveDir + "/frame." + std::to_string(frame) + ".obj";
		file.open(savePath);
		for (int i = 0; i < positions.cols(); i++)
		{
			file << "v " << positions.coeff(0, i) << " " << positions.coeff(1, i) << " " << positions.coeff(2, i) << "\n";
		}
		for (int i = 0; i < uvs.cols(); i++)
		{
			file << "vt " << uvs.coeff(0, i) << " " << uvs.coeff(1, i) << "\n";
		}
		for (int i = 0; i < triangles.size(); i++)
		{
			file << "f ";
			for (int j = 0; j < 3; j++)
			{
				file << triangles[i].v[j] + 1 << "/" << triangles[i].v[j] + 1 << " ";
			}
			file << "\n";
		}
		file.close();

		// saving metadata (no need to load the energy)
		std::string metadataSavepath = saveDir + "/metadata.txt." + std::to_string(frame); // to save energy and iteration
		file.open(metadataSavepath);
		file << std::to_string(frame) << "\n";
		file << std::to_string(energy) << "\n";
		file.close();
	}

	void saveSheet(const std::string fileName, const Sheet& sheet)
	{
		if (!std::filesystem::exists("./output"))
		{
			std::filesystem::create_directory("./output");
		}
		std::ofstream file("./output/" + fileName + ".obj");
		if (file.is_open())
		{
			const Matrix3N& positions = sheet.getPositionsN3();
			for (size_t i = 0; i < positions.cols(); i++)
			{
				file << "v " << positions(0, i) << " " << positions(1, i) << " " << positions(2, i) << std::endl;
			}
			const std::vector<trimesh::triangle_t>& triangles = sheet.getTriangles();
			for (size_t i = 0; i < triangles.size(); i++)
			{
				file << "f " << triangles[i].i() + 1 << " " << triangles[i].j() + 1 << " " << triangles[i].k() + 1 << std::endl;
			}
			file.close();
		}
	}
	void initCheckpointing(Sheet& sheet)
	{
		/**
		 * Setup the folder where checkpoints will be stored.
		 */
#ifndef ARGUS_CHECKPOINT
		std::cerr << "Error! Trying to create checkpoint folder when no checkpoints are requested! Exiting. " << std::endl;
		exit(0);
#endif
		std::string saveDir = argus::project_root + "checkpoints/" + sheet.getCheckpointDirName();
		std::cout << "Initiated checkpointing to " << saveDir << std::endl;
		boost::filesystem::create_directories(saveDir);
	}
	void saveCheckpoint(std::vector<std::string> data_to_save, std::string parent_directory, int frame, argus::Sheet& sheet, nlohmann::json state)
	{
		// Only save checkpoints if required
#ifdef ARGUS_CHECKPOINT
		std::map<std::string, int>& m = sheet.getParamsToSave();
		m.clear();
		for (auto it = data_to_save.begin(); it != data_to_save.end(); it++)
			sheet.addParamToSave(*it);
		saveSheetState(parent_directory + sheet.getCheckpointDirName(), frame, sheet, &m, state);
#endif
	}

	void saveSheetState(std::string saveDir, int frame, Sheet& sheet, std::map<std::string, int>* parametersToSave, nlohmann::json state)
	{
		// taken from https://github.com/AleksandarHaber/Save-and-Load-Eigen-Cpp-Matrices-Arrays-to-and-from-CSV-files/blob/master/source_file.cpp
		Matrix3N positions = sheet.getPositionsN3();
		Eigen::MatrixXi faces = sheet.getFaceMatrix();
#ifdef ARGUS_DEBUG
		std::cout << "Save Dir " << saveDir << std::endl;
#endif
		// make sure it is a directory
		if (saveDir[saveDir.length() - 1] != '/')
			saveDir = saveDir + "/";

		boost::filesystem::create_directories(saveDir);
		const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
		// create directory that contains the positions
		boost::filesystem::create_directories(saveDir + "positions/");
		boost::filesystem::create_directories(saveDir + "velocities/");
		std::string positionsSavepath = saveDir + "positions/positions.csv." + std::to_string(frame);
		std::string velSavepath = saveDir + "velocities/velocities.csv." + std::to_string(frame);
		if (frame == 0)
		{
			boost::filesystem::create_directories(saveDir + "restPositions/");
			boost::filesystem::create_directories(saveDir + "faces/");
			boost::filesystem::create_directories(saveDir + "clampedVertices/");
			if (sheet.useReferenceMesh())
			{
				boost::filesystem::create_directories(saveDir + "reference_mesh/");
				std::string refPosSavePath = saveDir + "reference_mesh/positions.csv";
				util::writeEigenData(refPosSavePath, sheet.getReferenceMeshPositions());
				std::string refFaceSavePath = saveDir + "reference_mesh/faces.csv";
				util::writeEigenData(refFaceSavePath, util::convertTrianglesToMatrix(sheet.getReferenceMeshTriangles()));
			}
			std::string restPositionsSavepath = saveDir + "restPositions/restpositions.csv";
			util::writeEigenData(restPositionsSavepath, sheet.getRestPositionsN3());
			std::string facesSavepath = saveDir + "faces/faces.csv";
			util::writeEigenData(facesSavepath, faces);
			std::string clampedVertsSavepath = saveDir + "clampedVertices/clamped_vertices.txt";
			util::writeEigenData(clampedVertsSavepath, sheet.getClampedVertices());
		}

		if (state.contains("collision_vertices"))
		{
			std::string collVertsSavepath = saveDir + "collisionVertices/";
			boost::filesystem::create_directories(collVertsSavepath);
			if (state["collision_vertices"] != 0)
			{
				Eigen::VectorXi v = state["collision_vertices"];
				util::writeEigenData(collVertsSavepath + "collision.csv." + std::to_string(frame), v);
			}
		}
		if ((*parametersToSave)["stvk_energy_per_face"])
		{
			sheet.computeGreenStrainPerFace();
			VectorN stvk_per_face = sheet.getStrainValuesPerFace();
			std::string strain_dir = saveDir + "stvk/";
			boost::filesystem::create_directories(strain_dir);
			util::writeEigenData(strain_dir + "strain.csv." + std::to_string(frame), stvk_per_face);
		}

		if ((*parametersToSave)["energy"])
		{
			std::string energy_folder = saveDir + "energies/";
			boost::filesystem::create_directories(energy_folder);
			std::ofstream outfile(energy_folder + "energies.txt", std::ios::app);
			if (outfile.is_open())
			{
				outfile << std::to_string(sheet.getEnergies()) << std::endl;
				outfile.close();
			}
			else
			{
				std::cerr << "Error in saving sheet energies\n"
				          << std::endl;
			}
		}

		if ((*parametersToSave)["positions"])
			util::writeEigenData(positionsSavepath, positions);

		if ((*parametersToSave)["velocities"])
			util::writeEigenData(velSavepath, sheet.getVelocitiesN3());
		if (parametersToSave)
		{
			if ((*parametersToSave)["initamps"])
			{
				VectorN initamplitudes = sheet.getAmplitudes();
				boost::filesystem::create_directories(saveDir + "initamplitudes/");
				std::string initampSavepath = saveDir + "initamplitudes/initamplitudes.csv";
				util::writeEigenData(initampSavepath, initamplitudes);
			}
			if ((*parametersToSave)["initdphisPerFace"])
			{
				Matrix2N initdphis = sheet.getDphis1VectorPerFace();
				boost::filesystem::create_directories(saveDir + "initdphisPerFace/");
				std::string initdphiSavepath = saveDir + "initdphisPerFace/initdphisPerFace.csv";
				util::writeEigenData(initdphiSavepath, initdphis);
			}

			if ((*parametersToSave)["amps"])
			{
				VectorN amplitudes = sheet.getAmplitudes();
				boost::filesystem::create_directories(saveDir + "amplitudes/");
				std::string ampSavepath = saveDir + "amplitudes/amplitudes.csv." + std::to_string(frame);
				util::writeEigenData(ampSavepath, amplitudes);
			}
			if ((*parametersToSave)["dphisPerFace"])
			{
				Matrix2N dphis = sheet.getDphis1VectorPerFace();
				boost::filesystem::create_directories(saveDir + "dphisPerFace/");
				std::string dphiSavepath = saveDir + "dphisPerFace/dphisPerFace.csv." + std::to_string(frame);
				util::writeEigenData(dphiSavepath, dphis);
			}
			if ((*parametersToSave)["dphisPerVertex"])
			{
				Matrix3N vtexdphis = sheet.getDphisPerVertex();
				boost::filesystem::create_directories(saveDir + "dphisPerVertex/");
				std::string dphiVtexSavepath = saveDir + "dphisPerVertex/dphisPerVertex.csv." + std::to_string(frame);
				util::writeEigenData(dphiVtexSavepath, vtexdphis);
			}
			if ((*parametersToSave)["phis"])
			{
				VectorN phis = sheet.getPhis();
				boost::filesystem::create_directories(saveDir + "phis/");
				std::string phiSavepath = saveDir + "phis/phis.csv." + std::to_string(frame);
				util::writeEigenData(phiSavepath, phis);
			}

			// wme.sheet.addParamToSave("cornerPhis");
			if ((*parametersToSave)["cornerPhis"])
			{
				Matrix3N cornerPhis = sheet.getCornerPhis();
				boost::filesystem::create_directories(saveDir + "cornerPhis/");
				std::string cornerphiSavepath = saveDir + "cornerPhis/cornerPhis.csv." + std::to_string(frame);
				util::writeEigenData(cornerphiSavepath, cornerPhis);
			}
			// wme.sheet.addParamToSave("wrinkleMeshPos");
			if ((*parametersToSave)["wrinkledMeshPos"])
			{
				Matrix3N wpos = sheet.getWrinkleMeshPositions();
				boost::filesystem::create_directories(saveDir + "wrinkledMeshPositions/");
				std::string wrinklemeshposSavepath = saveDir + "wrinkledMeshPositions/wrinkledMeshPositions.csv." + std::to_string(frame);
				util::writeEigenData(wrinklemeshposSavepath, wpos);
			}
			// wme.sheet.addParamToSave("wrinkleMeshFaces");
			if ((*parametersToSave)["wrinkledMeshFaces"])
			{
				Eigen::MatrixXi wfaces = sheet.getWrinkleMeshFaces();
				boost::filesystem::create_directories(saveDir + "wrinkledMeshFaces/");
				std::string wrinklemeshfacesSavepath = saveDir + "wrinkledMeshFaces/wrinkledMeshFaces.csv." + std::to_string(frame);
				util::writeEigenData(wrinklemeshfacesSavepath, wfaces);
			}
		}
	}

	bool loadSheetState(std::string loadDir, int frame, Sheet& sheet)
	{
		std::string ampSavepath = loadDir + "amplitudes/amplitudes.csv." + std::to_string(frame);
		std::string dphiSavepath = loadDir + "dphisPerFace/dphisPerFace.csv." + std::to_string(frame);
		std::string positionsSavepath = loadDir + "positions/positions.csv." + std::to_string(frame);
		std::string dphiPerVertexSavepath = loadDir + "dphisPerVertex/dphisPerVertex.csv";
		std::string cornerPhisSavepath = loadDir + "cornerPhis/cornerPhis.csv";
		std::string facesSavepath = loadDir + "faces/faces.csv";
		std::string restPositionssavepath = loadDir + "restPositions/restpositions.csv";
		std::string wrinkledMeshFacesSavepath = loadDir + "wrinkledMeshFaces/wrinkledMeshFaces.csv";
		std::string wrinkledMeshPositionsSavepath = loadDir + "wrinkledMeshPositions/wrinkledMeshPositions.csv";
		std::cout << positionsSavepath << std::endl;
		try
		{
			int nfaces = sheet.getTriangles().size();
			int nvertices = sheet.getPositionsN3().cols();
			argus::MatrixNN pos, amps, dphis, dphiPerVertex, cornerPhis, restPositions, wrinkledMeshPositions;
			Eigen::MatrixXi faces, wrinkledMeshFaces;
			bool status;
			status = util::openEigenData(positionsSavepath, &pos);
			if (status)
				sheet.getPositionsN3() = pos.reshaped(3, sheet.getPositionsN3().cols());
			else
				return false;
			status = util::openEigenData(ampSavepath, &amps);
			if (status)
				sheet.getAmplitudes() = amps.reshaped();
			status = util::openEigenData(dphiSavepath, &dphis);
			if (status)
				sheet.getDphis1VectorPerFace() = dphis.reshaped(2, sheet.getTriangles().size());
			status = util::openEigenData(dphiPerVertexSavepath, &dphiPerVertex);
			if (status)
			{
				sheet.getCheckpointParams()["dphisPerVertex"] = 1;
				sheet.getDphisPerVertex() = dphiPerVertex.reshaped(3, nvertices);
			}
			status = util::openEigenData(cornerPhisSavepath, &cornerPhis);
			if (status)
			{
				sheet.getCheckpointParams()["cornerPhis"] = 1;
				sheet.getCornerPhis() = cornerPhis.reshaped(3, nfaces);
			}
			status = util::openEigenData(facesSavepath, &faces);
			// TODO: add the functionality to convert from face matrix to face vector

			status = util::openEigenData(restPositionssavepath, &restPositions);
			if (status)
				sheet.getRestPositionsN3() = restPositions;

			status = util::openEigenData(wrinkledMeshFacesSavepath, &wrinkledMeshFaces);
			if (status)
			{
				sheet.getCheckpointParams()["wrinkledMeshFaces"] = 1;
				sheet.getWrinkleMeshFaces() = wrinkledMeshFaces;
			}
			status = util::openEigenData(wrinkledMeshPositionsSavepath, &wrinkledMeshPositions);
			if (status)
			{
				sheet.getCheckpointParams()["wrinkledMeshPositions"] = 1;
				sheet.getWrinkleMeshPositions() = wrinkledMeshPositions;
			}
			std::cout << "Loaded frame " << frame << std::endl;
			return true;
		}
		catch (...)
		{
			std::cout << "Failed to load frame " << frame << std::endl;
			throw "Error while loading sheet data!\n";
		}
		return true;
	}

};
