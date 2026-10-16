#ifndef SCENE_ANIMATOR_HPP
#define SCENE_ANIMATOR_HPP

#include "common.hpp"
#include "sheet.hpp"
#include "obstacles.hpp"
#include <nlohmann/json.hpp>
#include <Eigen/Geometry>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

using json = nlohmann::json;
namespace fs = std::filesystem;
namespace argus
{

	// Maps a Vertex ID to its 3D Position for a single frame
	using FrameData = std::unordered_map<int, Vector3>;

	// Maps a Frame Number to its corresponding FrameData
	using AnimationSequence = std::map<int, FrameData>;

	inline AnimationSequence parseKinematicTrack(const std::string& filepath)
	{
		AnimationSequence animSequence;
		std::ifstream file(filepath);

		if (!file.is_open())
		{
			std::cerr << "Error: Could not open tracking file: " << filepath << "\n";
			return animSequence;
		}

		std::string line;
		int currentFrame = -1;

		while (std::getline(file, line))
		{
			// Trim leading whitespace (handles cross-platform line endings gracefully)
			line.erase(0, line.find_first_not_of(" \t\r\n"));

			// Skip empty lines or comments
			if (line.empty() || line[0] == '#')
			{
				continue;
			}

			// Check if the line denotes a new frame
			if (line.rfind("Frame", 0) == 0)
			{
				std::istringstream iss(line);
				std::string prefix;
				iss >> prefix >> currentFrame;
			}
			// Otherwise, parse the vertex data
			else if (currentFrame != -1)
			{
				std::istringstream iss(line);
				int v_idx;
				Vector3 pos;

				// Extract ID, X, Y, Z
				if (iss >> v_idx >> pos(0) >> pos(1) >> pos(2))
				{
					animSequence[currentFrame][v_idx] = pos;
				}
			}
		}

		file.close();
		return animSequence;
	}

	struct TransformKeyframe
	{
		int frame;
		Eigen::Vector3d translation;
		Eigen::Quaterniond rotation;
		Eigen::Vector3d scale;

		bool hasTranslation = false;
		bool hasRotation = false;
		bool hasScale = false;
	};

	class AnimationTrack
	{
	public:
		std::string targetName;
		std::vector<int> vertexIndices;

		// NEW: Store the initial positions to avoid cumulative drift
		std::vector<Eigen::Vector3d> restPositions;

		std::vector<TransformKeyframe> keyframes;
		Eigen::Vector3d centerOfRotation;
		bool useCenterOfRotation = false;

		void addKeyframe(int frame, const Eigen::Vector3d& t, const Eigen::Quaterniond& r, const Eigen::Vector3d& s, bool ht, bool hr, bool hs)
		{
			TransformKeyframe k;
			k.frame = frame;
			k.translation = t;
			k.rotation = r;
			k.scale = s;
			k.hasTranslation = ht;
			k.hasRotation = hr;
			k.hasScale = hs;
			keyframes.push_back(k);

			std::sort(keyframes.begin(), keyframes.end(), [](const TransformKeyframe& a, const TransformKeyframe& b)
			    { return a.frame < b.frame; });
		}

		void evaluate(int currentFrame, Eigen::Vector3d& outT, Eigen::Quaterniond& outR, Eigen::Vector3d& outS) const
		{
			outT = Eigen::Vector3d::Zero();
			outR = Eigen::Quaterniond::Identity();
			outS = Eigen::Vector3d(1, 1, 1);

			if (keyframes.empty())
				return;

			if (currentFrame <= keyframes.front().frame)
			{
				const auto& kf = keyframes.front();
				if (kf.hasTranslation)
					outT = kf.translation;
				if (kf.hasRotation)
					outR = kf.rotation;
				if (kf.hasScale)
					outS = kf.scale;
				return;
			}

			if (currentFrame >= keyframes.back().frame)
			{
				const auto& kf = keyframes.back();
				if (kf.hasTranslation)
					outT = kf.translation;
				if (kf.hasRotation)
					outR = kf.rotation;
				if (kf.hasScale)
					outS = kf.scale;
				return;
			}

			const TransformKeyframe* p0 = &keyframes.front();
			const TransformKeyframe* p1 = &keyframes.back();

			for (size_t i = 0; i < keyframes.size() - 1; ++i)
			{
				if (currentFrame >= keyframes[i].frame && currentFrame < keyframes[i + 1].frame)
				{
					p0 = &keyframes[i];
					p1 = &keyframes[i + 1];
					break;
				}
			}

			double t = (double)(currentFrame - p0->frame) / (double)(p1->frame - p0->frame);

			Eigen::Vector3d t0 = p0->hasTranslation ? p0->translation : Eigen::Vector3d::Zero();
			Eigen::Vector3d t1 = p1->hasTranslation ? p1->translation : Eigen::Vector3d::Zero();
			outT = t0 * (1.0 - t) + t1 * t;

			Eigen::Vector3d s0 = p0->hasScale ? p0->scale : Eigen::Vector3d(1, 1, 1);
			Eigen::Vector3d s1 = p1->hasScale ? p1->scale : Eigen::Vector3d(1, 1, 1);
			outS = s0 * (1.0 - t) + s1 * t;

			Eigen::Quaterniond r0 = p0->hasRotation ? p0->rotation : Eigen::Quaterniond::Identity();
			Eigen::Quaterniond r1 = p1->hasRotation ? p1->rotation : Eigen::Quaterniond::Identity();
			outR = r0.slerp(t, r1);
		}
	};

	class SceneAnimator
	{
	private:
		std::vector<AnimationTrack> tracks;
		bool m_isActive = false;

	public:
		SceneAnimator()
		    : m_isActive(false)
		{
		}

		void load(const json& config)
		{
			tracks.clear();
			m_isActive = false;

			std::cout << "[SceneAnimator] Loading configuration..." << std::endl;

			if (!config.contains("animations"))
			{
				std::cerr << "[SceneAnimator] ERROR: JSON missing 'animations' array." << std::endl;
				return;
			}

			for (const auto& anim : config["animations"])
			{
				AnimationTrack track;

				if (anim.contains("target_obstacle"))
				{
					track.targetName = anim["target_obstacle"];
				}
				// --- 1. Load Vertex Indices (File OR Array) ---
				if (anim.contains("indices_path"))
				{
					std::string path = argus::project_root + anim["indices_path"].get<std::string>();
					std::ifstream file(path);
					if (file.is_open())
					{
						int idx;
						while (file >> idx)
						{
							track.vertexIndices.push_back(idx);
						}
					}
					else
					{
						std::cerr << "[SceneAnimator] Error opening indices file: " << path << std::endl;
					}
				}
				else if (anim.contains("indices")) // Check for direct array
				{
					for (auto& idx : anim["indices"])
					{
						track.vertexIndices.push_back(idx);
					}
				}

				// --- 2. Load Center (Applied if ANY indices were loaded) ---
				if (!track.vertexIndices.empty() && anim.contains("center"))
				{
					track.centerOfRotation = Eigen::Vector3d(
					    anim["center"][0], anim["center"][1], anim["center"][2]);
					track.useCenterOfRotation = true;
				}

				// if (anim.contains("indices_path"))
				// {
				// 	std::string path = argus::project_root + anim["indices_path"].get<std::string>();
				// 	std::ifstream file(path);
				// 	if (file.is_open())
				// 	{
				// 		int idx;
				// 		while (file >> idx)
				// 		{
				// 			track.vertexIndices.push_back(idx);
				// 		}
				// 	}
				// 	else
				// 	{
				// 		std::cerr << "[SceneAnimator] Error opening indices file: " << path << std::endl;
				// 	}

				// 	if (anim.contains("center"))
				// 	{
				// 		track.centerOfRotation = Eigen::Vector3d(
				// 		    anim["center"][0], anim["center"][1], anim["center"][2]);
				// 		track.useCenterOfRotation = true;
				// 	}
				// }
				// else if (anim.contains("vertex_indices"))
				// {

				// 	track.vertexIndices.push_back(idx);
				// }

				if (anim.contains("keyframes"))
				{
					for (const auto& kf : anim["keyframes"])
					{
						int f = kf["frame"];
						Eigen::Vector3d t(0, 0, 0);
						Eigen::Quaterniond r = Eigen::Quaterniond::Identity();
						Eigen::Vector3d s(1, 1, 1);
						bool ht = false, hr = false, hs = false;

						if (kf.contains("translation"))
						{
							t = Eigen::Vector3d(kf["translation"][0], kf["translation"][1], kf["translation"][2]);
							ht = true;
						}
						if (kf.contains("rotation_quat"))
						{
							r = Eigen::Quaterniond(kf["rotation_quat"][3], kf["rotation_quat"][0], kf["rotation_quat"][1], kf["rotation_quat"][2]);
							hr = true;
						}
						else if (kf.contains("rotation_euler"))
						{
							double x = kf["rotation_euler"][0];
							double y = kf["rotation_euler"][1];
							double z = kf["rotation_euler"][2];
							r = Eigen::AngleAxisd(x * M_PI / 180.0, Eigen::Vector3d::UnitX())
							    * Eigen::AngleAxisd(y * M_PI / 180.0, Eigen::Vector3d::UnitY())
							    * Eigen::AngleAxisd(z * M_PI / 180.0, Eigen::Vector3d::UnitZ());
							hr = true;
						}
						if (kf.contains("scale"))
						{
							s = Eigen::Vector3d(kf["scale"][0], kf["scale"][1], kf["scale"][2]);
							hs = true;
						}
						track.addKeyframe(f, t, r, s, ht, hr, hs);
					}
				}

				tracks.push_back(track);
			}

			if (!tracks.empty())
				m_isActive = true;
		}

		void update(int frame, Real dt, Sheet& sheet, std::vector<std::unique_ptr<Obstacle>>& obstacles)
		{
			if (!m_isActive)
				return;

			// Iterate using reference so we can modify 'track' (specifically restPositions)
			for (auto& track : tracks)
			{
				Eigen::Vector3d t_cur;
				Eigen::Quaterniond r_cur;
				Eigen::Vector3d s_cur;
				track.evaluate(frame, t_cur, r_cur, s_cur);

				if (!track.targetName.empty())
				{
					bool found = false;
					for (auto& obs : obstacles)
					{
						if (obs->getName() == track.targetName)
						{
							Vector3 final_t = t_cur.cast<Real>();
							Eigen::Quaternion final_r(r_cur.w(), r_cur.x(), r_cur.y(), r_cur.z());
							Vector3 final_s = s_cur.cast<Real>();
							obs->setTransform(final_t, final_r, final_s, dt);
							found = true;
						}
					}
					if (!found && frame == 0)
					{
						std::cerr << "[SceneAnimator] WARNING: Could not find obstacle named '"
						          << track.targetName << "' in scene!" << std::endl;
					}
				}

				// --- 2. Vertex/Sheet Animation ---
				if (!track.vertexIndices.empty())
				{
					// FIX: Capture Rest Positions on the first run
					if (track.restPositions.empty())
					{
						track.restPositions.reserve(track.vertexIndices.size());
						for (int idx : track.vertexIndices)
						{
							if (idx < sheet.getVertexCount())
							{
								// Store the initial position as Vector3d
								track.restPositions.push_back(sheet.getPositions().segment<3>(idx * 3).cast<double>());
							}
							else
							{
								track.restPositions.push_back(Eigen::Vector3d::Zero()); // Padding for invalid idx
							}
						}
					}

					// Prepare Transform
					Vector3 center = track.centerOfRotation.cast<Real>();
					Eigen::Quaternion<Real> rot = r_cur.cast<Real>();
					Vector3 trans = t_cur.cast<Real>(); // Support translation too

					// Apply transformation to REST positions
					for (size_t i = 0; i < track.vertexIndices.size(); ++i)
					{
						int idx = track.vertexIndices[i];
						if (idx < sheet.getVertexCount())
						{
							Vector3 pos_rest = track.restPositions[i].cast<Real>();

							// P_new = Center + Rot * (P_rest - Center) + Translation
							Vector3 newPos = center + rot * (pos_rest - center) + trans;

							sheet.getPositions().segment<3>(idx * 3) = newPos;

							// Important: If you are using a Physics Solver,
							// ensure these vertices have inverse_mass = 0,
							// otherwise the solver will fight this position update.
						}
					}
				}
			}
		}
	};
}
#endif