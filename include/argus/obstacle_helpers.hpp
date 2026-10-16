#ifndef OBSTACLE_HELPERS_HPP
#define OBSTACLE_HELPERS_HPP

#include "common.hpp"
#include <vector>
#include <string>
#include <unordered_map>
#include <fstream>
#include <iostream>

namespace argus
{
	class Obstacle;

	static probeWeights weights10 = probeWeights {

		{ { 0.333f, 0.333f, 0.333f },

		    { 0.666f, 0.167f, 0.167f },

		    { 0.167f, 0.666f, 0.167f },

		    { 0.167f, 0.167f, 0.666f },

		    { 0.500f, 0.500f, 0.000f },

		    { 0.500f, 0.000f, 0.500f },

		    { 0.000f, 0.500f, 0.500f },

		    { 1.0f, 0.0f, 0.0f },

		    { 0.0f, 1.0f, 0.0f },

		    { 0.0f, 0.0f, 1.0f },

		    { 0.0f, 0.0f, 0.0f },

		    { 0.0f, 0.0f, 0.0f },

		    { 0.0f, 0.0f, 0.0f },

		    { 0.0f, 0.0f, 0.0f },

		    { 0.0f, 0.0f, 0.0f },

		    { 0.0f, 0.0f, 0.0f },

		    { 0.0f, 0.0f, 0.0f },

		    { 0.0f, 0.0f, 0.0f },

		    { 0.0f, 0.0f, 0.0f },

		    { 0.0f, 0.0f, 0.0f },

		    { 0.0f, 0.0f, 0.0f },

		    { 0.0f, 0.0f, 0.0f },

		    { 0.0f, 0.0f, 0.0f },

		    { 0.0f, 0.0f, 0.0f },

		    { 0.0f, 0.0f, 0.0f } }

	};

	static probeWeights weights15 = probeWeights { { // Vertices (3)

		{ 1.000f, 0.000f, 0.000f },

		{ 0.000f, 1.000f, 0.000f },

		{ 0.000f, 0.000f, 1.000f },

		// Mid-edges (3)

		{ 0.500f, 0.500f, 0.000f },

		{ 0.500f, 0.000f, 0.500f },

		{ 0.000f, 0.500f, 0.500f },

		// Barycenter (1)

		{ 0.333f, 0.333f, 0.333f },

		// Inner "Ring" - closer to center (3)

		{ 0.400f, 0.400f, 0.200f },

		{ 0.400f, 0.200f, 0.400f },

		{ 0.200f, 0.400f, 0.400f },

		// Outer "Ring" - closer to vertices (3)

		{ 0.700f, 0.150f, 0.150f },

		{ 0.150f, 0.700f, 0.150f },

		{ 0.150f, 0.150f, 0.700f },

		// Mid-inner (2 additional)

		{ 0.250f, 0.250f, 0.500f },

		{ 0.500f, 0.250f, 0.250f },

		{ 0.0f, 0.0f, 0.0f },

		{ 0.0f, 0.0f, 0.0f },

		{ 0.0f, 0.0f, 0.0f },

		{ 0.0f, 0.0f, 0.0f },

		{ 0.0f, 0.0f, 0.0f },

		{ 0.0f, 0.0f, 0.0f },

		{ 0.0f, 0.0f, 0.0f },

		{ 0.0f, 0.0f, 0.0f },

		{ 0.0f, 0.0f, 0.0f },

		{ 0.0f, 0.0f, 0.0f } } };

	static probeWeights weights25 =

	    {

		    { // --- VERTICES (3) ---

		        { 1.000f, 0.000f, 0.000f },

		        { 0.000f, 1.000f, 0.000f },

		        { 0.000f, 0.000f, 1.000f },

		        // --- EDGES: 1/4 and 3/4 Points (6) ---

		        // These are critical for catching thin edges early

		        { 0.750f, 0.250f, 0.000f },

		        { 0.250f, 0.750f, 0.000f }, // Edge AB

		        { 0.750f, 0.000f, 0.250f },

		        { 0.250f, 0.000f, 0.750f }, // Edge AC

		        { 0.000f, 0.750f, 0.250f },

		        { 0.000f, 0.250f, 0.750f }, // Edge BC

		        // --- EDGES: 1/2 Mid-points (3) ---

		        { 0.500f, 0.500f, 0.000f },

		        { 0.500f, 0.000f, 0.500f },

		        { 0.000f, 0.500f, 0.500f },

		        // --- EDGES: 1/3 and 2/3 Points (6) ---

		        // More density to prevent lateral "pops"

		        { 0.666f, 0.333f, 0.000f },

		        { 0.333f, 0.666f, 0.000f },

		        { 0.666f, 0.000f, 0.333f },

		        { 0.333f, 0.000f, 0.666f },

		        { 0.000f, 0.666f, 0.333f },

		        { 0.000f, 0.333f, 0.666f },

		        // --- INTERIOR (7) ---

		        { 0.333f, 0.333f, 0.333f }, // Barycenter

		        { 0.450f, 0.450f, 0.100f },

		        { 0.450f, 0.100f, 0.450f },

		        { 0.100f, 0.450f, 0.450f },

		        { 0.200f, 0.200f, 0.600f },

		        { 0.200f, 0.600f, 0.200f },

		        { 0.600f, 0.200f, 0.200f } }

	    };

	struct SDFData
	{
		float distance;
		Eigen::Vector3f velocity;
	};

	struct IVec3
	{
		int x, y, z;
		bool operator==(const IVec3& o) const { return x == o.x && y == o.y && z == o.z; }
	};

	struct IVec3Hash
	{
		size_t operator()(const IVec3& v) const
		{
			return ((v.x * 73856093) ^ (v.y * 19349663) ^ (v.z * 83492791));
		}
	};

	class SDFVolume
	{
	private:
		std::unordered_map<IVec3, SDFData, IVec3Hash> sdfMap;
		float cellSize = 0.02f;
		int version = 1;

	public:
		SDFVolume() = default;
		explicit SDFVolume(const std::string& path)
		{
			if (!load(path))
			{
				std::cerr << "SDF Error: Could not load " << path << std::endl;
			}
		}

		bool exportToVTK(const std::string& path) const { throw "Not implemented!"; }
		void sample(const Eigen::Vector3d& p, float& outDist, Eigen::Vector3f& outVel) const
		{
			float invS = 1.0f / cellSize;
			Eigen::Vector3f pf = p.cast<float>() * invS;
			IVec3 base = { (int)std::floor(pf.x()), (int)std::floor(pf.y()), (int)std::floor(pf.z()) };
			Eigen::Vector3f d = pf - Eigen::Vector3f((float)base.x, (float)base.y, (float)base.z);

			float totalDist = 0.0f;
			Eigen::Vector3f totalVel = Eigen::Vector3f::Zero();
			float weightSum = 0.0f;

			for (int i = 0; i <= 1; ++i)
			{
				for (int j = 0; j <= 1; ++j)
				{
					for (int k = 0; k <= 1; ++k)
					{
						auto it = sdfMap.find({ base.x + i, base.y + j, base.z + k });
						if (it != sdfMap.end())
						{
							float w = (i ? d.x() : 1.0f - d.x()) * (j ? d.y() : 1.0f - d.y()) * (k ? d.z() : 1.0f - d.z());
							totalDist += it->second.distance * w;
							totalVel += it->second.velocity * w;
							weightSum += w;
						}
					}
				}
			}

			if (weightSum > 0.0001f)
			{
				outDist = totalDist / weightSum;
				outVel = totalVel / weightSum;
			}
			else
			{
				outDist = 999.0f;
				outVel = Eigen::Vector3f::Zero();
			}
		}

		double getDistance(const Eigen::Vector3d& p) const
		{
			float d;
			Eigen::Vector3f v;
			sample(p, d, v);
			return (double)d;
		}
		Eigen::Vector3f getVelocity(const Eigen::Vector3d& p) const
		{
			float d;
			Eigen::Vector3f v;
			sample(p, d, v);
			return v;
		}
		bool contains(const Eigen::Vector3d& p) const
		{
			IVec3 idx = { (int)std::floor(p.x() / cellSize), (int)std::floor(p.y() / cellSize), (int)std::floor(p.z() / cellSize) };
			return sdfMap.find(idx) != sdfMap.end();
		}

		bool load(const std::string& path)
		{
			std::ifstream in(path, std::ios::binary);
			if (!in.is_open())
				return false;
			sdfMap.clear();
			in.read((char*)&version, sizeof(int));
			if (version > 1000 || version < 1)
			{
				in.seekg(0);
				version = 1;
			}
			size_t count;
			in.read((char*)&count, sizeof(size_t));
			in.read((char*)&cellSize, sizeof(float));
			for (size_t i = 0; i < count; ++i)
			{
				IVec3 coord;
				SDFData data;
				in.read((char*)&coord, sizeof(IVec3));
				in.read((char*)&data.distance, sizeof(float));
				if (version == 2)
					in.read((char*)&data.velocity, sizeof(Eigen::Vector3f));
				else
					data.velocity = Eigen::Vector3f::Zero();
				sdfMap[coord] = data;
			}
			return true;
		}
	};

	inline void applyFriction(Real dt, Vector3& vel, const Real& kfriction, const Vector3& normal, const Vector3& v_body)
	{

		Vector3 v_rel = vel - v_body;
		// Decompose into normal and tangential components
		Real v_rel_norm = v_rel.dot(normal);
		Vector3 v_rel_tan = v_rel - v_rel_norm * normal;
		Real tan_norm = v_rel_tan.norm();

		Real friction_limit_vel = kfriction / dt;

		Vector3 v_rel_tan_new;
		if (tan_norm > friction_limit_vel)
			v_rel_tan_new = v_rel_tan.normalized() * (tan_norm - friction_limit_vel); // kinetic: reduce sliding
		else
			v_rel_tan_new = Vector3::Zero(); // static: stop sliding

		vel = v_body + v_rel_tan_new;
		if (v_rel_norm > 0)
			vel += v_rel_norm * normal;
	}

	void solveCollisionVertex(
	    const Obstacle* obs,
	    VectorN& X,
	    VectorN& V,
	    const VectorN& Xold,
	    Real dt,
	    Real threshold,
	    const Real& frictionK, unsigned int substep = 0,
	    std::vector<int>* contactVertices = NULL);

	void solveCollisionQuadrature(
	    const Obstacle* obs,
	    Real dt,
	    VectorN& X,
	    VectorN& V,
	    const VectorN& Xold,
	    const VectorN& W,
	    const std::vector<trimesh::triangle_t>& faces,
	    Real threshold,
	    const Real& frictionK,
	    const probeWeights& weights,
	    std::vector<int>* contactVertices,
	    int probe_points = 10);
}
#endif