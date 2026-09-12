#ifndef SELF_COLLISIONS_HPP
#define SELF_COLLISIONS_HPP
#include "common.hpp"
namespace argus
{
	namespace selfcollisionhandler
	{
		extern probeWeights weights;
		struct CollisionResult
		{
			bool collided = false;
			double u, v, w; // Barycentric coordinates
			Vector3 normal; // Direction to push the vertex (from triangle to vertex)
			double distance; // Current distance between vertex and triangle
			double t; // Time of impact (from CCD)
		};
		struct VirtualPoint
		{
			Eigen::Vector3d pos; // Current world position of the probe
			long vIdx[3]; // Global indices of the 3 vertices forming the parent triangle
			double b[3]; // Barycentric weights (alpha, beta, gamma)
			int faceId; // Index of the triangle that owns this probe
		};
		struct SpatialHasher
		{
			Real cellSize;
			size_t tableSize;
			// l≈max(Average Edge Length,Max Displacement per Frame)
			// size should be 2 x elements
			// NEW BUFFERS: Persistent storage to avoid frame-by-frame allocation
			std::vector<int> bucketCounts; // Stores how many IDs per bucket
			std::vector<int> bucketOffsets; // Stores where each bucket starts in flatEntries
			std::vector<int> flatEntries; // One giant contiguous list of all IDs
			// NEW EDGE BUFFERS
			std::vector<int> edgeBucketCounts;
			std::vector<int> edgeBucketOffsets;
			std::vector<int> edgeFlatEntries;
			SpatialHasher(Real l, size_t size = 1000003)
			    : cellSize(l)
			    , tableSize(size)
			{
				bucketCounts.resize(tableSize);
				bucketOffsets.resize(tableSize + 1);
				edgeBucketCounts.resize(tableSize);
				edgeBucketOffsets.resize(tableSize + 1);
			}

			// Maps 3D grid coordinates (i, j, k) to a hash table index
			size_t hashCoords(int i, int j, int k) const
			{
				return (static_cast<size_t>(i * 73856093) ^ static_cast<size_t>(j * 19349663) ^ static_cast<size_t>(k * 83492791)) % tableSize;
			}

			// Converts a world position to grid indices
			void posToIndices(const Vector3& p, int& i, int& j, int& k) const
			{
				i = static_cast<int>(std::floor(p.x() / cellSize));
				j = static_cast<int>(std::floor(p.y() / cellSize));
				k = static_cast<int>(std::floor(p.z() / cellSize));
			}
			// Maps 3D grid coordinates (i, j, k) to a hash table index
			size_t hashCoords(const Vector3& p) const
			{
				int i = static_cast<int>(std::floor(p.x() / cellSize));
				int j = static_cast<int>(std::floor(p.y() / cellSize));
				int k = static_cast<int>(std::floor(p.z() / cellSize));
				return (static_cast<size_t>(i * 73856093) ^ static_cast<size_t>(j * 19349663) ^ static_cast<size_t>(k * 83492791)) % tableSize;
			}
		};

		struct VertexAABB
		{
			Vector3 min_p;
			Vector3 max_p;
			// NEW: Store grid indices to avoid re-calculation
			int i_min, j_min, k_min;
			int i_max, j_max, k_max;
			// Creates an AABB that covers the triangle's path from Xold to Xnew
			void expandToSweptLine(const Vector3& vold, const Vector3& vnew)
			{
				min_p = vold.cwiseMin(vnew);
				max_p = vold.cwiseMax(vnew);
			}
			// In VertexAABB
			void expandToSweptLine(const Vector3& vold, const Vector3& vnew, double thickness)
			{
				min_p = vold.cwiseMin(vnew) - Vector3::Constant(thickness);
				max_p = vold.cwiseMax(vnew) + Vector3::Constant(thickness);
			}
			const bool contains(const Vector3& pos) const
			{
				bool isinsidex = false;
				bool isinsidey = false;
				bool isinsidez = false;
				if (pos(0) <= max_p(0) && pos(0) >= min_p(0))
					isinsidex = true;
				if (pos(1) <= max_p(1) && pos(1) >= min_p(1))
					isinsidey = true;
				if (pos(2) <= max_p(2) && pos(2) >= min_p(2))
					isinsidez = true;
				return (isinsidex && isinsidey && isinsidez);
			}
		};
		struct AABB
		{
			Vector3 min_p;
			Vector3 max_p;
			int i_min, j_min, k_min;
			int i_max, j_max, k_max;
			// Creates an AABB that covers the triangle's path from Xold to Xnew
			void expandToSweptTriangle(const Vector3& v0_old, const Vector3& v1_old, const Vector3& v2_old,
			    const Vector3& v0_new, const Vector3& v1_new, const Vector3& v2_new)
			{
				min_p = v0_old.cwiseMin(v1_old).cwiseMin(v2_old).cwiseMin(v0_new).cwiseMin(v1_new).cwiseMin(v2_new);
				max_p = v0_old.cwiseMax(v1_old).cwiseMax(v2_new).cwiseMax(v0_new).cwiseMax(v1_new).cwiseMax(v2_new);
			}

			// In AABB (Triangle)
			void expandToSweptTriangle(const Vector3& v0_old, const Vector3& v1_old, const Vector3& v2_old,
			    const Vector3& v0_new, const Vector3& v1_new, const Vector3& v2_new, double thickness, const SpatialHasher& hasher)
			{
				min_p = v0_old.cwiseMin(v1_old).cwiseMin(v2_old).cwiseMin(v0_new).cwiseMin(v1_new).cwiseMin(v2_new) - Vector3::Constant(thickness);
				max_p = v0_old.cwiseMax(v1_old).cwiseMax(v2_old).cwiseMax(v0_new).cwiseMax(v1_new).cwiseMax(v2_new) + Vector3::Constant(thickness);
				hasher.posToIndices(min_p - Vector3::Constant(thickness), i_min, j_min, k_min);
				hasher.posToIndices(max_p + Vector3::Constant(thickness), i_max, j_max, k_max);
			}
			const bool contains(const Vector3& pos) const
			{
				bool isinsidex = false;
				bool isinsidey = false;
				bool isinsidez = false;
				if (pos(0) <= max_p(0) && pos(0) >= min_p(0))
					isinsidex = true;
				if (pos(1) <= max_p(1) && pos(1) >= min_p(1))
					isinsidey = true;
				if (pos(2) <= max_p(2) && pos(2) >= min_p(2))
					isinsidez = true;
				return (isinsidex && isinsidey && isinsidez);
			}
			const bool checkOverlap(const VertexAABB& vaabb) const
			{
				if (vaabb.max_p.x() < min_p.x() || vaabb.min_p.x() > max_p.x())
					return false;
				if (vaabb.max_p.y() < min_p.y() || vaabb.min_p.y() > max_p.y())
					return false;
				if (vaabb.max_p.z() < min_p.z() || vaabb.min_p.z() > max_p.z())
					return false;
				return true;
			}
		};
		struct EdgeAABB
		{
			Vector3 min_p;
			Vector3 max_p;
			int i_min, j_min, k_min;
			int i_max, j_max, k_max;

			// Creates an AABB covering the edge's movement from Xold to Xnew
			void expandToSweptEdge(const Vector3& e0_old, const Vector3& e1_old,
			    const Vector3& e0_new, const Vector3& e1_new)
			{
				// Find min and max across all 4 points
				min_p = e0_old.cwiseMin(e1_old).cwiseMin(e0_new).cwiseMin(e1_new);
				max_p = e0_old.cwiseMax(e1_old).cwiseMax(e0_new).cwiseMax(e1_new);
			}

			// In EdgeAABB
			void expandToSweptEdge(const Vector3& e0_old, const Vector3& e1_old,
			    const Vector3& e0_new, const Vector3& e1_new, double thickness, const SpatialHasher& hasher)
			{
				min_p = e0_old.cwiseMin(e1_old).cwiseMin(e0_new).cwiseMin(e1_new) - Vector3::Constant(thickness);
				max_p = e0_old.cwiseMax(e1_old).cwiseMax(e0_new).cwiseMax(e1_new) + Vector3::Constant(thickness);
				hasher.posToIndices(min_p - Vector3::Constant(thickness), i_min, j_min, k_min);
				hasher.posToIndices(max_p + Vector3::Constant(thickness), i_max, j_max, k_max);
			}
			const bool contains(const Vector3& pos) const
			{
				bool isinsidex = false;
				bool isinsidey = false;
				bool isinsidez = false;
				if (pos(0) <= max_p(0) && pos(0) >= min_p(0))
					isinsidex = true;
				if (pos(1) <= max_p(1) && pos(1) >= min_p(1))
					isinsidey = true;
				if (pos(2) <= max_p(2) && pos(2) >= min_p(2))
					isinsidez = true;
				return (isinsidex && isinsidey && isinsidez);
			}
			const bool checkOverlap(const EdgeAABB& eaabb) const
			{
				if (eaabb.max_p.x() < min_p.x() || eaabb.min_p.x() > max_p.x())
					return false;
				if (eaabb.max_p.y() < min_p.y() || eaabb.min_p.y() > max_p.y())
					return false;
				if (eaabb.max_p.z() < min_p.z() || eaabb.min_p.z() > max_p.z())
					return false;
				return true;
			}
		};

		void hashSweptTriangle(int tri_id,
		    const Vector3& v0_old, const Vector3& v1_old, const Vector3& v2_old,
		    const Vector3& v0_new, const Vector3& v1_new, const Vector3& v2_new,
		    const SpatialHasher& hasher, std::vector<std::vector<int>>& table, AABB& box);

		void hashSweptEdge(int edge_id,
		    const Vector3& e0_old, const Vector3& e1_old,
		    const Vector3& e0_new, const Vector3& e1_new,
		    const SpatialHasher& hasher, std::vector<std::vector<int>>& edgeTable, EdgeAABB& box);
		void hashSweptTriangle(int tri_id,
		    const Vector3& v0_old, const Vector3& v1_old, const Vector3& v2_old,
		    const Vector3& v0_new, const Vector3& v1_new, const Vector3& v2_new,
		    const SpatialHasher& hasher, std::vector<std::vector<int>>& table,
		    double thickness);
		void hashSweptEdge(int edge_id,
		    const Vector3& e0_old, const Vector3& e1_old,
		    const Vector3& e0_new, const Vector3& e1_new,
		    const SpatialHasher& hasher, std::vector<std::vector<int>>& edgeTable,
		    double thickness);

		void broadPhaseVertexFace(std::vector<int>& collisionPairs, const VectorN& X, const VectorN& Xold, const std::vector<trimesh::triangle_t>& triangles, const std::vector<AABB>& aabss, const SpatialHasher& hasher);
		void broadPhaseEdgeEdge(std::vector<int>& collisionPairs, const std::vector<trimesh::edge_t>& edges, const std::vector<EdgeAABB>& aabss, const SpatialHasher& hasher);
		void resolveSelfCollisionsVertexFace(VectorN& X, const VectorN& Xold, const VectorN& massInv, const std::vector<trimesh::triangle_t>& faces, const std::vector<trimesh::edge_t>& edges);
		void applyPBDCorrectionEdgeEdge(int va_id, int vb_id, int vc_id, int vd_id,
		    const CollisionResult& res,
		    const VectorN& invMass,
		    double thickness,
		    VectorN& deltaX, VectorN& weightX, bool repulsion = false);
		CollisionResult narrowPhaseEdgeEdge(const int& e1_v1, const int& e1_v2, const int& e2_v1, const int& e2_v2,
		    const Vector3& Pa_old, const Vector3& Pb_old,
		    const Vector3& Pc_old, const Vector3& Pd_old,
		    const Vector3& Pa_new, const Vector3& Pb_new,
		    const Vector3& Pc_new, const Vector3& Pd_new,
		    double thickness);
		CollisionResult runDiscreteEdgeEdgeTest(const Vector3& Pa, const Vector3& Pb,
		    const Vector3& Pc, const Vector3& Pd,
		    double thickness);
		CollisionResult narrowPhaseVertexFace(const int& vid, const int& fid, const Vector3& P_old, const Vector3& P0_old, const Vector3& P1_old, const Vector3& P2_old, const Vector3& P_new, const Vector3& P0_new, const Vector3& P1_new, const Vector3& P2_new);
		void applyPBDCorrectionVertexFace(int vid, int v0_id, int v1_id, int v2_id,
		    const CollisionResult& res,
		    const VectorN& invMass,
		    double thickness,
		    VectorN& deltaX, VectorN& weightX, bool repulsion = false);
		void solveAllSelfCollisions(
		    const std::vector<trimesh::triangle_t>& faces,
		    const std::vector<trimesh::edge_t>& edges,
		    const VectorN& Xold,
		    VectorN& Xnew,
		    const VectorN& massInv,
		    double thickness, bool doccd = true);
		// Main Solver
		void solveAllSelfCollisionsRepulsions(
		    const std::vector<trimesh::triangle_t>& faces,
		    const std::vector<trimesh::edge_t>& edges,
		    const VectorN& Xold,
		    VectorN& Xnew,
		    const VectorN& massInv,
		    double thickness,
		    bool doccd = true);

		// Point-Triangle distance (returns barycentric coords of closest point on triangle)
		double pointTriangleDistance(
		    const Vector3& p,
		    const Vector3& a, const Vector3& b, const Vector3& c,
		    double& u, double& v, double& w);

		// CCD via Bisection
		CollisionResult solveProbeTriangleCCD(
		    const Vector3& p_old, const Vector3& p_new,
		    const trimesh::triangle_t& tri,
		    const VectorN& Xold, const VectorN& Xnew,
		    double thickness);
		void solveAllSelfCollisionsRepulsions(
		    const std::vector<trimesh::triangle_t>& faces,
		    const VectorN& Xold,
		    VectorN& Xnew,
		    VectorN& V,
		    const VectorN& massInv,
		    PrecomputedSelfCollisionData& staticData, // h is now inside here
		    Real stiffness,
		    Real damping,
		    Real dt);
	}
}
#endif
