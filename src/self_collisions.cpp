#include "halfedge/trimesh_types.h"
#include <argus/self_collisions.hpp>
#include <omp.h>
#define MAX_COLLISION_BUFFER 128
namespace argus
{
	namespace selfcollisionhandler
	{
		probeWeights weights = probeWeights { { //
			// Vertices (3)
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
		CollisionResult runDiscreteProximityTest(const Vector3& P,
		    const Vector3& P0,
		    const Vector3& P1,
		    const Vector3& P2,
		    Real thickness)
		{
			CollisionResult result;
			result.collided = false; // Initialize to false

			// 1. Barycentric Projection Math (Your existing logic is perfect)
			Vector3 v0 = P1 - P0;
			Vector3 v1 = P2 - P0;
			Vector3 v2 = P - P0;

			Real d00 = v0.dot(v0);
			Real d01 = v0.dot(v1);
			Real d11 = v1.dot(v1);
			Real d20 = v2.dot(v0);
			Real d21 = v2.dot(v1);
			Real denom = d00 * d11 - d01 * d01;

			if (std::abs(denom) < 1e-12)
				return result;

			result.v = (d11 * d20 - d01 * d21) / denom;
			result.w = (d00 * d21 - d01 * d20) / denom;
			result.u = 1.0 - result.v - result.w;

			// 2. The "Inside Triangle" Check
			const Real bar_eps = 1e-4;
			if (result.u >= -bar_eps && result.v >= -bar_eps && result.w >= -bar_eps)
			{
				// 3. Distance calculation
				Vector3 closestPoint = result.u * P0 + result.v * P1 + result.w * P2;
				Vector3 diff = P - closestPoint;
				result.distance = diff.norm();

				if (result.distance < thickness)
				{
					result.collided = true;
					// Normal points from triangle TO vertex
					if (result.distance > 1e-12)
					{
						result.normal = diff / result.distance;
					}
					else
					{
						// If distance is zero, we use the geometric normal of the triangle
						Vector3 edge1 = P1 - P0;
						Vector3 edge2 = P2 - P0;
						Vector3 triNormal = edge1.cross(edge2);
						Real nLen = triNormal.norm();

						if (nLen > 1e-12)
						{
							result.normal = triNormal / nLen;
						}
						else
						{
							// Absolute fallback: just push up
							result.normal = Vector3(0, 1, 0);
						}
					}
				}
			}
			return result;
		}
		// Returns the earliest root t in [0, 1], or -1.0 if no root exists
		Real solveCubic(Real a, Real b, Real c, Real d)
		{
			const Real eps = 1e-11;

			// 1. If 'a' is near zero, it's actually a quadratic: bt^2 + ct + d = 0
			if (std::abs(a) < eps)
			{
				if (std::abs(b) < eps)
				{
					// Linear: ct + d = 0 -> t = -d/c
					if (std::abs(c) < eps)
						return -1.0;
					Real t = -d / c;
					return (t >= 0.0 && t <= 1.0) ? t : -1.0;
				}
				// Quadratic formula
				Real disc = c * c - 4.0 * b * d;
				if (disc < 0)
					return -1.0;
				Real sqrt_disc = std::sqrt(disc);
				Real t1 = (-c - sqrt_disc) / (2.0 * b);
				Real t2 = (-c + sqrt_disc) / (2.0 * b);

				Real best_t = 2.0; // placeholder
				if (t1 >= 0.0 && t1 <= 1.0)
					best_t = std::min(best_t, t1);
				if (t2 >= 0.0 && t2 <= 1.0)
					best_t = std::min(best_t, t2);
				return (best_t <= 1.0) ? best_t : -1.0;
			}

			// 2. Standard Cubic: Solve using Cardano's Method or Root Finding
			// For simplicity and robustness, many solvers use an analytical root finder
			// but here is the logic for finding the earliest root t in [0, 1]

			// Normalize to t^3 + At^2 + Bt + C = 0
			Real A = b / a;
			Real B = c / a;
			Real C = d / a;

			Real Q = (3.0 * B - A * A) / 9.0;
			Real R = (9.0 * A * B - 27.0 * C - 2.0 * A * A * A) / 54.0;
			Real D = Q * Q * Q + R * R; // Discriminant

			Real roots[3];
			int num_roots = 0;

			if (D > 0)
			{ // One real root
				Real S = std::cbrt(R + std::sqrt(D));
				Real T = std::cbrt(R - std::sqrt(D));
				roots[0] = -A / 3.0 + (S + T);
				num_roots = 1;
			}
			else if (D == 0)
			{ // All real, at least two equal
				Real S = std::cbrt(R);
				roots[0] = -A / 3.0 + 2.0 * S;
				roots[1] = -A / 3.0 - S;
				num_roots = 2;
			}
			else
			{ // Three distinct real roots
				Real theta = std::acos(R / std::sqrt(-Q * Q * Q));
				Real sqrt_neg_Q = std::sqrt(-Q);
				roots[0] = 2.0 * sqrt_neg_Q * std::cos(theta / 3.0) - A / 3.0;
				roots[1] = 2.0 * sqrt_neg_Q * std::cos((theta + 2.0 * M_PI) / 3.0) - A / 3.0;
				roots[2] = 2.0 * sqrt_neg_Q * std::cos((theta + 4.0 * M_PI) / 3.0) - A / 3.0;
				num_roots = 3;
			}

			// 3. Find the smallest root in [0, 1]
			Real min_t = 2.0;
			for (int i = 0; i < num_roots; ++i)
			{
				if (roots[i] >= 0.0 && roots[i] <= 1.0)
				{
					if (roots[i] < min_t)
						min_t = roots[i];
				}
			}

			return (min_t <= 1.0) ? min_t : -1.0;
		}
		Real checkCubicExtrema(Real a, Real b, Real c, Real d)
		{
			// f(t) = at^3 + bt^2 + ct + d
			// f'(t) = 3at^2 + 2bt + c
			Real da = 3.0 * a;
			Real db = 2.0 * b;
			Real dc = c;

			const Real eps = 1e-11;
			if (std::abs(da) < eps)
			{
				// Derivative is linear: 2bt + c = 0 -> t = -c / 2b
				if (std::abs(db) < eps)
					return -1.0;
				Real t = -dc / db;
				if (t > 0.0 && t < 1.0)
				{
					Real val = a * t * t * t + b * t * t + c * t + d;
					if (val * d < 0)
						return t;
				}
				return -1.0;
			}

			Real discriminant = db * db - 4.0 * da * dc;
			if (discriminant < 0)
				return -1.0;

			Real sqrt_d = std::sqrt(discriminant);
			Real roots[2];
			roots[0] = (-db - sqrt_d) / (2.0 * da);
			roots[1] = (-db + sqrt_d) / (2.0 * da);

			// We want the earliest t in (0, 1) that causes a sign change relative to d
			Real earliest_t = 2.0;
			for (int i = 0; i < 2; ++i)
			{
				Real t = roots[i];
				if (t > 0.0 && t < 1.0)
				{
					Real val = a * t * t * t + b * t * t + c * t + d;
					// Sign change check: if val has different sign than d (the start volume)
					if (val * d < 0)
					{
						if (t < earliest_t)
							earliest_t = t;
					}
				}
			}

			return (earliest_t < 1.0) ? earliest_t : -1.0;
		}
		void hashSweptEdge(int edge_id,
		    const Vector3& e0_old, const Vector3& e1_old,
		    const Vector3& e0_new, const Vector3& e1_new,
		    const SpatialHasher& hasher, std::vector<std::vector<int>>& edgeTable,
		    Real thickness)
		{
			EdgeAABB box;
			box.expandToSweptEdge(e0_old, e1_old, e0_new, e1_new);

			int i_min, j_min, k_min, i_max, j_max, k_max;
			hasher.posToIndices(box.min_p - Vector3::Constant(thickness), i_min, j_min, k_min);
			hasher.posToIndices(box.max_p + Vector3::Constant(thickness), i_max, j_max, k_max);

			for (int i = i_min; i <= i_max; ++i)
			{
				for (int j = j_min; j <= j_max; ++j)
				{
					for (int k = k_min; k <= k_max; ++k)
					{
						size_t h = hasher.hashCoords(i, j, k);

#pragma omp critical(edge_table_lock)
						{
							edgeTable[h].push_back(edge_id);
						}
					}
				}
			}
		}
		void hashSweptTriangle(int tri_id,
		    const Vector3& v0_old, const Vector3& v1_old, const Vector3& v2_old,
		    const Vector3& v0_new, const Vector3& v1_new, const Vector3& v2_new,
		    const SpatialHasher& hasher, std::vector<std::vector<int>>& table,
		    Real thickness)
		{
			// Create a temporary AABB to find the grid range
			AABB box;
			box.expandToSweptTriangle(v0_old, v1_old, v2_old, v0_new, v1_new, v2_new);

			// Inflate indices by thickness manually for the loop
			int i_min, j_min, k_min, i_max, j_max, k_max;
			hasher.posToIndices(box.min_p - Vector3::Constant(thickness), i_min, j_min, k_min);
			hasher.posToIndices(box.max_p + Vector3::Constant(thickness), i_max, j_max, k_max);

			for (int i = i_min; i <= i_max; ++i)
			{
				for (int j = j_min; j <= j_max; ++j)
				{
					for (int k = k_min; k <= k_max; ++k)
					{
						size_t h = hasher.hashCoords(i, j, k);

#pragma omp critical(face_table_lock)
						{
							table[h].push_back(tri_id);
						}
					}
				}
			}
		}
		void hashSweptTriangle(int tri_id,
		    const Vector3& v0_old, const Vector3& v1_old, const Vector3& v2_old,
		    const Vector3& v0_new, const Vector3& v1_new, const Vector3& v2_new,
		    const SpatialHasher& hasher, std::vector<std::vector<int>>& table, AABB& box)
		{
			box.expandToSweptTriangle(v0_old, v1_old, v2_old, v0_new, v1_new, v2_new);

			int i_min, j_min, k_min, i_max, j_max, k_max;
			hasher.posToIndices(box.min_p, i_min, j_min, k_min);
			hasher.posToIndices(box.max_p, i_max, j_max, k_max);

			// Iterate through all cells the swept volume might touch
			for (int i = i_min; i <= i_max; ++i)
			{
				for (int j = j_min; j <= j_max; ++j)
				{
					for (int k = k_min; k <= k_max; ++k)
					{
						size_t h = hasher.hashCoords(i, j, k);

						// Optimization: Ensure tri_id isn't added multiple times to same bucket
						// (Though with a large table, this is rare)
						table[h].push_back(tri_id);
					}
				}
			}
		}
		void hashSweptEdge(int edge_id,
		    const Vector3& e0_old, const Vector3& e1_old,
		    const Vector3& e0_new, const Vector3& e1_new,
		    const SpatialHasher& hasher, std::vector<std::vector<int>>& edgeTable, EdgeAABB& box)
		{

			box.expandToSweptEdge(e0_old, e1_old, e0_new, e1_new);

			int i_min, j_min, k_min, i_max, j_max, k_max;
			hasher.posToIndices(box.min_p, i_min, j_min, k_min);
			hasher.posToIndices(box.max_p, i_max, j_max, k_max);

			// Populate the specialized edge hash table
			for (int i = i_min; i <= i_max; ++i)
			{
				for (int j = j_min; j <= j_max; ++j)
				{
					for (int k = k_min; k <= k_max; ++k)
					{
						size_t h = hasher.hashCoords(i, j, k);
						edgeTable[h].push_back(edge_id);
					}
				}
			}
		}
		void broadPhaseVertexFace(std::vector<int>& collisionPairs, const VectorN& X, const VectorN& Xold, const std::vector<trimesh::triangle_t>& triangles, const std::vector<AABB>& aabss, const SpatialHasher& hasher)
		{
			int nverts = X.size() / 3;
			int nfaces = triangles.size();
			// possible pairs each vertex hits each face
			const int MAX_COLL = MAX_COLLISION_BUFFER;
			std::vector<int> counts(nverts, 0);

			std::vector<int> flatBuffer(nverts * MAX_COLL, -1); // Stores face IDs
			bool maxlimitexceed = false;
#ifdef USE_OMP
#pragma omp parallel for schedule(static)
#endif
			for (int vid = 0; vid < nverts; vid++)
			{
				// 1. Get positions from X_prev (start of step) and X (predicted)
				Vector3 p_old = Xold.segment<3>(3 * vid);
				Vector3 p_new = X.segment<3>(3 * vid);

				VertexAABB vaabb;
				vaabb.expandToSweptLine(p_old, p_new);

				// 2. Query ALL cells overlapping the vertex path AABB
				int i_min, j_min, k_min, i_max, j_max, k_max;
				hasher.posToIndices(vaabb.min_p, i_min, j_min, k_min);
				hasher.posToIndices(vaabb.max_p, i_max, j_max, k_max);

				int local_counter = 0;
				for (int i = i_min; i <= i_max; ++i)
				{
					for (int j = j_min; j <= j_max; ++j)
					{
						for (int k = k_min; k <= k_max; ++k)
						{
							size_t h = hasher.hashCoords(i, j, k);

							// 1. Get the contiguous memory range for this specific hash bucket
							int start = hasher.bucketOffsets[h];
							int end = hasher.bucketOffsets[h + 1];

							// 2. Iterate through the flat array instead of a vector of vectors
							for (int entryIdx = start; entryIdx < end; ++entryIdx)
							{
								// Retrieve the Face ID directly from the contiguous flatEntries array
								int fid = hasher.flatEntries[entryIdx];

								// 3. Adjacency Filter (Remains the same)
								if (triangles[fid].v[0] == vid || triangles[fid].v[1] == vid || triangles[fid].v[2] == vid)
									continue;

								// 4. Swept AABB check (Remains the same)
								if (aabss[fid].checkOverlap(vaabb))
								{
									if (local_counter < MAX_COLL)
									{
										// Check for duplicates
										bool already_added = false;
										for (int k_idx = 0; k_idx < local_counter; ++k_idx)
										{
											if (flatBuffer[vid * MAX_COLL + k_idx] == fid)
											{
												already_added = true;
												break;
											}
										}

										if (!already_added)
										{
											flatBuffer[vid * MAX_COLL + local_counter] = fid;
											local_counter++;
										}
									}
									else
									{
										maxlimitexceed = true;
									}
								}
							}
						}
					}
				}
				counts[vid] = local_counter;
			}
			collisionPairs.swap(flatBuffer);
			if (maxlimitexceed)
				std::cerr << "warning maximum self-collision face-vertex pairs exceeded!" << std::endl;
		}
		/* void broadPhaseEdgeEdge(std::vector<int>& collisionPairs, const std::vector<trimesh::edge_t>& edges, const std::vector<EdgeAABB>& aabss, const SpatialHasher& hasher)
		{
		    int nedges = edges.size();
		    const int MAX_COLL_EE = MAX_COLLISION_BUFFER; // Edges can sometimes have more candidates
		    std::vector<int> eeCounts(nedges, 0);
		    std::vector<int> eeFlatBuffer(nedges * MAX_COLL_EE, -1);

#ifdef USE_OMP
#pragma omp parallel for schedule(static)
#endif
		    for (int ei_id = 0; ei_id < nedges; ei_id++)
		    {
		        const auto& edge_i = edges[ei_id];
		        const EdgeAABB& box_i = aabss[ei_id]; // Swept AABB of edge i

		        // 1. Get the range of cells this edge overlaps
		        int i_min, j_min, k_min, i_max, j_max, k_max;
		        hasher.posToIndices(box_i.min_p, i_min, j_min, k_min);
		        hasher.posToIndices(box_i.max_p, i_max, j_max, k_max);

		        int counter = 0;
		        // 2. Query all buckets in the AABB range
		        for (int i = i_min; i <= i_max; ++i)
		        {
		            for (int j = j_min; j <= j_max; ++j)
		            {
		                for (int k = k_min; k <= k_max; ++k)
		                {
		                    size_t h = hasher.hashCoords(i, j, k);

		                    // 3. Check against edges stored in this bucket
		                    for (int ej_id : table[h])
		                    {

		                        // A. Avoid checking the edge against itself
		                        // B. Use ID comparison (ei < ej) to avoid checking the same pair twice
		                        if (ei_id >= ej_id)
		                            continue;

		                        // C. Neighbor Filter: Skip if edges share a vertex
		                        if (edge_i.v[0] == edges[ej_id].v[0] || edge_i.v[0] == edges[ej_id].v[1] || edge_i.v[1] == edges[ej_id].v[0] || edge_i.v[1] == edges[ej_id].v[1])
		                        {
		                            continue;
		                        }

		                        // D. Mid-Phase: Swept AABB overlap
		                        if (box_i.checkOverlap(aabss[ej_id]))
		                        {
		                            if (counter < MAX_COLL_EE)
		                            {
		                                // Check for duplicates (an edge can be in multiple shared buckets)
		                                bool already_added = false;
		                                for (int k_idx = 0; k_idx < counter; ++k_idx)
		                                {
		                                    if (eeFlatBuffer[ei_id * MAX_COLL_EE + k_idx] == ej_id)
		                                    {
		                                        already_added = true;
		                                        break;
		                                    }
		                                }

		                                if (!already_added)
		                                {
		                                    eeFlatBuffer[ei_id * MAX_COLL_EE + counter] = ej_id;
		                                    counter++;
		                                }
		                            }
		                            else
		                            {
		                                std::cerr << "Error in edge edge broadphase collision check! Max collisions reached!";
		                            }
		                        }
		                    }
		                }
		            }
		        }
		        eeCounts[ei_id] = counter;
		    }
		    collisionPairs.swap(eeFlatBuffer);
		} */
		void broadPhaseEdgeEdge(std::vector<int>& collisionPairs,
		    const std::vector<trimesh::edge_t>& edges,
		    const std::vector<EdgeAABB>& aabss,
		    const SpatialHasher& hasher) // <-- CHANGE: Pass hasher directly
		{
			int nedges = edges.size();
			const int MAX_COLL_EE = MAX_COLLISION_BUFFER;
			std::vector<int> eeFlatBuffer(nedges * MAX_COLL_EE, -1);
			bool maxlimitexceed = false;

#ifdef USE_OMP
#pragma omp parallel for schedule(static)
#endif
			for (int ei_id = 0; ei_id < nedges; ei_id++)
			{
				const auto& edge_i = edges[ei_id];
				const EdgeAABB& box_i = aabss[ei_id];

				int i_min, j_min, k_min, i_max, j_max, k_max;
				hasher.posToIndices(box_i.min_p, i_min, j_min, k_min);
				hasher.posToIndices(box_i.max_p, i_max, j_max, k_max);

				int counter = 0;
				for (int i = i_min; i <= i_max; ++i)
				{
					for (int j = j_min; j <= j_max; ++j)
					{
						for (int k = k_min; k <= k_max; ++k)
						{
							size_t h = hasher.hashCoords(i, j, k);

							// >>> CHANGE START: Replacing "for (int ej_id : table[h])" <<<
							int start = hasher.edgeBucketOffsets[h];
							int end = hasher.edgeBucketOffsets[h + 1];

							for (int entryIdx = start; entryIdx < end; ++entryIdx)
							{
								int ej_id = hasher.edgeFlatEntries[entryIdx];
								// >>> CHANGE END <<<

								// A. Avoid checking the edge against itself & double counting
								if (ei_id >= ej_id)
									continue;

								// C. Neighbor Filter: Skip if edges share a vertex
								if (edge_i.v[0] == edges[ej_id].v[0] || edge_i.v[0] == edges[ej_id].v[1] || edge_i.v[1] == edges[ej_id].v[0] || edge_i.v[1] == edges[ej_id].v[1])
								{
									continue;
								}

								// D. Mid-Phase: Swept AABB overlap
								if (box_i.checkOverlap(aabss[ej_id]))
								{
									if (counter < MAX_COLL_EE)
									{
										bool already_added = false;
										for (int k_idx = 0; k_idx < counter; ++k_idx)
										{
											if (eeFlatBuffer[ei_id * MAX_COLL_EE + k_idx] == ej_id)
											{
												already_added = true;
												break;
											}
										}

										if (!already_added)
										{
											eeFlatBuffer[ei_id * MAX_COLL_EE + counter] = ej_id;
											counter++;
										}
									}
									else
									{
										maxlimitexceed = true;
									}
								}
							}
						}
					}
				}
			}
			collisionPairs.swap(eeFlatBuffer);
			if (maxlimitexceed)
				std::cerr << "warning maximum self-collision edge pairs exceeded!" << std::endl;
		}
		void applyPBDCorrectionVertexFace(int vid, int v0_id, int v1_id, int v2_id,
		    const CollisionResult& res,
		    const VectorN& invMass,
		    Real thickness,
		    VectorN& deltaX, VectorN& weightX, bool repulsion)
		{
			if (!res.collided)
				return;

			Real C = res.distance - thickness;
			if (C >= 0.0)
				return;

			Real w_v = invMass[vid];
			Real w0 = invMass[v0_id], w1 = invMass[v1_id], w2 = invMass[v2_id];

			// Effective triangle mass at the impact point
			Real w_tri = (w0 * res.u * res.u) + (w1 * res.v * res.v) + (w2 * res.w * res.w);
			Real sum_w = w_v + w_tri;

			if (sum_w < 1e-12)
				return;

			Real lambda = C / sum_w;

			/* // Apply displacements to the predicted positions X
			X.segment<3>(3 * vid) -= w_v * lambda * res.normal;
			X.segment<3>(3 * v0_id) += w0 * res.u * lambda * res.normal;
			X.segment<3>(3 * v1_id) += w1 * res.v * lambda * res.normal;
			X.segment<3>(3 * v2_id) += w2 * res.w * lambda * res.normal; */
			Vector3 d_v = -w_v * lambda * res.normal;
			Vector3 d_v0 = w0 * res.u * lambda * res.normal;
			Vector3 d_v1 = w1 * res.v * lambda * res.normal;
			Vector3 d_v2 = w2 * res.w * lambda * res.normal;

			// Use ATOMIC to add to the global accumulation buffers.
			// This is much faster than 'critical' because it uses hardware-level locks.
			for (int i = 0; i < 3; ++i)
			{
#pragma omp atomic
				deltaX[3 * vid + i] += d_v[i];
#pragma omp atomic
				deltaX[3 * v0_id + i] += d_v0[i];
#pragma omp atomic
				deltaX[3 * v1_id + i] += d_v1[i];
#pragma omp atomic
				deltaX[3 * v2_id + i] += d_v2[i];
			}
			Real ww = 1;
			if (repulsion)
				ww = 2;
// Track how many constraints are affecting these vertices
#pragma omp atomic
			weightX[vid] += ww;
#pragma omp atomic
			weightX[v0_id] += ww;
#pragma omp atomic
			weightX[v1_id] += ww;
#pragma omp atomic
			weightX[v2_id] += ww;
		}
		CollisionResult narrowPhaseVertexFace(const int& vid, const int& fid, const Vector3& P_old, const Vector3& P0_old, const Vector3& P1_old, const Vector3& P2_old, const Vector3& P_new, const Vector3& P0_new, const Vector3& P1_new, const Vector3& P2_new)
		{
			// 1. Relative positions at t=0
			Vector3 u1 = P1_old - P0_old;
			Vector3 u2 = P2_old - P0_old;
			Vector3 u3 = P_old - P0_old;

			// 2. Relative velocities (differences in displacement)
			Vector3 v_p = P_new - P_old;
			Vector3 v_p0 = P0_new - P0_old;
			Vector3 v_p1 = P1_new - P1_old;
			Vector3 v_p2 = P2_new - P2_old;

			Vector3 w1 = v_p1 - v_p0;
			Vector3 w2 = v_p2 - v_p0;
			Vector3 w3 = v_p - v_p0;

			// 3. Calculate coefficients
			// a = w3 . (w1 x w2)
			Real a = w3.dot(w1.cross(w2));

			// b = u3 . (w1 x w2) + w3 . (u1 x w2 + w1 x u2)
			Real b = u3.dot(w1.cross(w2)) + w3.dot(u1.cross(w2) + w1.cross(u2));

			// c = w3 . (u1 x u2) + u3 . (u1 x w2 + w1 x u2)
			Real c = w3.dot(u1.cross(u2)) + u3.dot(u1.cross(w2) + w1.cross(u2));

			// d = u3 . (u1 x u2)
			Real d = u3.dot(u1.cross(u2));

			const Real eps = 1e-10;

			// Case 1: Signed volumes at start and end
			Real vol_start = d;
			Real vol_end = a + b + c + d;

			// Check if they are already coplanar or stay on the same side
			bool sign_change = (vol_start * vol_end <= 0.0);
			CollisionResult res;
			// Case 2: Degeneracy Check (Moving in the same plane)
			// If the relative velocity vectors are almost coplanar with the triangle
			bool is_degenerate = (std::abs(a) < eps && std::abs(b) < eps && std::abs(c) < eps);
			if (is_degenerate)
			{
				// FALLBACK: If they are moving in the same plane, the cubic solver
				// is unstable. Use a discrete proximity test for the current frame.
				// TODO: use sheet thickness here
				res = runDiscreteProximityTest(P_new, P0_new, P1_new, P2_new, 1e-4);

				// if (res.collided)
				// {
				// 	// Now you have all the data needed for PBD:
				// 	// res.u, res.v, res.w, res.normal, res.distance
				// 	// applyPBDCorrection(vid, fid, res);
				// }
			}
			else if (sign_change)
			{
				// IMPACT LIKELY: The vertex crossed the plane of the triangle.
				// Proceed to the cubic root finder.
				Real t = solveCubic(a, b, c, d);

				// If analytical solver fails but signs differ, Bisection is GUARANTEED to find the root.
				if (t < 0.0 && (vol_start * vol_end < 0))
				{
					Real t0 = 0.0, t1 = 1.0;
					for (int i = 0; i < 32; ++i)
					{ // 32 iterations gives near-Real precision
						Real mid = (t0 + t1) * 0.5;
						Real f_mid = a * mid * mid * mid + b * mid * mid + c * mid + d;
						if (f_mid * (a * t0 * t0 * t0 + b * t0 * t0 + c * t0 + d) < 0)
							t1 = mid;
						else
							t0 = mid;
					}
					t = (t0 + t1) * 0.5;
				}
				if (t >= 0.0 && t <= 1.0)
				{
					// Interpolate positions to time t
					Vector3 Pt = P_old + t * (P_new - P_old);
					Vector3 P0t = P0_old + t * (P0_new - P0_old);
					Vector3 P1t = P1_old + t * (P1_new - P1_old);
					Vector3 P2t = P2_old + t * (P2_new - P2_old);

					// Call your proximity test at time t (with a very small thickness/epsilon)
					res = runDiscreteProximityTest(Pt, P0t, P1t, P2t, 1e-6);

					// if (res.collided)
					// {
					// 	// Impact confirmed at time t!
					// 	// Resolve using the original P_new and the normal calculated at time t
					// 	// applyPBDCorrection(vid, fid, res);
					// }
				}
			}
			else
			{
				// EXTREME VALUES CHECK: Even if there is no sign change,
				// the vertex could have crossed and come back (rare but possible).
				// Usually, we check the local extrema of the cubic in [0, 1].
				Real t = checkCubicExtrema(a, b, c, d);

				if (t >= 0.0 && t <= 1.0)
				{
					// Interpolate positions to time t
					Vector3 Pt = P_old + t * (P_new - P_old);
					Vector3 P0t = P0_old + t * (P0_new - P0_old);
					Vector3 P1t = P1_old + t * (P1_new - P1_old);
					Vector3 P2t = P2_old + t * (P2_new - P2_old);

					// Call your proximity test at time t (with a very small thickness/epsilon)
					res = runDiscreteProximityTest(Pt, P0t, P1t, P2t, 1e-6);

					// if (res.collided)
					// {
					// 	// Impact confirmed at time t!
					// 	// Resolve using the original P_new and the normal calculated at time t
					// 	// applyPBDCorrection(vid, fid, res);
					// }
				}
			}
			return res;
		}

		CollisionResult runDiscreteEdgeEdgeTest(const Vector3& Pa, const Vector3& Pb,
		    const Vector3& Pc, const Vector3& Pd,
		    Real thickness)
		{
			CollisionResult res;
			res.collided = false;

			Vector3 u = Pb - Pa;
			Vector3 v = Pd - Pc;
			Vector3 w = Pa - Pc;
			Real a = u.dot(u), b = u.dot(v), c = v.dot(v), d = u.dot(w), e = v.dot(w);
			Real D = a * c - b * b;
			Real sc, tc;

			if (D < 1e-12)
			{ // parallel
				sc = 0.0;
				tc = (b > c ? d / b : e / c);
			}
			else
			{
				sc = (b * e - c * d) / D;
				tc = (a * e - b * d) / D;
			}

			// Clamp parameters to [0, 1] segment range
			sc = std::max(0.0, std::min(1.0, sc));
			tc = std::max(0.0, std::min(1.0, tc));

			Vector3 Q1 = Pa + sc * u;
			Vector3 Q2 = Pc + tc * v;
			Vector3 diff = Q1 - Q2;
			res.distance = diff.norm();

			if (res.distance < thickness)
			{
				res.collided = true;
				res.u = sc; // Param for Edge 1
				res.v = tc; // Param for Edge 2
				res.normal = (res.distance > 1e-12) ? diff / res.distance : u.cross(v).normalized();
			}
			return res;
		}
		CollisionResult narrowPhaseEdgeEdge(const int& e1_v1, const int& e1_v2, const int& e2_v1, const int& e2_v2,
		    const Vector3& Pa_old, const Vector3& Pb_old,
		    const Vector3& Pc_old, const Vector3& Pd_old,
		    const Vector3& Pa_new, const Vector3& Pb_new,
		    const Vector3& Pc_new, const Vector3& Pd_new,
		    Real thickness)
		{
			// Relative vectors at t=0
			Vector3 u1 = Pb_old - Pa_old;
			Vector3 u2 = Pd_old - Pc_old;
			Vector3 u3 = Pc_old - Pa_old;

			// Relative velocities
			Vector3 va = Pa_new - Pa_old;
			Vector3 vb = Pb_new - Pb_old;
			Vector3 vc = Pc_new - Pc_old;
			Vector3 vd = Pd_new - Pd_old;

			Vector3 w1 = vb - va;
			Vector3 w2 = vd - vc;
			Vector3 w3 = vc - va;

			// Cubic coefficients for EE coplanarity
			Real a = w1.dot(w2.cross(w3));
			Real b = u1.dot(w2.cross(w3)) + w1.dot(u2.cross(w3) + w2.cross(u3));
			Real c = w1.dot(u2.cross(u3)) + u1.dot(u2.cross(w3) + w2.cross(u3));
			Real d = u1.dot(u2.cross(u3));

			Real vol_start = d;
			Real vol_end = a + b + c + d;
			bool sign_change = (vol_start * vol_end <= 0.0);

			Real t_hit = -1.0;
			if (sign_change)
			{
				t_hit = solveCubic(a, b, c, d);
				// If analytical solver fails but signs differ, Bisection is GUARANTEED to find the root.
				if (t_hit < 0.0 && (vol_start * vol_end < 0))
				{
					Real t0 = 0.0, t1 = 1.0;
					for (int i = 0; i < 32; ++i)
					{ // 32 iterations gives near-Real precision
						Real mid = (t0 + t1) * 0.5;
						Real f_mid = a * mid * mid * mid + b * mid * mid + c * mid + d;
						if (f_mid * (a * t0 * t0 * t0 + b * t0 * t0 + c * t0 + d) < 0)
							t1 = mid;
						else
							t0 = mid;
					}
					t_hit = (t0 + t1) * 0.5;
				}
			}
			else
			{
				t_hit = checkCubicExtrema(a, b, c, d);
			}

			CollisionResult res;
			if (t_hit >= 0.0 && t_hit <= 1.0)
			{
				// Interpolate to t_hit
				Vector3 Pat = Pa_old + t_hit * (Pa_new - Pa_old);
				Vector3 Pbt = Pb_old + t_hit * (Pb_new - Pb_old);
				Vector3 Pct = Pc_old + t_hit * (Pc_new - Pc_old);
				Vector3 Pdt = Pd_old + t_hit * (Pd_new - Pd_old);

				// Run segment-segment distance test at t_hit
				res = runDiscreteEdgeEdgeTest(Pat, Pbt, Pct, Pdt, thickness);
				res.t = t_hit;
			}
			return res;
		}
		void applyPBDCorrectionEdgeEdge(int va_id, int vb_id, int vc_id, int vd_id,
		    const CollisionResult& res,
		    const VectorN& invMass,
		    Real thickness,
		    VectorN& deltaX, VectorN& weightX, bool repulsion)
		{
			if (!res.collided)
				return;

			Real C = res.distance - thickness;
			if (C >= 0.0)
				return;

			Real wa = invMass[va_id], wb = invMass[vb_id];
			Real wc = invMass[vc_id], wd = invMass[vd_id];

			// s and t from the test are stored in res.u and res.v
			Real s = res.u;
			Real t = res.v;

			// Effective mass for EE constraint
			Real w1 = wa * (1.0 - s) * (1.0 - s) + wb * s * s;
			Real w2 = wc * (1.0 - t) * (1.0 - t) + wd * t * t;
			Real sum_w = w1 + w2;

			if (sum_w < 1e-12)
				return;

			Real lambda = C / sum_w;
			// Calculate displacement vectors for all 4 vertices
			Vector3 d_va = -wa * (1.0 - s) * lambda * res.normal;
			Vector3 d_vb = -wb * s * lambda * res.normal;
			Vector3 d_vc = wc * (1.0 - t) * lambda * res.normal;
			Vector3 d_vd = wd * t * lambda * res.normal;

			// Atomic updates for positions (3 components each)
			for (int i = 0; i < 3; ++i)
			{
#pragma omp atomic
				deltaX[3 * va_id + i] += d_va[i];
#pragma omp atomic
				deltaX[3 * vb_id + i] += d_vb[i];
#pragma omp atomic
				deltaX[3 * vc_id + i] += d_vc[i];
#pragma omp atomic
				deltaX[3 * vd_id + i] += d_vd[i];
			}
			Real ww = 1;
			if (repulsion)
				ww = 2;
// Increment weights for all 4 vertices
#pragma omp atomic
			weightX[va_id]
			    += ww;
#pragma omp atomic
			weightX[vb_id] += ww;
#pragma omp atomic
			weightX[vc_id] += ww;
#pragma omp atomic
			weightX[vd_id] += ww;
			/* // Apply displacements
			X.segment<3>(3 * va_id) -= wa * (1.0 - s) * lambda * res.normal;
			X.segment<3>(3 * vb_id) -= wb * s * lambda * res.normal;

			X.segment<3>(3 * vc_id) += wc * (1.0 - t) * lambda * res.normal;
			X.segment<3>(3 * vd_id) += wd * t * lambda * res.normal; */
		}
		void solveAllSelfCollisions(
		    const std::vector<trimesh::triangle_t>& faces,
		    const std::vector<trimesh::edge_t>& edges,
		    const VectorN& Xold,
		    VectorN& Xnew,
		    const VectorN& massInv,
		    Real thickness, bool doccd)
		{
			const int nverts = Xnew.size() / 3;
			const int nfaces = faces.size();
			const int nedges = edges.size();
			const int MAX_COLL = MAX_COLLISION_BUFFER;

			// Cell size should be large enough to contain the movement + thickness
			// to avoid an object spanning too many cells.
			SpatialHasher hasher(thickness * 4.0, 1000003 /* nfaces * 2 */);
			std::fill(hasher.bucketCounts.begin(), hasher.bucketCounts.end(), 0);

			// 1. Create AABB Buffers
			std::vector<AABB> faceAABBs(nfaces);
			std::vector<EdgeAABB> edgeAABBs(nedges);

// 2. Hash and Expand AABBs (Parallel)
// Note: hashSwept functions must contain internal #pragma omp critical blocks
#ifdef USE_OMP
#pragma omp parallel for
#endif
			for (int i = 0; i < nfaces; ++i)
			{
				const auto& f = faces[i];
				Vector3 v0o = Xold.segment<3>(3 * f.v[0]), v1o = Xold.segment<3>(3 * f.v[1]), v2o = Xold.segment<3>(3 * f.v[2]);
				Vector3 v0n = Xnew.segment<3>(3 * f.v[0]), v1n = Xnew.segment<3>(3 * f.v[1]), v2n = Xnew.segment<3>(3 * f.v[2]);

				faceAABBs[i].expandToSweptTriangle(v0o, v1o, v2o, v0n, v1n, v2n, thickness, hasher);

				// // Updated to pass thickness

				for (int x = faceAABBs[i].i_min; x <= faceAABBs[i].i_max; ++x)
				{
					for (int y = faceAABBs[i].j_min; y <= faceAABBs[i].j_max; ++y)
					{
						for (int z = faceAABBs[i].k_min; z <= faceAABBs[i].k_max; ++z)
						{
							size_t h = hasher.hashCoords(x, y, z);
// ATOMIC: No lock required, very fast!
#pragma omp atomic
							hasher.bucketCounts[h]++;
						}
					}
				}
			}
			hasher.bucketOffsets[0] = 0;
			for (int i = 0; i < hasher.tableSize; ++i)
			{
				hasher.bucketOffsets[i + 1] = hasher.bucketOffsets[i] + hasher.bucketCounts[i];
			}

			// Resize the flat entries to hold every single instance found in the counts
			hasher.flatEntries.resize(hasher.bucketOffsets[hasher.tableSize]);

			// Create a temporary copy of offsets to use as 'write pointers'
			std::vector<int> writePtrs = hasher.bucketOffsets;
#ifdef USE_OMP
#pragma omp parallel for
#endif
			for (int i = 0; i < nfaces; ++i)
			{

				for (int x = faceAABBs[i].i_min; x <= faceAABBs[i].i_max; ++x)
				{
					for (int y = faceAABBs[i].j_min; y <= faceAABBs[i].j_max; ++y)
					{
						for (int z = faceAABBs[i].k_min; z <= faceAABBs[i].k_max; ++z)

						{
							size_t h = hasher.hashCoords(x, y, z);

							int pos;
// ATOMIC CAPTURE: Gets a unique slot for this thread
#pragma omp atomic capture
							pos = writePtrs[h]++;

							hasher.flatEntries[pos] = i; // Lock-free write!
						}
					}
				}
			}

			// Reset edge counts
			std::fill(hasher.edgeBucketCounts.begin(), hasher.edgeBucketCounts.end(), 0);
#ifdef USE_OMP
#pragma omp parallel for
#endif
			for (int i = 0; i < nedges; ++i)
			{
				const auto& e = edges[i];
				Vector3 e0o = Xold.segment<3>(3 * e.v[0]), e1o = Xold.segment<3>(3 * e.v[1]);
				Vector3 e0n = Xnew.segment<3>(3 * e.v[0]), e1n = Xnew.segment<3>(3 * e.v[1]);

				edgeAABBs[i].expandToSweptEdge(e0o, e1o, e0n, e1n, thickness, hasher);

				// Updated to pass thickness

				for (int x = edgeAABBs[i].i_min; x <= edgeAABBs[i].i_max; ++x)
				{
					for (int y = edgeAABBs[i].j_min; y <= edgeAABBs[i].j_max; ++y)
					{
						for (int z = edgeAABBs[i].k_min; z <= edgeAABBs[i].k_max; ++z)

						{
							size_t h = hasher.hashCoords(x, y, z);
#pragma omp atomic
							hasher.edgeBucketCounts[h]++;
						}
					}
				}
			}
			hasher.edgeBucketOffsets[0] = 0;
			for (int i = 0; i < hasher.tableSize; ++i)
			{
				hasher.edgeBucketOffsets[i + 1] = hasher.edgeBucketOffsets[i] + hasher.edgeBucketCounts[i];
			}

			// Resize the flat entries to hold every edge instance
			hasher.edgeFlatEntries.resize(hasher.edgeBucketOffsets[hasher.tableSize]);

			// Create write pointers for edges
			std::vector<int> edgeWritePtrs = hasher.edgeBucketOffsets;
#ifdef USE_OMP
#pragma omp parallel for
#endif
			for (int i = 0; i < nedges; ++i)
			{

				for (int x = edgeAABBs[i].i_min; x <= edgeAABBs[i].i_max; ++x)
				{
					for (int y = edgeAABBs[i].j_min; y <= edgeAABBs[i].j_max; ++y)
					{
						for (int z = edgeAABBs[i].k_min; z <= edgeAABBs[i].k_max; ++z)

						{
							size_t h = hasher.hashCoords(x, y, z);

							int pos;
#pragma omp atomic capture
							pos = edgeWritePtrs[h]++;

							hasher.edgeFlatEntries[pos] = i;
						}
					}
				}
			}

			// 3. Broadphase Candidate Detection

			// new
			std::vector<int> vfPairs;
			std::vector<int> eePairs;

			// Pass the hasher directly instead of the faceTable/edgeTable
			broadPhaseVertexFace(vfPairs, Xnew, Xold, faces, faceAABBs, hasher);
			broadPhaseEdgeEdge(eePairs, edges, edgeAABBs, hasher);

			// 4. Vertex-Face Narrow Phase & Correction
			// Initialize to zero every frame
			VectorN deltaX = VectorN::Zero(Xnew.size());
			VectorN weightX = VectorN::Zero(nverts);

			if (doccd)
			{ // Global ccd->proximity check
#ifdef USE_OMP
#pragma omp parallel for schedule(dynamic)
#endif
				for (int vid = 0; vid < nverts; vid++)
				{
					const Vector3 p_o = Xold.segment<3>(3 * vid);
					// CHANGE: Use a local constant start position.
					// We don't update Xnew here, so we don't "refresh" p_n anymore.
					const Vector3 p_n_start = Xnew.segment<3>(3 * vid);

					for (int c = 0; c < MAX_COLL; ++c)
					{
						int fid = vfPairs[vid * MAX_COLL + c];
						if (fid == -1)
							break;

						const auto& tri = faces[fid];
						if (vid == tri.v[0] || vid == tri.v[1] || vid == tri.v[2])
							continue;

						// CACHE: Triangle positions
						const Vector3 t0o = Xold.segment<3>(3 * tri.v[0]), t1o = Xold.segment<3>(3 * tri.v[1]), t2o = Xold.segment<3>(3 * tri.v[2]);
						const Vector3 t0n = Xnew.segment<3>(3 * tri.v[0]), t1n = Xnew.segment<3>(3 * tri.v[1]), t2n = Xnew.segment<3>(3 * tri.v[2]);

						// 2. CCD
						CollisionResult res = narrowPhaseVertexFace(vid, fid, p_o, t0o, t1o, t2o, p_n_start, t0n, t1n, t2n);

						if (res.collided)
						{
							// CHANGE: Call the accumulation version
							applyPBDCorrectionVertexFace(vid, tri.v[0], tri.v[1], tri.v[2],
							    res, massInv, thickness, deltaX, weightX);
						}

						// 3. Proximity Fallback
						CollisionResult resProx = runDiscreteProximityTest(p_n_start, t0n, t1n, t2n, thickness);
						if (resProx.collided)
						{
							// CHANGE: Call the accumulation version
							applyPBDCorrectionVertexFace(vid, tri.v[0], tri.v[1], tri.v[2],
							    resProx, massInv, thickness, deltaX, weightX);
						}
					}
				}
				// 5. Edge-Edge Narrow Phase & Correction
#ifdef USE_OMP
#pragma omp parallel for schedule(dynamic)
#endif
				for (int eid = 0; eid < nedges; eid++)
				{
					const auto& e1 = edges[eid];
					// CACHE: Edge 1 positions once per outer loop
					const Vector3 e1v0o = Xold.segment<3>(3 * e1.v[0]), e1v1o = Xold.segment<3>(3 * e1.v[1]);
					Vector3 e1v0n = Xnew.segment<3>(3 * e1.v[0]), e1v1n = Xnew.segment<3>(3 * e1.v[1]);

					for (int c = 0; c < MAX_COLL; ++c)
					{
						int ejid = eePairs[eid * MAX_COLL + c];
						if (ejid == -1)
							break;

						const auto& e2 = edges[ejid];
						if (eid == ejid)
							continue;

						// Neighbor-check
						if (e1.v[0] == e2.v[0] || e1.v[0] == e2.v[1] || e1.v[1] == e2.v[0] || e1.v[1] == e2.v[1])
							continue;

						// >>> CHANGE: CACHE Edge 2 positions locally to avoid repeated .segment<3> calls <<<
						const Vector3 e2v0o = Xold.segment<3>(3 * e2.v[0]), e2v1o = Xold.segment<3>(3 * e2.v[1]);
						Vector3 e2v0n = Xnew.segment<3>(3 * e2.v[0]), e2v1n = Xnew.segment<3>(3 * e2.v[1]);

						// 2. CCD Logic
						// >>> CHANGE: Pass cached old and new positions to the CCD solver <<<
						CollisionResult res = narrowPhaseEdgeEdge(e1.v[0], e1.v[1], e2.v[0], e2.v[1],
						    e1v0o, e1v1o, e2v0o, e2v1o,
						    e1v0n, e1v1n, e2v0n, e2v1n,
						    thickness);

						if (res.collided)
						{

							applyPBDCorrectionEdgeEdge(e1.v[0], e1.v[1], e2.v[0], e2.v[1], res, massInv, thickness, deltaX, weightX);
						}

						// 3. Proximity Fallback
						// >>> CHANGE: Use cached positions for final safety check <<<
						CollisionResult resProx = runDiscreteEdgeEdgeTest(e1v0n, e1v1n, e2v0n, e2v1n, thickness);

						if (resProx.collided)
						{

							applyPBDCorrectionEdgeEdge(e1.v[0], e1.v[1], e2.v[0], e2.v[1], resProx, massInv, thickness, deltaX, weightX);
						}
					}
				}
#ifdef USE_OMP
#pragma omp parallel for
#endif
				for (int i = 0; i < nverts; i++)
				{
					if (weightX[i] > 0.0)
					{
						Xnew.segment<3>(3 * i) += (deltaX.segment<3>(3 * i) / weightX[i]) * 1.5;
						deltaX.segment<3>(3 * i).setZero();
						weightX[i] = 0;
					}
				}
			}
			// Global repulsion
#ifdef USE_OMP
#pragma omp parallel for schedule(dynamic)
#endif
			for (int vid = 0; vid < nverts; vid++)
			{
				const Vector3 p_o = Xold.segment<3>(3 * vid);
				// CHANGE: Use a local constant start position.
				// We don't update Xnew here, so we don't "refresh" p_n anymore.
				const Vector3 p_n_start = Xnew.segment<3>(3 * vid);

				for (int c = 0; c < MAX_COLL; ++c)
				{
					int fid = vfPairs[vid * MAX_COLL + c];
					if (fid == -1)
						break;

					const auto& tri = faces[fid];
					if (vid == tri.v[0] || vid == tri.v[1] || vid == tri.v[2])
						continue;

					// CACHE: Triangle positions
					const Vector3 t0o = Xold.segment<3>(3 * tri.v[0]), t1o = Xold.segment<3>(3 * tri.v[1]), t2o = Xold.segment<3>(3 * tri.v[2]);
					const Vector3 t0n = Xnew.segment<3>(3 * tri.v[0]), t1n = Xnew.segment<3>(3 * tri.v[1]), t2n = Xnew.segment<3>(3 * tri.v[2]);

					// 1. Repulsion
					Real repulsionDist = thickness * 2.0;
					CollisionResult resRepel = runDiscreteProximityTest(p_n_start, t0n, t1n, t2n, repulsionDist);

					if (resRepel.collided)
					{
						// CHANGE: Call the accumulation version of the correction
						applyPBDCorrectionVertexFace(vid, tri.v[0], tri.v[1], tri.v[2],
						    resRepel, massInv, repulsionDist, deltaX, weightX, true);
					}
				}
			}
#ifdef USE_OMP
#pragma omp parallel for schedule(dynamic)
#endif
			for (int eid = 0; eid < nedges; eid++)
			{
				const auto& e1 = edges[eid];
				// CACHE: Edge 1 positions once per outer loop
				const Vector3 e1v0o = Xold.segment<3>(3 * e1.v[0]), e1v1o = Xold.segment<3>(3 * e1.v[1]);
				Vector3 e1v0n = Xnew.segment<3>(3 * e1.v[0]), e1v1n = Xnew.segment<3>(3 * e1.v[1]);

				for (int c = 0; c < MAX_COLL; ++c)
				{
					int ejid = eePairs[eid * MAX_COLL + c];
					if (ejid == -1)
						break;

					const auto& e2 = edges[ejid];
					if (eid == ejid)
						continue;

					// Neighbor-check
					if (e1.v[0] == e2.v[0] || e1.v[0] == e2.v[1] || e1.v[1] == e2.v[0] || e1.v[1] == e2.v[1])
						continue;

					// >>> CHANGE: CACHE Edge 2 positions locally to avoid repeated .segment<3> calls <<<
					const Vector3 e2v0o = Xold.segment<3>(3 * e2.v[0]), e2v1o = Xold.segment<3>(3 * e2.v[1]);
					Vector3 e2v0n = Xnew.segment<3>(3 * e2.v[0]), e2v1n = Xnew.segment<3>(3 * e2.v[1]);

					// 1. Repulsion Logic
					Real repulsionDist = thickness * 2.0;
					// >>> CHANGE: Use cached e1 and e2 positions <<<
					CollisionResult resRepel = runDiscreteEdgeEdgeTest(e1v0n, e1v1n, e2v0n, e2v1n, repulsionDist);

					if (resRepel.collided)
					{

						applyPBDCorrectionEdgeEdge(e1.v[0], e1.v[1], e2.v[0], e2.v[1], resRepel, massInv, repulsionDist, deltaX, weightX, true);
					}
				}
			}
#ifdef USE_OMP
#pragma omp parallel for
#endif
			for (int i = 0; i < nverts; i++)
			{
				if (weightX[i] > 0.0)
				{
					Xnew.segment<3>(3 * i) += (deltaX.segment<3>(3 * i) / weightX[i]) * 1.5;
					deltaX.segment<3>(3 * i).setZero();
					weightX[i] = 0;
				}
			}
		}
		void solveAllSelfCollisionsRepulsions(
		    const std::vector<trimesh::triangle_t>& faces,
		    const std::vector<trimesh::edge_t>& edges,
		    const VectorN& Xold,
		    VectorN& Xnew,
		    const VectorN& massInv,
		    double thickness, bool doccd)
		{
			SpatialHasher hasher(thickness * 4.0);
			std::vector<std::vector<int>> triTable(hasher.tableSize);

			// 1. Hash Triangles
			for (int i = 0; i < (int)faces.size(); ++i)
			{
				hashSweptTriangle(i,
				    Xold.segment<3>(faces[i].v[0] * 3), Xold.segment<3>(faces[i].v[1] * 3), Xold.segment<3>(faces[i].v[2] * 3),
				    Xnew.segment<3>(faces[i].v[0] * 3), Xnew.segment<3>(faces[i].v[1] * 3), Xnew.segment<3>(faces[i].v[2] * 3),
				    hasher, triTable, thickness);
			}

			// 2. Resolve Probes
			for (int i = 0; i < (int)faces.size(); ++i)
			{
				const auto& f_p = faces[i];
				for (int w_idx = 0; w_idx < 25; ++w_idx)
				{
					Real wu = weights.weights[w_idx][0], wv = weights.weights[w_idx][1], ww = weights.weights[w_idx][2];
					if (wu + wv + ww < 1e-6)
						continue;

					Vector3 p_old = wu * Xold.segment<3>(f_p.v[0] * 3) + wv * Xold.segment<3>(f_p.v[1] * 3) + ww * Xold.segment<3>(f_p.v[2] * 3);
					Vector3 p_new = wu * Xnew.segment<3>(f_p.v[0] * 3) + wv * Xnew.segment<3>(f_p.v[1] * 3) + ww * Xnew.segment<3>(f_p.v[2] * 3);

					size_t cell = hasher.hashCoords(p_new);
					for (int target_id : triTable[cell])
					{
						// Topology Filter
						bool skip = false;
						for (int a = 0; a < 3; ++a)
							for (int b = 0; b < 3; ++b)
								if (f_p.v[a] == faces[target_id].v[b])
									skip = true;
						if (skip)
							continue;

						CollisionResult res;
						if (doccd)
						{
							double low = 0.0, high = 1.0;
							for (int s = 0; s < 8; ++s)
							{
								double mid = (low + high) * 0.5;
								Vector3 pm = (1. - mid) * p_old + mid * p_new;
								Vector3 am = (1. - mid) * Xold.segment<3>(faces[target_id].v[0] * 3) + mid * Xnew.segment<3>(faces[target_id].v[0] * 3);
								Vector3 bm = (1. - mid) * Xold.segment<3>(faces[target_id].v[1] * 3) + mid * Xnew.segment<3>(faces[target_id].v[1] * 3);
								Vector3 cm = (1. - mid) * Xold.segment<3>(faces[target_id].v[2] * 3) + mid * Xnew.segment<3>(faces[target_id].v[2] * 3);
								double d = pointTriangleDistance(pm, am, bm, cm, res.u, res.v, res.w);
								if (d < thickness)
								{
									high = mid;
									res.collided = true;
									res.distance = d;
									res.t = mid;
									res.normal = (bm - am).cross(cm - am).normalized();
									if ((pm - am).dot(res.normal) < 0)
										res.normal *= -1.0;
								}
								else
									low = mid;
							}
						}

						if (res.collided)
						{
							Vector3 impulse = res.normal * (thickness - res.distance) * 0.5;
							for (int k = 0; k < 3; ++k)
								Xnew.segment<3>(f_p.v[k] * 3) += weights.weights[w_idx][k] * massInv(f_p.v[k] * 3) * impulse;
							double tw[3] = { res.u, res.v, res.w };
							for (int k = 0; k < 3; ++k)
								Xnew.segment<3>(faces[target_id].v[k] * 3) -= tw[k] * massInv(faces[target_id].v[k] * 3) * impulse;
						}
					}
				}
			}
		}

		// Robust Point-Triangle Distance Implementation
		double pointTriangleDistance(const Vector3& p, const Vector3& a, const Vector3& b, const Vector3& c, double& u, double& v, double& w)
		{
			Vector3 ab = b - a, ac = c - a, ap = p - a;
			double d1 = ab.dot(ap), d2 = ac.dot(ap);
			if (d1 <= 0.0 && d2 <= 0.0)
			{
				u = 1;
				v = 0;
				w = 0;
				return (p - a).norm();
			}
			Vector3 bp = p - b;
			double d3 = ab.dot(bp), d4 = ac.dot(bp);
			if (d3 >= 0.0 && d4 <= d3)
			{
				u = 0;
				v = 1;
				w = 0;
				return (p - b).norm();
			}
			double vc = d1 * d4 - d3 * d2;
			if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
			{
				v = d1 / (d1 - d3);
				u = 1.0 - v;
				w = 0;
				return (p - (a + v * ab)).norm();
			}
			Vector3 cp = p - c;
			double d5 = ab.dot(cp), d6 = ac.dot(cp);
			if (d6 >= 0.0 && d5 <= d6)
			{
				u = 0;
				v = 0;
				w = 1;
				return (p - c).norm();
			}
			double vb = d5 * d2 - d1 * d6;
			if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
			{
				w = d2 / (d2 - d6);
				u = 1.0 - w;
				v = 0;
				return (p - (a + w * ac)).norm();
			}
			double va = d3 * d6 - d5 * d4;
			if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
			{
				w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
				v = 1.0 - w;
				u = 0;
				return (p - (b + w * (c - b))).norm();
			}
			double denom = 1.0 / (va + vb + vc);
			v = vb * denom;
			w = vc * denom;
			u = 1.0 - v - w;
			return (p - (a + v * ab + w * ac)).norm();
		}

		// 		void solveAllSelfCollisionsRepulsions(
		// 		    const std::vector<trimesh::triangle_t>& faces,
		// 		    const VectorN& Xold,
		// 		    VectorN& Xnew,
		// 		    VectorN& V,
		// 		    const VectorN& massInv,
		// 		    PrecomputedSelfCollisionData& staticData,
		// 		    Real stiffness,
		// 		    Real damping,
		// 		    Real dt)
		// 		{
		// 			const int nVerts = Xnew.rows() / 3;
		// 			const int numPoints = staticData.pointCloudBuffer.size();
		// 			const Real h = staticData.gridCellSize;
		// 			const Real invH = 1.0 / h;

		// // --- PHASE 1: GATHER ---
		// #pragma omp parallel for
		// 			for (int i = 0; i < numPoints; ++i)
		// 			{
		// 				auto& pt = staticData.pointCloudBuffer[i];
		// 				pt.pos = pt.b[0] * Xnew.segment<3>(3 * pt.vIdx[0]) + pt.b[1] * Xnew.segment<3>(3 * pt.vIdx[1]) + pt.b[2] * Xnew.segment<3>(3 * pt.vIdx[2]);
		// 			}

		// 			// --- PHASE 2: SPATIAL HASH ---
		// 			const int tableSize = numPoints;
		// 			std::vector<int> head(tableSize, -1);
		// 			std::vector<int> next(numPoints, -1);
		// 			auto hashFunc = [&](const Vector3& pos)
		// 			{
		// 				long x = (long)std::floor(pos.x() * invH);
		// 				long y = (long)std::floor(pos.y() * invH);
		// 				long z = (long)std::floor(pos.z() * invH);
		// 				return static_cast<size_t>((x * 73856093) ^ (y * 19349663) ^ (z * 83492791)) % tableSize;
		// 			};
		// 			for (int i = 0; i < numPoints; ++i)
		// 			{
		// 				size_t hIdx = hashFunc(staticData.pointCloudBuffer[i].pos);
		// 				next[i] = head[hIdx];
		// 				head[hIdx] = i;
		// 			}

		// 			// --- PHASE 3: RESOLUTION ---
		// 			int numThreads = omp_get_max_threads();
		// 			std::vector<std::vector<Real>> threadDelta(numThreads, std::vector<Real>(nVerts * 3, 0.0));
		// 			Real effective_stiffness = 1.0 - std::exp(-stiffness * dt);

		// #pragma omp parallel
		// 			{
		// 				int tid = omp_get_thread_num();
		// 				Real* localDelta = threadDelta[tid].data();
		// #pragma omp for
		// 				for (int i = 0; i < numPoints; ++i)
		// 				{
		// 					const auto& pi = staticData.pointCloudBuffer[i];
		// 					long cx = (long)std::floor(pi.pos.x() * invH);
		// 					long cy = (long)std::floor(pi.pos.y() * invH);
		// 					long cz = (long)std::floor(pi.pos.z() * invH);

		// 					for (long dx = -1; dx <= 1; ++dx)
		// 					{
		// 						for (long dy = -1; dy <= 1; ++dy)
		// 						{
		// 							for (long dz = -1; dz <= 1; ++dz)
		// 							{
		// 								size_t hIdx = static_cast<size_t>(((cx + dx) * 73856093) ^ ((cy + dy) * 19349663) ^ ((cz + dz) * 83492791)) % tableSize;
		// 								int j = head[hIdx];
		// 								while (j != -1)
		// 								{
		// 									if (i < j)
		// 									{
		// 										const auto& pj = staticData.pointCloudBuffer[j];
		// 										if (pi.faceId != pj.faceId)
		// 										{
		// 											// Adjacency Check
		// 											bool shared = false;
		// 											for (int a = 0; a < 3; ++a)
		// 												for (int b = 0; b < 3; ++b)
		// 													if (pi.vIdx[a] == pj.vIdx[b])
		// 														shared = true;

		// 											if (!shared)
		// 											{
		// 												Vector3 diff = pi.pos - pj.pos;
		// 												Real distSq = diff.squaredNorm();
		// 												Real minGap = pi.radius + pj.radius;

		// 												if (distSq < minGap * minGap && distSq > 1e-12)
		// 												{
		// 													Real dist = std::sqrt(distSq);
		// 													Vector3 normal = diff / dist;
		// 													// CLAMP overlap to 2% of cell size to prevent huge jumps
		// 													Real overlap = std::min(minGap - dist, h * 0.02);

		// 													Real w1 = 0, w2 = 0;
		// 													for (int k = 0; k < 3; ++k)
		// 													{
		// 														w1 += pi.b[k] * pi.b[k] * massInv(3 * pi.vIdx[k]);
		// 														w2 += pj.b[k] * pj.b[k] * massInv(3 * pj.vIdx[k]);
		// 													}

		// 													if (w1 + w2 > 1e-12)
		// 													{
		// 														// Stiffness Impulse
		// 														// Real p = (effective_stiffness * overlap) / (w1 + w2);
		// 														Real p = (stiffness * overlap) / (w1 + w2);
		// 														Vector3 impulse = p * normal;

		// 														// FRICTION BLOCK (Integrated & Stable)
		// 														const Real mu = 0.001;
		// 														Vector3 vi = pi.b[0] * V.segment<3>(3 * pi.vIdx[0]) + pi.b[1] * V.segment<3>(3 * pi.vIdx[1]) + pi.b[2] * V.segment<3>(3 * pi.vIdx[2]);
		// 														Vector3 vj = pj.b[0] * V.segment<3>(3 * pj.vIdx[0]) + pj.b[1] * V.segment<3>(3 * pj.vIdx[1]) + pj.b[2] * V.segment<3>(3 * pj.vIdx[2]);
		// 														Vector3 v_rel = vi - vj;
		// 														Vector3 v_tangent = v_rel - v_rel.dot(normal) * normal;
		// 														if (v_tangent.norm() > 1e-8)
		// 														{
		// 															// Real tangent_disp = v_tangent.norm() * dt;
		// 															// impulse -= std::min(mu * p, tangent_disp / (w1 + w2)) * v_tangent.normalized();
		// 															impulse -= std::min(mu * p, v_tangent.norm() * 0.1 / (w1 + w2)) * v_tangent.normalized();
		// 														}

		// 														for (int k = 0; k < 3; ++k)
		// 														{
		// 															Vector3 n1 = (massInv(3 * pi.vIdx[k]) * pi.b[k]) * impulse;
		// 															Vector3 n2 = (massInv(3 * pj.vIdx[k]) * pj.b[k]) * impulse;
		// 															localDelta[3 * pi.vIdx[k] + 0] += n1.x();
		// 															localDelta[3 * pi.vIdx[k] + 1] += n1.y();
		// 															localDelta[3 * pi.vIdx[k] + 2] += n1.z();
		// 															localDelta[3 * pj.vIdx[k] + 0] -= n2.x();
		// 															localDelta[3 * pj.vIdx[k] + 1] -= n2.y();
		// 															localDelta[3 * pj.vIdx[k] + 2] -= n2.z();
		// 														}
		// 													}
		// 												}
		// 											}
		// 										}
		// 									}
		// 									j = next[j];
		// 								}
		// 							}
		// 						}
		// 					}
		// 				}
		// 			}
		// 			Real dampingFactor = std::exp(-damping * dt);
		// // --- PHASE 4: ACCUMULATE ---
		// #pragma omp parallel for
		// 			for (int v = 0; v < nVerts; ++v)
		// 			{
		// 				Vector3 d(0, 0, 0);
		// 				for (int t = 0; t < numThreads; ++t)
		// 				{
		// 					d.x() += threadDelta[t][3 * v];
		// 					d.y() += threadDelta[t][3 * v + 1];
		// 					d.z() += threadDelta[t][3 * v + 2];
		// 				}

		// 				// Hard clamp of movement to 10% of cell size
		// 				if (d.norm() > h * 0.1)
		// 					d = d.normalized() * (h * 0.1);

		// 				Xnew.segment<3>(3 * v) += d;
		// 				// Real dfactor = 0.95; // original
		// 				// Real dfactor = 1;
		// 				// VELOCITY FILTER: Apply a damping factor (0.95) to the update
		// 				if (dt > 1e-8)
		// 				{
		// 					Vector3 newV = (Xnew.segment<3>(3 * v) - Xold.segment<3>(3 * v)) / dt;
		// 					V.segment<3>(3 * v) = newV * dampingFactor; //(1 - damping);
		// 				}
		// 			}
		// 		}
		// 	}
		void solveAllSelfCollisionsRepulsions(
		    const std::vector<trimesh::triangle_t>& faces,
		    const VectorN& Xold,
		    VectorN& Xnew,
		    VectorN& V,
		    const VectorN& massInv,
		    PrecomputedSelfCollisionData& staticData,
		    Real stiffness,
		    Real damping,
		    Real dt)
		{
			const int nVerts = Xnew.rows() / 3;
			const int numPoints = staticData.pointCloudBuffer.size();
			const Real h = staticData.gridCellSize;
			const Real invH = 1.0 / h;

// --- PHASE 1: GATHER ---
#pragma omp parallel for
			for (int i = 0; i < numPoints; ++i)
			{
				auto& pt = staticData.pointCloudBuffer[i];
				pt.pos = pt.b[0] * Xnew.segment<3>(3 * pt.vIdx[0]) + pt.b[1] * Xnew.segment<3>(3 * pt.vIdx[1]) + pt.b[2] * Xnew.segment<3>(3 * pt.vIdx[2]);
			}

			// --- PHASE 2: SPATIAL HASH ---
			const int tableSize = numPoints;
			std::vector<int> head(tableSize, -1);
			std::vector<int> next(numPoints, -1);
			auto hashFunc = [&](const Vector3& pos)
			{
				long x = (long)std::floor(pos.x() * invH);
				long y = (long)std::floor(pos.y() * invH);
				long z = (long)std::floor(pos.z() * invH);
				return static_cast<size_t>((x * 73856093) ^ (y * 19349663) ^ (z * 83492791)) % tableSize;
			};
			for (int i = 0; i < numPoints; ++i)
			{
				size_t hIdx = hashFunc(staticData.pointCloudBuffer[i].pos);
				next[i] = head[hIdx];
				head[hIdx] = i;
			}

			// --- PHASE 3: RESOLUTION ---
			int numThreads = omp_get_max_threads();
			std::vector<std::vector<Real>> threadDelta(numThreads, std::vector<Real>(nVerts * 3, 0.0));
			Real effective_stiffness = 1.0 - std::exp(-stiffness * dt);
#pragma omp parallel for
			for (int i = 0; i < numPoints; ++i)
			{
				auto& p = staticData.pointCloudBuffer[i];
				// Precalculate the weighted inverse mass for this point
				p.wInv = p.b[0] * p.b[0] * massInv(3 * p.vIdx[0]) + p.b[1] * p.b[1] * massInv(3 * p.vIdx[1]) + p.b[2] * p.b[2] * massInv(3 * p.vIdx[2]);
			}
			const Real* Vdata = V.data();
#pragma omp parallel
			{
				int tid = omp_get_thread_num();
				Real* localDelta = threadDelta[tid].data();
#pragma omp for
				for (int i = 0; i < numPoints; ++i)
				{
					const auto& pi = staticData.pointCloudBuffer[i];
					long cx = (long)std::floor(pi.pos.x() * invH);
					long cy = (long)std::floor(pi.pos.y() * invH);
					long cz = (long)std::floor(pi.pos.z() * invH);

					for (long dx = -1; dx <= 1; ++dx)
					{
						for (long dy = -1; dy <= 1; ++dy)
						{
							for (long dz = -1; dz <= 1; ++dz)
							{
								size_t hIdx = static_cast<size_t>(((cx + dx) * 73856093) ^ ((cy + dy) * 19349663) ^ ((cz + dz) * 83492791)) % tableSize;
								int j = head[hIdx];
								while (j != -1)
								{
									if (i < j)
									{
										const auto& pj = staticData.pointCloudBuffer[j];
										if (pi.faceId != pj.faceId)
										{
											// Adjacency Check
											// bool shared = false;
											// for (int a = 0; a < 3; ++a)
											// 	for (int b = 0; b < 3; ++b)
											// 		if (pi.vIdx[a] == pj.vIdx[b])
											// 			shared = true;
											bool shared = (pi.vIdx[0] == pj.vIdx[0] || pi.vIdx[0] == pj.vIdx[1] || pi.vIdx[0] == pj.vIdx[2] || pi.vIdx[1] == pj.vIdx[0] || pi.vIdx[1] == pj.vIdx[1] || pi.vIdx[1] == pj.vIdx[2] || pi.vIdx[2] == pj.vIdx[0] || pi.vIdx[2] == pj.vIdx[1] || pi.vIdx[2] == pj.vIdx[2]);

											if (!shared)
											{
												Vector3 diff = pi.pos - pj.pos;
												Real distSq = diff.squaredNorm();
												Real minGap = pi.radius + pj.radius;

												if (distSq < minGap * minGap && distSq > 1e-12)
												{
													Real dist = std::sqrt(distSq);
													Vector3 normal = diff / dist;
													// CLAMP overlap to 2% of cell size to prevent huge jumps
													Real overlap = std::min(minGap - dist, h * 0.02);

													// Real w1 = 0, w2 = 0;
													// for (int k = 0; k < 3; ++k)
													// {
													// 	w1 += pi.b[k] * pi.b[k] * massInv(3 * pi.vIdx[k]);
													// 	w2 += pj.b[k] * pj.b[k] * massInv(3 * pj.vIdx[k]);
													// }

													// if (w1 + w2 > 1e-12)
													Real wSum = pi.wInv + pj.wInv;
													if (wSum > 1e-12)
													{
														// Stiffness Impulse
														// Real p = (effective_stiffness * overlap) / (w1 + w2);
														Real p = (stiffness * overlap) / (wSum);
														Vector3 impulse = p * normal;

														// FRICTION BLOCK (Integrated & Stable)
														const Real mu = 0.001;
														// Vector3 vi = pi.b[0] * V.segment<3>(3 * pi.vIdx[0]) + pi.b[1] * V.segment<3>(3 * pi.vIdx[1]) + pi.b[2] * V.segment<3>(3 * pi.vIdx[2]);
														// Vector3 vj = pj.b[0] * V.segment<3>(3 * pj.vIdx[0]) + pj.b[1] * V.segment<3>(3 * pj.vIdx[1]) + pj.b[2] * V.segment<3>(3 * pj.vIdx[2]);
														Vector3 vi = pi.b[0] * Vector3::Map(Vdata + 3 * pi.vIdx[0]) + pi.b[1] * Vector3::Map(Vdata + 3 * pi.vIdx[1]) + pi.b[2] * Vector3::Map(Vdata + 3 * pi.vIdx[2]);
														Vector3 vj = pj.b[0] * Vector3::Map(Vdata + 3 * pj.vIdx[0]) + pj.b[1] * Vector3::Map(Vdata + 3 * pj.vIdx[1]) + pj.b[2] * Vector3::Map(Vdata + 3 * pj.vIdx[2]);
														Vector3 v_rel = vi - vj;
														Vector3 v_tangent = v_rel - v_rel.dot(normal) * normal;
														if (v_tangent.norm() > 1e-8)
														{
															// Real tangent_disp = v_tangent.norm() * dt;
															// impulse -= std::min(mu * p, tangent_disp / (w1 + w2)) * v_tangent.normalized();
															impulse -= std::min(mu * p, v_tangent.norm() * 0.1 / (wSum)) * v_tangent.normalized();
														}

														for (int k = 0; k < 3; ++k)
														{
															Vector3 n1 = (massInv(3 * pi.vIdx[k]) * pi.b[k]) * impulse;
															Vector3 n2 = (massInv(3 * pj.vIdx[k]) * pj.b[k]) * impulse;
															localDelta[3 * pi.vIdx[k] + 0] += n1.x();
															localDelta[3 * pi.vIdx[k] + 1] += n1.y();
															localDelta[3 * pi.vIdx[k] + 2] += n1.z();
															localDelta[3 * pj.vIdx[k] + 0] -= n2.x();
															localDelta[3 * pj.vIdx[k] + 1] -= n2.y();
															localDelta[3 * pj.vIdx[k] + 2] -= n2.z();
														}
													}
												}
											}
										}
									}
									j = next[j];
								}
							}
						}
					}
				}
			}
			Real dampingFactor = std::exp(-damping * dt);
// --- PHASE 4: ACCUMULATE ---
#pragma omp parallel for
			for (int v = 0; v < nVerts; ++v)
			{
				Vector3 d(0, 0, 0);
				for (int t = 0; t < numThreads; ++t)
				{
					d.x() += threadDelta[t][3 * v];
					d.y() += threadDelta[t][3 * v + 1];
					d.z() += threadDelta[t][3 * v + 2];
				}

				// Hard clamp of movement to 10% of cell size
				if (d.norm() > h * 0.1)
					d = d.normalized() * (h * 0.1);

				Xnew.segment<3>(3 * v) += d;
				// Real dfactor = 0.95; // original
				// Real dfactor = 1;
				// VELOCITY FILTER: Apply a damping factor (0.95) to the update
				if (dt > 1e-8)
				{
					Vector3 newV = (Xnew.segment<3>(3 * v) - Xold.segment<3>(3 * v)) / dt;
					V.segment<3>(3 * v) = newV * dampingFactor; //(1 - damping);
				}
			}
		}
	}

}
