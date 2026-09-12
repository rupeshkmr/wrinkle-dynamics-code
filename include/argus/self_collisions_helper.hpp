#ifndef SELF_COLLISIONS_HELPER
#define SELF_COLLISIONS_HELPER
#include "common.hpp"

namespace argus
{
	struct scParams
	{
		bool enable_self_collisions;
		Real self_collision_maxsubsteps;
		Real threshold;
		bool doccd;
		Real damping;
		Real stiffness;
	};
	void initializeStaticWeights(PrecomputedSelfCollisionData& staticData);
	void precomputeRestRadiiAndGrid(
	    const std::vector<trimesh::triangle_t>& faces,
	    const VectorN& Xrest,
	    PrecomputedSelfCollisionData& staticData,
	    Real safetyMultiplier = 0.55, // 0.5 makes them touch, 0.55 ensures overlap
	    Real minThickness = 0.0001); // Absolute floor (e.g., 1mm)
	// void initSelfCollision(const std::vector<trimesh::triangle_t>& faces, const VectorN& Xrest);
	void initSelfCollision(const std::vector<trimesh::triangle_t>& faces, const VectorN& Xrest, PrecomputedSelfCollisionData& selfCollisionData);

	void resolveSelfCollisionsSubstepped(
	    const std::vector<trimesh::triangle_t>& faces,
	    const VectorN& xt, // Position at the very start of the frame
	    VectorN& y, // Predicted position (modified to resolved position)
	    VectorN& v, // Velocity (modified)
	    int maxSubsteps,
	    const VectorN& invInertia,
	    PrecomputedSelfCollisionData& sd,
	    Real stiffness,
	    Real damping,
	    Real totalDt);
}

#endif