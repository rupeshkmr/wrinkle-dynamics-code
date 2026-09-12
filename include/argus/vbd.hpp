#ifndef VBD_HPP
#define VBD_HPP

#include "common.hpp"
#include "time_integrator.hpp"
#include "obstacles.hpp"
#include "wrinkles_vbd.hpp"
#include "optimizers.hpp"
namespace argus
{

	class VBD : public TimeIntegrator
	{
		SharedVector m_B; // provides a container of VectorN
		Real m_kdamp_pos;
		bool m_sim_friction;
		unsigned int m_pos_vbd_iterations;
		PrecomputedSelfCollisionData selfCollisionData;
		nlohmann::json m_stats, m_config;
		std::string m_sdfCollider;

		void initSelfCollision(const std::vector<trimesh::triangle_t>& faces, const VectorN& Xrest);
		struct
		{
			bool enable_self_collisions;
			Real self_collision_maxsubsteps;
			Real threshold;
			bool doccd;
			Real stiffness;
			Real damping;
		} selfCollisionParams;

	public:
		void init(std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs, nlohmann::json& params) override;
		void initMeshes(std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obstacles);
		void resolveSelfCollisionsSubstepped(
		    const std::vector<trimesh::triangle_t>& faces,
		    const VectorN& xt, // Position at start of frame
		    VectorN& y, // Target position
		    VectorN& v, // Frame velocity
		    int maxSubsteps,
		    const VectorN& invInertia,
		    PrecomputedSelfCollisionData& sd,
		    Real stiffness,
		    Real damping,
		    Real totalDt);
		nlohmann::json getStats() override { return m_stats; }

		void step(Real timestep, std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs, nlohmann::json& state, unsigned int frame = 0) override;
		void initVBD(std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects);
		void vbdUpdate(Real timestep, std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs);
		void vbdUpdatePos(Real timestep, std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs, nlohmann::json& state, unsigned int frame);
		VectorN getAdaptiveAccelerationForVBD(DynamicsObject& object, VectorN& var, bool pos);
		void resolveCollisions(std::vector<std::unique_ptr<Obstacle>>& obs, VectorN& updated_pos, const VectorN& oldPositions, VectorN& velocities, Real timestep, nlohmann::json& state, unsigned int frame, std::vector<unsigned int>& collVerts);
		// void resolveCollisions(std::vector<std::unique_ptr<Obstacle>>& obs, VectorN& updated_pos, const VectorN& oldPositions, VectorN& velocities, Real timestep, nlohmann::json& state, unsigned int frame);
		void getCollisionFreeState(std::vector<std::unique_ptr<Obstacle>>& obstacles, VectorN& X, const VectorN& Xold, Real timestep, std::vector<int>* collVerts, unsigned int frame);
		void getCollisionFreeState(std::vector<std::unique_ptr<Obstacle>>& obstacles, VectorN& X, const VectorN& Xold, const std::vector<trimesh::triangle_t>& triangles, Real timestep, std::vector<int>* collVerts = NULL, unsigned int frame = 0);
		void getCollisionFreeState(std::vector<std::unique_ptr<Obstacle>>& obstacles, VectorN& X, VectorN& V, const VectorN& Xold, const VectorN& W, const std::vector<trimesh::triangle_t>& triangles, Real timestep, std::vector<int>* collVerts, unsigned int frame);
		void resolveCollisions(std::vector<std::unique_ptr<Obstacle>>& obstacles, VectorN& updated_pos, const VectorN& oldPositions, VectorN& velocities, const VectorN& W, const std::vector<trimesh::triangle_t>& triangles, Real timestep, nlohmann::json& state, unsigned int frame, std::vector<unsigned int>& collVerts);

		void addExternalForces(std::vector<std::unique_ptr<DynamicsObject>>& objects, Real timestep);
		void applyConstraints(Real timestep, std::vector<std::unique_ptr<DynamicsObject>>& objects);
		void updatePos(Real timestep, std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs, nlohmann::json& state, unsigned int frame);
		~VBD();
	};
}

#endif /* VBD */
