#ifndef SYMPLECTIC_EULER_HPP
#define SYMPLECTIC_EULER_HPP
#include "common.hpp"
#include "time_integrator.hpp"
#include "self_collisions_helper.hpp"

namespace argus
{
	class SymplecticEuler : public TimeIntegrator
	{
		Real m_kdamp;
		int m_substeps, m_collision_iterations;
		std::string m_sdfCollider;
		nlohmann::json m_stats, m_config;
		struct scParams selfCollisionParams;
		PrecomputedSelfCollisionData selfCollisionData;
		std::vector<std::unique_ptr<Eigen::SimplicialLDLT<Eigen::SparseMatrix<Real>>>> solvers;
		std::vector<Eigen::SparseMatrix<Real>> MassMatrix;

	public:
		// void step(Real timeStep,
		//     std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects) override;
		void init(std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obstacles, nlohmann::json& params) override;
		// {
		// 	m_kdamp = 0;
		// 	if (params.contains("kdamp_pos"))
		// 		m_kdamp = params["kdamp_pos"].get<Real>();
		// 	m_substeps = 1;
		// 	if (params.contains("substeps"))
		// 		m_substeps = params["substeps"].get<Real>();
		// };
		void resolveCollisions(std::vector<std::unique_ptr<Obstacle>>& obs, VectorN& updated_pos, const VectorN& oldPositions, VectorN& velocities, const VectorN& W, const std::vector<trimesh::triangle_t>& triangles, Real timestep, nlohmann::json& state, unsigned int frame, std::vector<unsigned int>& collVerts, unsigned int substep = 0);
		void getCollisionFreeState(std::vector<std::unique_ptr<Obstacle>>& obstacles, VectorN& X, VectorN& V, const VectorN& Xold, const VectorN& W, const std::vector<trimesh::triangle_t>& triangles, Real timestep, std::vector<int>* collVerts, unsigned int frame, unsigned int substep = 0);

		void step(Real timeStep, std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs, nlohmann::json& state, unsigned int frame = 0) override;
		void updatePositions(int objid, Real timeStep, DynamicsObject& object, std::vector<std::unique_ptr<argus::Obstacle>>& obs, nlohmann::json& state, unsigned int frame);
		nlohmann::json getStats() override { return m_stats; }
	};
}

#endif /* SYMPLECTIC_EULER_HPP */
