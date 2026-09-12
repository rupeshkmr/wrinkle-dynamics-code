#ifndef SYMPLECTIC_EULER_HPP_WRINKLES
#define SYMPLECTIC_EULER_HPP_WRINKLES

#include "common.hpp"
#include "self_collisions_helper.hpp"
#include "optimizers.hpp"
#include "time_integrator.hpp"

namespace argus
{
	class SymplecticEulerWrinkles : public TimeIntegrator
	{
		Real m_kdamp_pos, m_kdamp_amps, m_substeps;
		std::unique_ptr<Optimizer> m_optimizer;
		std::unique_ptr<Optimizer> m_boundedOptimizer;
		std::string m_sdfCollider;
		int m_collision_iterations;
		nlohmann::json m_stats, m_config;
		struct scParams selfCollisionParams;
		PrecomputedSelfCollisionData selfCollisionData;
		std::vector<Eigen::SparseMatrix<Real>> m_ampInvMassMatrices;
		std::vector<std::unique_ptr<Eigen::SimplicialLDLT<Eigen::SparseMatrix<Real>>>> solvers;
		std::vector<Eigen::SparseMatrix<Real>> MassMatrix;

	public:
		void init(std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obstacles, nlohmann::json& params) override;
		// void step(Real timeStep,
		//     std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects) override;
		void step(Real timeStep, std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs, nlohmann::json& state, unsigned int frame = 0) override;
		void initWrinkleMesh(std::vector<std::unique_ptr<DynamicsObject>>& objects, bool amps, bool dphis);
		void optimize(std::vector<std::unique_ptr<DynamicsObject>>& objects, bool amps, bool dphis, bool positions);
		void resolveCollisions(std::vector<std::unique_ptr<Obstacle>>& obs, VectorN& updated_pos, const VectorN& oldPositions, VectorN& velocities, const VectorN& W, const std::vector<trimesh::triangle_t>& triangles, Real timestep, nlohmann::json& state, unsigned int frame, std::vector<unsigned int>& collVerts, unsigned int substep);
		void getCollisionFreeState(std::vector<std::unique_ptr<Obstacle>>& obstacles, VectorN& X, VectorN& V, const VectorN& Xold, const VectorN& W, const std::vector<trimesh::triangle_t>& triangles, Real timestep, std::vector<int>* collVerts, unsigned int frame, unsigned int substep);
		void updatePositions(int objid, Real timeStep, VectorN& x, VectorN& v, const VectorN& xt, const VectorN& forces, DynamicsObject& object, std::vector<std::unique_ptr<argus::Obstacle>>& obs, nlohmann::json& state, unsigned int frame, int step);
		nlohmann::json getStats() override { return m_stats; }
	};
}

#endif /* SYMPLECTIC_EULER_HPP FOR WRINKLES */
