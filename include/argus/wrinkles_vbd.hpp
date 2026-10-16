#ifndef WRINKLE_VBD_HPP
#define WRINKLE_VBD_HPP

#include "common.hpp"
#include "self_collisions_helper.hpp"
#include "time_integrator.hpp"
#include "obstacles.hpp"
#include "optimizers.hpp"
namespace argus
{

	class WrinklesVBD : public TimeIntegrator
	{
		SharedVector m_B; // provides a container of VectorN
		Real m_kdamp_pos, m_kdamp_amps;
		bool m_sim_friction;
		std::string m_ampPosUpdate;
		unsigned int m_amp_vbd_iterations;
		unsigned int m_pos_vbd_iterations;
		unsigned int m_collision_iterations;
		std::unique_ptr<Optimizer> m_optimizer;
		std::unique_ptr<Optimizer> m_boundedOptimizer;
		PrecomputedSelfCollisionData selfCollisionData;

		nlohmann::json m_stats, m_config;
		std::string m_sdfCollider;
		struct scParams selfCollisionParams;

	public:
		// void init(std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, nlohmann::json& params) override;
		void init(std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obstacles, nlohmann::json& params) override;

		void step(Real timestep, std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs, nlohmann::json& state, unsigned int frame = 0) override;
		nlohmann::json getStats() override { return m_stats; }
		void initVBD(std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects);
		void vbdUpdate(Real timestep, std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs);
		void vbdUpdateAmpPosCombined(Real timestep, std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs, nlohmann::json& state, unsigned int frame);
		void vbdUpdateAmpPosInSameIter(Real timestep, std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs, nlohmann::json& state, unsigned int frame);
		void vbdUpdateAmpsDphis(Real timestep, std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs);
		void vbdUpdateAmps(Real timestep, std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs);
		void vbdUpdatePos(Real timestep, std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs, nlohmann::json& state, unsigned int frame);
		VectorN getAdaptiveAccelerationForVBD(DynamicsObject& object, VectorN& var, bool amp, bool pos);
		void resolveCollisions(std::vector<std::unique_ptr<Obstacle>>& obstacles, VectorN& updated_pos, const VectorN& oldPositions, VectorN& velocities, Real timestep, nlohmann::json& state, unsigned int frame, std::vector<unsigned int>& collVerts);
		void resolveCollisions(std::vector<std::unique_ptr<Obstacle>>& obstacles, VectorN& updated_pos, const VectorN& oldPositions, VectorN& velocities, const VectorN& W, const std::vector<trimesh::triangle_t>& triangles, Real timestep, nlohmann::json& state, unsigned int frame, std::vector<unsigned int>& collVerts);
		void getCollisionFreeState(std::vector<std::unique_ptr<Obstacle>>& obstacles, VectorN& X, const VectorN& Xold, Real timestep, std::vector<int>* collVerts = NULL, unsigned int frame = 0);
		void getCollisionFreeState(std::vector<std::unique_ptr<Obstacle>>& obstacles, VectorN& X, const VectorN& Xold, const std::vector<trimesh::triangle_t>& triangles, Real timestep, std::vector<int>* collVerts = NULL, unsigned int frame = 0);
		void getCollisionFreeState(std::vector<std::unique_ptr<Obstacle>>& obstacles, VectorN& X, VectorN& V, const VectorN& Xold, const VectorN& W, const std::vector<trimesh::triangle_t>& triangles, Real timestep, std::vector<int>* collVerts, unsigned int frame);

		void addExternalForces(std::vector<std::unique_ptr<DynamicsObject>>& objects, Real timestep);
		void applyConstraints(Real timestep, std::vector<std::unique_ptr<DynamicsObject>>& objects);
		void initWrinkleMesh(std::vector<std::unique_ptr<DynamicsObject>>& objects, bool amps = false, bool dphis = false);
		void optimize(std::vector<std::unique_ptr<DynamicsObject>>& objects, bool amps = false, bool dphis = false, bool positions = false);
		void optimizeDphisParallel(std::vector<std::unique_ptr<DynamicsObject>>& objects);
		void optimizeDphisParallel(std::vector<std::unique_ptr<DynamicsObject>>& objects, const std::vector<int>& fids);
		void updatePos(Real timestep, std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs, nlohmann::json& state, unsigned int frame);
		void Ehat(int vid);
		void initMeshes(std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obstacles);

		/*     def Ehat(self, ip):
		        a = ip[1]
		        adofs = self.sheet.activeDofs
		        self.sheet.updateActiveDofs([1, 2])
		        ip[1] = a
		        ip[2] = self.omega(a)
		        e = self.sheet.E(ip)
		        self.sheet.updateActiveDofs(adofs)
		        return e

		    def gradEhat(self, ip):
		        a = ip[1]
		        adofs = self.sheet.activeDofs
		        self.sheet.updateActiveDofs([1, 2])
		        ip[0] = self.x[0]
		        ip[1] = a
		        ip[2] = self.omega(a)
		        g = self.sheet.Ea(ip) + self.sheet.Ef(ip) * self.gradOmega(a)
		        self.sheet.updateActiveDofs(adofs)
		        return np.array([0, g, 0]) */

		~WrinklesVBD();
	};
}

#endif /* WRINKLES_VBD */
