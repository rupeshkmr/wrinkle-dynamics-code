#ifndef WRINKLES_SOLVER_HPP
#define WRINKLES_SOLVER_HPP
#include "static_solver.hpp"
#include "optimizers.hpp"

namespace argus
{
	class WrinkleSolver : public StaticSolver
	{
		std::unique_ptr<Optimizer> m_optimizer;
		std::unique_ptr<Optimizer> m_boundedOptimizer;

	public:
		void init(std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, nlohmann::json& params) override;
		void solve(std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs, nlohmann::json& state) override;
		void initWrinkleMesh(std::vector<std::unique_ptr<DynamicsObject>>& objects, bool amps, bool dphis);
		void optimize(std::vector<std::unique_ptr<DynamicsObject>>& objects, bool amps, bool dphis, bool positions);
		void optimizeDphisParallel(std::vector<std::unique_ptr<DynamicsObject>>& objects);
		void optimizeDphisGlobal(std::vector<std::unique_ptr<DynamicsObject>>& objects);
		void optimizeAmplitudesGlobal(std::vector<std::unique_ptr<DynamicsObject>>& objects);
		void optimizeAmplitudesNew(std::vector<std::unique_ptr<DynamicsObject>>& objects);
		void optimizeAmplitudesVBD(std::vector<std::unique_ptr<DynamicsObject>>& objects);
		VectorN getAdjacentDphis(int vid, std::unique_ptr<DynamicsObject>& dynamicObject);
		void setAdjacentDphis(int vid, const VectorN& dphis, std::unique_ptr<DynamicsObject>& dynamicObject);
		VectorN getNewDphi(VectorN* gradDphis, Real avals, unsigned int vid, std::unique_ptr<DynamicsObject>& dynamicObject);
	};
}
#endif
