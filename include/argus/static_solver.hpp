#ifndef STATIC_SOLVER_HPP
#define STATIC_SOLVER_HPP

#include "common.hpp"
#include "dynamics_object.hpp"
#include "obstacles.hpp"
#include <vector>
#include <memory>
#include <variant>

namespace argus
{
	class StaticSolver
	{
	public:
		/** Must be called after all objects have been added to world
		 * else first frame computations will take more time
		 * assuming lbfgs solver, however the parameters will be similar for other static solvers
		 * the lower and upper bounds are for box constrained solvers
		 */
		// Must be called after all objects have been added to world
		// else first frame computations will take more time
		virtual void init(std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, nlohmann::json& params) { };
		virtual void solve(std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs, nlohmann::json& state) = 0;
		virtual ~StaticSolver() = default;
	};
}

#endif /* STATIC SOLVER */
