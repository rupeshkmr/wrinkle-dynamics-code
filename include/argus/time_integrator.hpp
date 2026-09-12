#ifndef TIME_INTEGRATOR_HPP
#define TIME_INTEGRATOR_HPP

#include "common.hpp"
#include "dynamics_object.hpp"
#include "obstacles.hpp"
#include <vector>
#include <memory>
#include <variant>

namespace argus
{
	class TimeIntegrator
	{
	public:
		// Must be called after all objects have been added to world
		// else first frame computations will take more time
		virtual void init(std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obstacles, nlohmann::json& params) { throw "Not Implemented"; };
		// virtual void step(Real timeStep,
		//     std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects)
		virtual void step(Real timeStep, std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs, nlohmann::json& state, unsigned int frame = 0) = 0;
		virtual nlohmann::json getStats() { throw "Not Implemented"; }
		virtual ~TimeIntegrator() = default;
	};
}

#endif /* TIME_INTEGRATOR_HPP */
