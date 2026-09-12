#ifndef WORLD_HPP
#define WORLD_HPP

#include "common.hpp"

#include <vector>
#include <memory>
#include "obstacles.hpp"
#include "dynamics_object.hpp"
#include "time_integrator.hpp"
#include "static_solver.hpp"
namespace argus
{
	// Shift this to some other file after we have more
	// integrators
	enum class Integrator
	{
		SymplecticEuler = 0,
		WrinklesVBD = 1,
		VBD = 2,
		SymplecticEulerWrinkles = 3,
		SymplecticEulerStaticWrinkles = 4,
	};

	// The idea here is to have solvers that can be used by the world for static solve.
	enum class Solver
	{
		WrinkleSolver = 0,
	};

	class World
	{
		// Currently assuming world to be made of only
		// square sheets and all are located at origin
		// World takes ownership of sheets
		// What about ordering information?
		std::vector<std::unique_ptr<DynamicsObject>> m_DynamicObjects;
		std::vector<std::unique_ptr<Obstacle>> m_Obstacles;

		// Assuming only ExplicitEuler Integrator
		std::unique_ptr<TimeIntegrator> m_Integrator;
		// Static solver
		std::unique_ptr<StaticSolver> m_Solver;
		// Optionally Store state
		nlohmann::json m_State;

	public:
		DynamicsObject& addDynamic(std::unique_ptr<DynamicsObject> sheet);
		Obstacle& addObstacle(std::unique_ptr<Obstacle> obstacle);
		Obstacle& replaceObstacle(std::unique_ptr<Obstacle> obstacle);
		std::vector<std::unique_ptr<Obstacle>>& getObstacles() { return m_Obstacles; }
		void setIntegrator(std::string name);
		void setSolver(std::string name);
		void setIntegrator(Integrator integrator);
		void setSolver(Solver solver);
		// Call before the first step()
		// and after making any changes to the world
		void prepare(nlohmann::json& params);
		void prepare(Real tolerance, int max_iter = 100) { throw "Static solver is not implemented!"; } // for static solver
		void solve(); // time integrator
		void step(Real deltaTime, unsigned int frame = 0); // time integrator
		void saveCheckpoint(std::vector<std::string> data_to_save, std::string parent_directory, int frame);
		nlohmann::json getTimeStepStats() { return m_Integrator->getStats(); }
	};
}

#endif /* WORLD_HPP */
