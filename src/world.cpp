#include <argus/world.hpp>
#include <argus/symplectic_euler.hpp>
#include <argus/wrinkles_solver.hpp>
#include <argus/wrinkles_vbd.hpp>
#include <argus/symplectic_euler_wrinkles.hpp>
#include <argus/symplectic_euler_static_wrinkles.hpp>
#include <argus/vbd.hpp>
#include <argus/sheet_state_io.hpp>

namespace argus
{
	DynamicsObject& World::addDynamic(std::unique_ptr<DynamicsObject> sheet)
	{
		m_DynamicObjects.push_back(std::move(sheet));
		return *m_DynamicObjects[m_DynamicObjects.size() - 1];
	}
	Obstacle& World::addObstacle(std::unique_ptr<Obstacle> obstacle)
	{
		m_Obstacles.push_back(std::move(obstacle));
		return *m_Obstacles[m_Obstacles.size() - 1];
	}
	Obstacle& World::replaceObstacle(std::unique_ptr<Obstacle> obstacle)
	{
		m_Obstacles.clear();
		m_Obstacles.push_back(std::move(obstacle));
		return *m_Obstacles[m_Obstacles.size() - 1];
	}

	void World::setIntegrator(std::string name)
	{
		if (name == "wrinklesvbd")
			setIntegrator(argus::Integrator::WrinklesVBD);
		else if (name == "vbd")
			setIntegrator(argus::Integrator::VBD);
		else if (name == "symplecticeulerwrinkles")
			setIntegrator(argus::Integrator::SymplecticEulerWrinkles);
		else if (name == "symplecticeulerstaticwrinkles")
			setIntegrator(argus::Integrator::SymplecticEulerStaticWrinkles);
		else if (name == "symplecticeuler")
			setIntegrator(argus::Integrator::SymplecticEuler);
		else
			throw "Integrator not passed!";
	}

	void World::setIntegrator(Integrator integrator)
	{
		switch (integrator)
		{
		case Integrator::SymplecticEuler:
			m_Integrator.reset(new argus::SymplecticEuler());
			break;
		case Integrator::WrinklesVBD:
			m_Integrator.reset(new argus::WrinklesVBD());
			break;
		case Integrator::VBD:
			m_Integrator.reset(new argus::VBD());
			break;
		case Integrator::SymplecticEulerWrinkles:
			m_Integrator.reset(new argus::SymplecticEulerWrinkles());
			break;
		case Integrator::SymplecticEulerStaticWrinkles:
			m_Integrator.reset(new argus::SymplecticEulerStaticWrinkles());
			break;
		}
	}

	void World::setSolver(std::string name)
	{
		if (name == "wrinklesolver")
			setSolver(argus::Solver::WrinkleSolver);
		else
			throw "Solver not passed!";
	}

	void World::setSolver(Solver solver)
	{
		switch (solver)
		{
		case Solver::WrinkleSolver:
			m_Solver.reset(new argus::WrinkleSolver());
			break;
		}
	}

	void World::prepare(nlohmann::json& params)
	{
		// push all cloth objects that are inside the sdf to outside
		if (m_Integrator)
			m_Integrator->init(m_DynamicObjects, m_Obstacles, params);
		if (m_Solver)
			m_Solver->init(m_DynamicObjects, params);
	}

	void World::step(Real deltaTime, unsigned int frame)
	{
		if (m_Integrator)
		{
			m_Integrator->step(deltaTime, m_DynamicObjects, m_Obstacles, m_State, frame);
		}
	}

	void World::solve()
	{
		if (m_Solver)
		{
			m_Solver->solve(m_DynamicObjects, m_Obstacles, m_State);
		}

	} // static solver

	void World::saveCheckpoint(std::vector<std::string> data_to_save, std::string parent_directory, int frame)
	{
		for (auto&& object : m_DynamicObjects)
		{
			object->writeCheckpoint(data_to_save, parent_directory, frame, m_State);
		}
	}
}
