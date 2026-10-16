#include <argus/dynamics_object.hpp>
#include <argus/symplectic_euler.hpp>
#include <argus/util.hpp>

namespace argus
{
	// void SymplecticEuler::step(Real timeStep, std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects)

	void SymplecticEuler::getCollisionFreeState(std::vector<std::unique_ptr<Obstacle>>& obstacles, VectorN& X, VectorN& V, const VectorN& Xold, const VectorN& W, const std::vector<trimesh::triangle_t>& triangles, Real timestep, std::vector<int>* collVerts, unsigned int frame, unsigned int substep)
	{
		/**
		 * Project updated_pos to collision free states.
		 */
		int nobs = obstacles.size();
		for (int iter = 0; iter < m_collision_iterations; iter++)
			for (int i = 0; i < nobs; i++)
			{
				if (m_sdfCollider == "basic")
					obstacles.at(i)->project(X, V, Xold, timestep, frame, collVerts, substep);
				else if (m_sdfCollider == "quadrature")
					obstacles.at(i)->project(timestep, X, V, Xold, W, triangles, frame, collVerts, &m_config);
				else
					throw "not implemented!";
			}
	}

	void SymplecticEuler::resolveCollisions(std::vector<std::unique_ptr<Obstacle>>& obs, VectorN& updated_pos, const VectorN& oldPositions, VectorN& velocities, const VectorN& W, const std::vector<trimesh::triangle_t>& triangles, Real timestep, nlohmann::json& state, unsigned int frame, std::vector<unsigned int>& collVerts, unsigned int substep)
	{
		/**
		 * Project sheet positions to collision free states after resolving collisions from all obstacles.
		 */
		if (obs.size() != 0)
		{
			unsigned int coliters = 10;
#ifdef ARGUS_CHECKPOINT
			if (simConf.contains("Obstacle Collision iterations | All objects") == false)
				simConf["Obstacle Collision iterations | All objects"] = coliters;
#endif
			// create a copy if normal info is required.
			std::vector<int> collisionVertices;

			for (int i = 0; i < coliters; i++)
			{
#ifdef ARGUS_CHECKPOINT
				getCollisionFreeState(obs, updated_pos, velocities, oldPositions, W, triangles, timestep, &collisionVertices, frame, substep);
#else
				getCollisionFreeState(obs, updated_pos, velocities, oldPositions, W, triangles, timestep, NULL, frame, substep);
#endif
			}
			// velocities = (updated_pos - oldPositions) / timestep;
#ifdef ARGUS_CHECKPOINT
			// store state
			std::set<int> cset(collisionVertices.begin(), collisionVertices.end());
			std::vector<unsigned int> cvec(cset.begin(), cset.end());
			collVerts = cvec;

			if (collisionVertices.size() == 0)
			{
				state["collision_vertices"] = collisionVertices;
			}
			else
			{
				Eigen::VectorXi collverts;
				collverts.resize(cvec.size());
				for (unsigned int i = 0; i < cvec.size(); i++)
					collverts(i) = cvec.at(i);
				state["collision_vertices"] = collverts;
			}
#endif
		}
	}
	void SymplecticEuler::init(std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obstacles, nlohmann::json& params)
	{
		m_kdamp = 0;
		m_substeps = 1;
		m_config = params;
		m_stats["optimizeDphis"] = 0;
		m_stats["collisions"] = 0;
		m_stats["self_collisions"] = 0;
		m_stats["update"] = 0;
		m_stats["total"] = 0;
		for (int i = 0; i < dynamicObjects.size(); i++)
		{
			int nverts = dynamicObjects[i]->getVertexCount();
			solvers.emplace_back(std::make_unique<Eigen::SimplicialLDLT<Eigen::SparseMatrix<Real>>>());
			// setup the mass matrices
			Eigen::SparseMatrix<Real> M(3 * nverts, 3 * nverts);
			std::vector<Eigen::Triplet<Real>> triplets;
			triplets.reserve(3 * nverts);
			for (int j = 0; j < nverts; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					// setup triplet
					triplets.push_back(Eigen::Triplet<Real>(3 * j + k, 3 * j + k, dynamicObjects[i]->getInertia()[j]));
				}
			}
			M.setFromTriplets(triplets.begin(), triplets.end());
			M.makeCompressed();
			MassMatrix.emplace_back(M);
		}

		if (params.contains("substeps"))
			m_substeps = params["substeps"].get<Real>();

		if (params.contains("kdamp_pos"))
			m_kdamp = params["kdamp_pos"].get<Real>();
		if (params.contains("obstacle_collision_maxiter"))
		{
			m_collision_iterations = params["obstacle_collision_maxiter"].get<unsigned int>();
		}
		else
		{
			m_collision_iterations = 1;
		}
		if (params.contains("sdf_collider"))
			m_sdfCollider = params["sdf_collider"].get<std::string>();
		else
			m_sdfCollider = "basic";

		if (params.contains("self_collisions_maxsubsteps"))
			selfCollisionParams.self_collision_maxsubsteps = params["self_collisions_maxsubsteps"].get<int>();
		else
			selfCollisionParams.self_collision_maxsubsteps = 1;

		if (params.contains("self_collisions"))
			selfCollisionParams.enable_self_collisions = params["self_collisions"].get<bool>();
		else
			selfCollisionParams.enable_self_collisions = false;
		if (params.contains("self_collision_ccd"))
		{
			selfCollisionParams.doccd = params["self_collision_ccd"].get<bool>();
		}
		else
			selfCollisionParams.doccd = false;
		if (params.contains("self_collision_stiffness"))
		{
			selfCollisionParams.stiffness = params["self_collision_stiffness"].get<Real>();
		}
		else
			selfCollisionParams.stiffness = 0.01;
		if (params.contains("self_collision_damping"))
		{
			selfCollisionParams.damping = params["self_collision_damping"].get<Real>();
		}
		else
			selfCollisionParams.damping = 0.01;

		if (params.contains("self_collision_threshold"))
		{
			selfCollisionParams.threshold = params["self_collision_threshold"].get<Real>();
		}
		else
			selfCollisionParams.threshold = 0.01;

		initSelfCollision(dynamicObjects[0]->getTriangles(), dynamicObjects[0]->getRestPositions(), selfCollisionData);
	}

	void SymplecticEuler::step(Real timeStep, std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs, nlohmann::json& state, unsigned int frame)
	{
		Real subdt = timeStep / m_substeps;
		Real substep_kdamp = 1.0 - std::pow(1.0 - m_kdamp, 1.0 / m_substeps);
		auto ttotal1 = argus::util::captureCurrentTime();

		int nobs = obs.size();
		std::vector<unsigned int> collverts;
		int counter = 0;
		for (auto&& object : dynamicObjects)
		{
			// std::cout << "Time step " << timeStep << std::endl;
			// std::cout << "Damping force " << m_kdamp << std::endl;
			VectorN x, xt, xtsub, vnew, vt;
			int nverts = object->getVertexCount();
			xt = object->getPositions();
			// vt = object->getVelocities();
			x = xt;
			vnew = object->getVelocities();
			xtsub.resize(3 * nverts);

			auto massInverse = object->getInvInertia().asDiagonal();
			VectorN invInertia = object->getInvInertia();
			if (frame == 0 || frame == 1)
			{

				// set the solvers
				Real alpha1 = 0, alpha2 = m_kdamp;
				std::cout << " using alpha2 in laplacian damping " << alpha2 << std::endl;
				Eigen::SparseMatrix<Real> MM = MassMatrix[counter];
				Eigen::SparseMatrix<Real> A = MM + subdt * (alpha1 * MM) + subdt * (alpha2 * object->getDampingLaplacian());
				// A.diagonal().array() += 1e-8;
				A.makeCompressed();
				solvers[counter]->compute(A);
			}
			updatePositions(counter, subdt, *object, obs, state, frame);

			counter++;
		}
		auto ttotal2 = argus::util::captureCurrentTime();
		m_stats["total"] = m_stats.value("total", 0LL) + std::chrono::duration_cast<std::chrono::milliseconds>(ttotal2 - ttotal1).count();
	}

	void SymplecticEuler::updatePositions(int objid, Real timeStep, DynamicsObject& object, std::vector<std::unique_ptr<argus::Obstacle>>& obs, nlohmann::json& state, unsigned int frame)
	{

		auto massInverse = object.getInvInertia().asDiagonal();
		const VectorN& invInertia = object.getInvInertia();
		Map<VectorN> forces = object.getForces();
		VectorN x = object.getPositions();
		VectorN xtsub, vnew;
		vnew = object.getVelocities();
		int nobs = obs.size();
		auto tupdate1 = argus::util::captureCurrentTime();
		Real substep_kdamp = 1.0 - std::pow(1.0 - m_kdamp, 1.0 / m_substeps);

		for (int substep = 0; substep < m_substeps; substep++)
		{
			xtsub = x;
			object.getPositions() = x;
			Map<VectorN> forces = object.getForces();

			// apply damping
			if (solvers[objid]->info() == Eigen::Success)
			{
				VectorN rhs = MassMatrix[objid] * vnew + timeStep * forces;
				vnew = solvers[objid]->solve(rhs);
				object.applyConstraints(Map<VectorN>(vnew.data(), vnew.size()));
			}
			else
			{
				std::cerr << "Laplacian damping solver failed to factorize! FAlling back to ether\n";
				VectorN Fn = massInverse * forces;
				vnew += timeStep * Fn;
				vnew *= (1 - substep_kdamp);
			}

			// use when damping is not needed
			/* 		VectorN Fn
			            = massInverse * forces;
			        vnew += timeStep * Fn; */

			x += timeStep * vnew;
			auto tupdate2 = argus::util::captureCurrentTime();

			m_stats["update"] = m_stats.value("update", 0LL) + std::chrono::duration_cast<std::chrono::milliseconds>(tupdate2 - tupdate1).count();
			auto tselfcollisions1 = argus::util::captureCurrentTime();

			if (selfCollisionParams.enable_self_collisions)
			{
				resolveSelfCollisionsSubstepped(
				    object.getTriangles(),
				    xtsub, // Position at the very start of the frame
				    x, // Predicted position (modified to resolved position)
				    vnew, // Velocity (modified)
				    selfCollisionParams.self_collision_maxsubsteps,
				    invInertia,
				    selfCollisionData,
				    selfCollisionParams.stiffness, // stiffness
				    selfCollisionParams.damping, // damping
				    timeStep);
				vnew = (x - xtsub) / timeStep; // Is it required here?
			}
			auto tselfcollisions2 = argus::util::captureCurrentTime();
			m_stats["self_collisions"] = m_stats.value("self_collisions", 0LL) + std::chrono::duration_cast<std::chrono::milliseconds>(tselfcollisions2 - tselfcollisions1).count();
			auto tcollisions1 = argus::util::captureCurrentTime();
			if (nobs > 0)
			{
				std::vector<unsigned int> collverts;
				resolveCollisions(obs, x, xtsub, vnew, invInertia, object.getTriangles(), timeStep, state, frame, collverts, substep);
			}
			auto tcollisions2 = argus::util::captureCurrentTime();
			m_stats["collisions"] = m_stats.value("collisions", 0LL) + std::chrono::duration_cast<std::chrono::milliseconds>(tcollisions2 - tcollisions1).count();
			object.getPositions() = x;

			object.applyDynamicConstraints(timeStep);
			x = object.getPositions();

			vnew = (x - xtsub) / timeStep;

			// auto tupdate3 = argus::util::captureCurrentTime();

			// x = xtsub + timeStep * vnew;
			// auto tupdate4 = argus::util::captureCurrentTime();
			// m_stats["update"] = m_stats.value("update", 0LL) + std::chrono::duration_cast<std::chrono::milliseconds>(tupdate4 - tupdate3).count();
		}
		object.getVelocities() = vnew;
		object.getPositions() = x;
	}
}
