#include <argus/common.hpp>
#include <argus/symplectic_euler_static_wrinkles.hpp>
#include <argus/lbfgsb_wrapper.hpp>
#include <argus/lbfgs_wrapper.hpp>
#include <argus/util.hpp>
namespace argus
{
	void SymplecticEulerStaticWrinkles::initWrinkleMesh(std::vector<std::unique_ptr<DynamicsObject>>& objects, bool amps, bool dphis)
	{
		if (amps == false && dphis == false)
			throw "Error! Calling initWrinkleMesh without chosing a parameter.\n";
		for (auto&& object : objects)
		{
			// TODO: Remove copying here
			VectorN a;
			Matrix2N w;
			object->getInitialWrinkleParameters(0.1, a, w, amps);
			if (amps)
				object->getAmpVec() = a;
			if (dphis)
				object->getDphis() = w;
		}
	}
	void SymplecticEulerStaticWrinkles::optimize(std::vector<std::unique_ptr<DynamicsObject>>& objects, bool amps, bool dphis, bool positions)
	{
		if (!(amps || dphis || positions))
			throw "Error, enable atleast one sheet parameter to be optimized!\n";
		if (amps)
		{
			for (auto&& object : objects)
			{
				unsigned int vertexCount, faceCount;
				vertexCount = object->getVertexCount();
				faceCount = object->getFaceCount();
				std::function<Real(VectorN & X, VectorN & deriv)> func = [&object, vertexCount, faceCount](VectorN& X, VectorN& deriv) -> Real
				{
					VectorN ip, fderiv;
					ip.resize(4 * vertexCount + 2 * faceCount);
					fderiv.resize(4 * vertexCount + 2 * faceCount);
					ip(Eigen::seq(0, vertexCount - 1)) = X;
					ip(Eigen::seq(vertexCount, vertexCount + 2 * faceCount - 1)) = object->getDphis().reshaped();
					ip(Eigen::seq(vertexCount + 2 * faceCount, 4 * vertexCount + 2 * faceCount - 1)) = object->getPositions().reshaped();
					Real e = object->getEnergies(ip, &fderiv);
					deriv = fderiv(Eigen::seq(0, vertexCount - 1));
					object->clampDofs(deriv);
					return e;
				};

				VectorN lb, ub;
				ub.resize(vertexCount);
				lb.resize(vertexCount);
				// set the bounds
				for (int i = 0; i < vertexCount; i++)
				{
					ub(i) = std::numeric_limits<argus::Real>::infinity();
					lb(i) = 0;
				}
				VectorN X = object->getAmpVec();
				if (m_boundedOptimizer)
				{
					try
					{
						m_boundedOptimizer->solve(func, X, &lb, &ub);
					}
					catch (const char* err)
					{
						std::cerr << "Error while optimizing for amps! " << err << std::endl;
					}
					catch (...)
					{
						std::cerr << "Other error while optimizing for amps! " << std::endl;
					}
				}
				else
					throw "Error bounded optimizer is not initialized before calling optimize!\n";
				object->getAmpVec() = X;
			}
		}

		if (dphis)
		{
			// optimizeDphisParallel(objects);
			for (auto&& object : objects)
			{
				object->optimizeParallelEnergies();
			}
		}

		if (positions)
			throw "Static optimizer not implemented for positions!\n";
	}

	void SymplecticEulerStaticWrinkles::getCollisionFreeState(std::vector<std::unique_ptr<Obstacle>>& obstacles, VectorN& X, VectorN& V, const VectorN& Xold, const VectorN& W, const std::vector<trimesh::triangle_t>& triangles, Real timestep, std::vector<int>* collVerts, unsigned int frame, unsigned int substep)
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

	void SymplecticEulerStaticWrinkles::resolveCollisions(std::vector<std::unique_ptr<Obstacle>>& obs, VectorN& updated_pos, const VectorN& oldPositions, VectorN& velocities, const VectorN& W, const std::vector<trimesh::triangle_t>& triangles, Real timestep, nlohmann::json& state, unsigned int frame, std::vector<unsigned int>& collVerts, unsigned int substep)
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
	void SymplecticEulerStaticWrinkles::init(std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obstacles, nlohmann::json& params)
	{
		m_kdamp_amps = 0;
		m_kdamp_pos = 0;
		m_substeps = 1;
		m_config = params;
		m_stats["optimizeDphis"] = 0;
		m_stats["collisions"] = 0;
		m_stats["self_collisions"] = 0;
		m_stats["symplectic_update"] = 0;
		m_stats["misc"] = 0;
		m_stats["total_time_taken_by_time_integrator"] = 0;
		m_stats["call_to_pos_update"] = 0;
		// set mass matrix
		for (auto&& it : dynamicObjects)
		{
			VectorN massInv = it->getInvInertia();
			int nverts = it->getVertexCount();
			Eigen::SparseMatrix<Real> ampmassinv(nverts, nverts);

			for (int i = 0; i < nverts; i++)
			{
				ampmassinv.insert(i, i) = massInv(3 * i); // std::min(std::min(massInv(3 * i), massInv(3 * i + 1)), massInv(3 * i + 2));
			}
			ampmassinv.makeCompressed();
			m_ampInvMassMatrices.emplace_back(ampmassinv);
		}
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
		if (params.contains("kdamp_amps"))
			m_kdamp_amps = params["kdamp_amps"].get<Real>();
		if (params.contains("kdamp_pos"))
			m_kdamp_pos = params["kdamp_pos"].get<Real>();
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
		Real tolerance = 1e-6;
		int maxiter = 100;
		if (params.contains("bounded_optimizer"))
		{
			if (params.contains("bounded_optimizer_tolerance"))
				tolerance = params["bounded_optimizer_tolerance"].get<Real>();
			if (params.contains("bounded_optimizer_maxiter"))
				maxiter = params["bounded_optimizer_maxiter"].get<Real>();
			m_boundedOptimizer = std::make_unique<LBFGSBWrapper>();
			m_boundedOptimizer->init(tolerance, maxiter);
		}
		if (params.contains("optimizer"))
		{
			if (params.contains("optimizer_tolerance"))
				tolerance = params["optimizer_tolerance"].get<Real>();
			if (params.contains("optimizer_maxiter"))
				maxiter = params["optimizer_maxiter"].get<Real>();
			m_optimizer = std::make_unique<LBFGSWrapper>();
			m_optimizer->init(tolerance, maxiter);
		}
		initSelfCollision(dynamicObjects[0]->getTriangles(), dynamicObjects[0]->getRestPositions(), selfCollisionData);
		initWrinkleMesh(dynamicObjects, true, true);
		optimize(dynamicObjects, true, true, false);
	}

	void SymplecticEulerStaticWrinkles::step(Real timeStep, std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs, nlohmann::json& state, unsigned int frame)
	{
		// m_stats["optimizeDphis"] = 0;

		auto ttotal1 = argus::util::captureCurrentTime();

		int nobs = obs.size();
		int counter = 0;
		for (auto&& object : dynamicObjects)
		{
			// m_stats["misc"] = 0;
			auto tmisc1 = argus::util::captureCurrentTime();

			int nverts = object->getVertexCount();
			int nfaces = object->getFaceCount();
			VectorN x, xt, v;
			x.resize(3 * nverts);
			xt.resize(3 * nverts);
			v.resize(3 * nverts);
			VectorN xold = object->getPositions();
			// VectorN invInertia = object->getInvInertia();
			const Eigen::SparseMatrix<Real>& invInertiaAmps = m_ampInvMassMatrices[counter];
			std::vector<unsigned int> collverts;
			// auto massInverse = invInertia.asDiagonal();
			Real subdt = timeStep / m_substeps;

			if (frame == 0 || frame == 1)
			{
				// set the solvers
				Real alpha1 = 0.0, alpha2 = m_kdamp_pos;
				Eigen::SparseMatrix<Real> MM = MassMatrix[counter];
				Eigen::SparseMatrix<Real> A = MM + /* subdt * (alpha1 * MM) + */ subdt * (alpha2 * object->getDampingLaplacian());
				A.diagonal().array() += 1e-8;
				A.makeCompressed();
				solvers[counter]->compute(A);
			}
			x = object->getPositions();
			auto tmisc2 = argus::util::captureCurrentTime();
			m_stats["misc"] = m_stats.value("misc", 0LL) + std::chrono::duration_cast<std::chrono::milliseconds>(tmisc2 - tmisc1).count();
			// curtain

			// std::vector<int> tV = { 0, 1, 5, 10, 15, 26, 31, 42, 47, 82, 87, 98, 103, 138, 143, 154, 159, 290, 307, 315, 367, 375, 532, 555, 607, 612, 699, 739, 993, 1003, 1059, 1067, 1088 };

			for (int substep = 0; substep < m_substeps; substep++)
			{
				auto tupdate1 = argus::util::captureCurrentTime();

				VectorN forces = object->getWrinkleForces();
				// update amps
				// VectorN ampForces = 2.0 * forces.segment(0, nverts);
				// object->getAmpVelocities() += subdt * invInertiaAmps * ampForces;
				// object->getAmpVelocities() *= (1 - m_kdamp_amps);
				// object->getAmpVec() += subdt * object->getAmpVelocities();
				// update pos
				VectorN posForces = forces.segment(nverts + 2 * nfaces, 3 * nverts);
				// object->getVelocities() += subdt * (massInverse * posForces);
				// object->getVelocities() *= (1 - m_kdamp_pos);
				// object->getPositions() += subdt * object->getVelocities();
				xt = x;
				v = object->getVelocities();
				auto tupdate2 = argus::util::captureCurrentTime();
				m_stats["symplectic_update"] = m_stats.value("symplectic_update", 0LL) + std::chrono::duration_cast<std::chrono::milliseconds>(tupdate2 - tupdate1).count();
				// m_stats["call_to_pos_update"]
				auto cpu1 = argus::util::captureCurrentTime();

				updatePositions(counter, subdt, x, v, xt, posForces, *object, obs, state, frame, substep);
				auto cpu2 = argus::util::captureCurrentTime();
				m_stats["call_to_pos_update"] = m_stats.value("call_to_pos_update", 0LL) + std::chrono::duration_cast<std::chrono::milliseconds>(cpu2 - cpu1).count();

				auto tmisc1 = argus::util::captureCurrentTime();

				object->getPositions() = x;
				/* for (int i = 0; i < tV.size(); i++)
				{
				    int vviid = tV[i];
				    object->getPositions().segment<2>(3 * vviid + 1) = xold.segment<2>(3 * vviid + 1);
				} */
				object->getVelocities() = v;
				object->applyDynamicConstraints(subdt);
				auto tmisc2 = argus::util::captureCurrentTime();
				m_stats["misc"] = m_stats.value("misc", 0LL) + std::chrono::duration_cast<std::chrono::milliseconds>(tmisc2 - tmisc1).count();

				object->applyTopologicalConstraints();
				// update dphis
				auto toptimizeDphis1 = argus::util::captureCurrentTime();
				// auto oldDphis = object->getDphis();
				initWrinkleMesh(dynamicObjects, true, true);
				// object->getDphis() = (object->getDphis() + oldDphis) / 2.0;
				optimize(dynamicObjects, true, true, false);
				auto toptimizeDphis2 = argus::util::captureCurrentTime();
				m_stats["optimizeDphis"] = m_stats.value("optimizeDphis", 0LL) + std::chrono::duration_cast<std::chrono::milliseconds>(toptimizeDphis2 - toptimizeDphis1).count();
			}
			counter++;
		}
		auto ttotal2 = argus::util::captureCurrentTime();
		m_stats["total_time_taken_by_time_integrator"] = m_stats.value("total_time_taken_by_time_integrator", 0LL) + std::chrono::duration_cast<std::chrono::milliseconds>(ttotal2 - ttotal1).count();
	}
	void SymplecticEulerStaticWrinkles::updatePositions(int objid, Real timeStep,
	    VectorN& x, VectorN& v, const VectorN& xt, const VectorN& forces,
	    DynamicsObject& object, std::vector<std::unique_ptr<argus::Obstacle>>& obs,
	    nlohmann::json& state, unsigned int frame, int step)
	{
		// m_stats["symplectic_update"] = m_stats.value("symplectic_update", 0LL) + std::chrono::duration_cast<std::chrono::milliseconds>(tupdate2 - tupdate1).count();
		auto tsupdate1 = argus::util::captureCurrentTime();

		auto& massInverse = object.getInvInertia().asDiagonal();
		auto invInertia = object.getInvInertia();
		int nobs = obs.size();

		// apply damping
		if (solvers[objid]->info() == Eigen::Success)
		{
			VectorN rhs = MassMatrix[objid] * v + timeStep * forces;
			v = solvers[objid]->solve(rhs);
			object.applyConstraints(Map<VectorN>(v.data(), v.size()));
		}
		else
		{
			std::cerr << "Laplacian damping solver failed to factorize! FAlling back to ether\n";
			VectorN Fn = massInverse * forces;
			v += timeStep * Fn;
			v *= (1 - m_kdamp_pos);
		}

		x += timeStep * v;
		auto tsupdate2 = argus::util::captureCurrentTime();
		m_stats["symplectic_update"] = m_stats.value("symplectic_update", 0LL) + std::chrono::duration_cast<std::chrono::milliseconds>(tsupdate2 - tsupdate1).count();

		auto tselfcollisions1 = argus::util::captureCurrentTime();

		if (selfCollisionParams.enable_self_collisions)
		{
			resolveSelfCollisionsSubstepped(
			    object.getTriangles(),
			    xt, // Position at the very start of the frame
			    x, // Predicted position (modified to resolved position)
			    v, // Velocity (modified)
			    selfCollisionParams.self_collision_maxsubsteps,
			    invInertia,
			    selfCollisionData,
			    selfCollisionParams.stiffness, // stiffness
			    selfCollisionParams.damping, // damping
			    timeStep);
			v = (x - xt) / timeStep; // Is it required here?
		}
		auto tselfcollisions2 = argus::util::captureCurrentTime();
		m_stats["self_collisions"] = m_stats.value("self_collisions", 0LL) + std::chrono::duration_cast<std::chrono::milliseconds>(tselfcollisions2 - tselfcollisions1).count();

		auto tcollisions1 = argus::util::captureCurrentTime();
		if (nobs > 0)
		{
			std::vector<unsigned int> collverts;
			resolveCollisions(obs, x, xt, v, invInertia, object.getTriangles(), timeStep, state, frame, collverts, step);
		}
		auto tcollisions2 = argus::util::captureCurrentTime();
		m_stats["collisions"] = m_stats.value("collisions", 0LL) + std::chrono::duration_cast<std::chrono::milliseconds>(tcollisions2 - tcollisions1).count();
	}
}
