#include <argus/wrinkles_vbd.hpp>
#include <argus/lbfgsb_wrapper.hpp>
#include <argus/lbfgs_wrapper.hpp>
#include <argus/self_collisions.hpp>
#include <argus/util.hpp>
#include <dlib/optimization.h>
namespace argus
{
	WrinklesVBD::~WrinklesVBD() = default;
	void WrinklesVBD::init(std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obstacles, nlohmann::json& params)
	{
		/**
		 * Initialize the timestepper
		 * Current config allows global dynamic damping parameters for all sheets.
		 * In future we may want multiple sheets with different parameters.
		 */

		// default values
		m_kdamp_amps = 0;
		m_kdamp_pos = 0;
		m_amp_vbd_iterations = 5;
		m_pos_vbd_iterations = 5;
		m_sim_friction = false;
		m_config = params;
		// setup stats
		//
		m_stats["vbdUpdateAmps"] = 0;
		m_stats["vbdUpdatePos"] = 0;
		m_stats["optimizeDphis"] = 0;
		m_stats["constraint"] = 0;
		m_stats["init"] = 0;
		m_ampPosUpdate = params.value("ampsposupdate", "separate");
		if (params.contains("obstacle_collision_maxiter"))
		{
			m_collision_iterations = params["obstacle_collision_maxiter"].get<unsigned int>();
		}
		else
		{
			m_collision_iterations = 1;
		}
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

		if (params.contains("sdf_collider"))
			m_sdfCollider = params["sdf_collider"].get<std::string>();
		else
			m_sdfCollider = "basic";
		if (params.contains("kdamp_amps"))
			m_kdamp_amps = params["kdamp_amps"].get<Real>();
		if (params.contains("kdamp_pos"))
			m_kdamp_pos = params["kdamp_pos"].get<Real>();
		if (params.contains("sim_friction"))
		{
			m_sim_friction = params["sim_friction"].get<unsigned int>();
		}
#ifdef ARGUS_CHECKPOINT
		argus::simConf["Friction Enabled"] = m_sim_friction;
#endif
		if (params.contains("bounded_optimizer"))
		{
			Real tolerance = 1e-6;
			int maxiter = 100;
			if (params.contains("bounded_optimizer_tolerance"))
				tolerance = params["bounded_optimizer_tolerance"].get<Real>();
			if (params.contains("bounded_optimizer_maxiter"))
				maxiter = params["bounded_optimizer_maxiter"].get<Real>();
			m_boundedOptimizer = std::make_unique<LBFGSBWrapper>();
			m_boundedOptimizer->init(tolerance, maxiter);
		}
		if (params.contains("optimizer"))
		{
			Real tolerance = 1e-6;
			int maxiter = 100;
			if (params.contains("optimizer_tolerance"))
				tolerance = params["optimizer_tolerance"].get<Real>();
			if (params.contains("optimizer_maxiter"))
				maxiter = params["optimizer_maxiter"].get<Real>();
			m_optimizer = std::make_unique<LBFGSWrapper>();
			m_optimizer->init(tolerance, maxiter);
		}
		if (params.contains("pos_vbd_iter"))
			m_pos_vbd_iterations = params["pos_vbd_iter"].get<unsigned int>();

		if (params.contains("amp_vbd_iter"))
			m_amp_vbd_iterations = params["amp_vbd_iter"].get<unsigned int>();

		// move all vertices that are inside the obstacles
		initMeshes(dynamicObjects, obstacles);
		initSelfCollision(dynamicObjects[0]->getTriangles(), dynamicObjects[0]->getRestPositions(), selfCollisionData);
		initVBD(dynamicObjects);
		initWrinkleMesh(dynamicObjects, true, true);
		optimize(dynamicObjects, true, true, false);
	}
	void WrinklesVBD::initMeshes(std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obstacles)
	{
		for (auto&& object : dynamicObjects)
		{
			Eigen::VectorXi cv = object->getClampedVertices();
			int nverts = cv.size();
			VectorN X, Xold, V;
			X.resize(3 * nverts);
			V.resize(3 * nverts);
			V.setZero();
			Xold.resize(3 * nverts);
			for (int i = 0; i < nverts; i++)
			{
				X.segment<3>(3 * i) = object->getPositions().segment<3>(3 * cv(i));
			}
			Xold = X;
			for (auto&& obs : obstacles)
				obs->project(X, V, Xold, 1, 0, NULL);
			for (int i = 0; i < nverts; i++)
			{
				object->getPositions().segment<3>(3 * cv(i)) = X.segment<3>(3 * i);
			}
		}
	}
	void WrinklesVBD::initVBD(std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects)
	{
		for (auto&& object : dynamicObjects)
		{
			object->requireVertexColoringInfo();
		}
	}
	void WrinklesVBD::addExternalForces(std::vector<std::unique_ptr<DynamicsObject>>& objects, Real timestep)
	{
		// cylinder twist 16x16 new
		// working for cylinder twist example
		Real theta = 2.5 * timestep;
		for (int i = 0; i <= 8; i++)
		{
			// rotate about z
			Matrix33 rotation;
			rotation << cos(theta), -sin(theta), 0,
			    sin(theta), cos(theta), 0,
			    0, 0, 1;
			Vector3 pos = objects.at(0)->getPositions()(Eigen::seq(3 * i, 3 * i + 2));
			objects.at(0)->getPositions()(Eigen::seq(3 * i, 3 * i + 2)) = rotation * pos;
			objects.at(0)->getVelocities()(Eigen::seq(3 * i, 3 * i + 2)) = (rotation * pos - pos) / timestep;
		}
		for (int i = 85; i <= 113; i += 4)
		{
			// rotate about z
			Matrix33 rotation;
			rotation << cos(theta), -sin(theta), 0,
			    sin(theta), cos(theta), 0,
			    0, 0, 1;
			Vector3 pos = objects.at(0)->getPositions()(Eigen::seq(3 * i, 3 * i + 2));
			objects.at(0)->getPositions()(Eigen::seq(3 * i, 3 * i + 2)) = rotation * pos;
			objects.at(0)->getVelocities()(Eigen::seq(3 * i, 3 * i + 2)) = (rotation * pos - pos) / timestep;
		}
		// for highres
		// Real theta = 0.05 * timestep;
		// for (int i = 1; i <= 8; i++)
		// {
		// 	// rotate about z
		// 	Matrix33 rotation;
		// 	rotation << cos(theta), -sin(theta), 0,
		// 	    sin(theta), cos(theta), 0,
		// 	    0, 0, 1;
		// 	Vector3 pos = objects.at(0)->getPositions()(Eigen::seq(3 * i, 3 * i + 2));
		// 	objects.at(0)->getPositions()(Eigen::seq(3 * i, 3 * i + 2)) = rotation * pos;
		// }
		// for (int i = 85; i <= 113; i += 4)
		// {
		// 	// rotate about z
		// 	Matrix33 rotation;
		// 	rotation << cos(theta), -sin(theta), 0,
		// 	    sin(theta), cos(theta), 0,
		// 	    0, 0, 1;
		// 	Vector3 pos = objects.at(0)->getPositions()(Eigen::seq(3 * i, 3 * i + 2));
		// 	objects.at(0)->getPositions()(Eigen::seq(3 * i, 3 * i + 2)) = rotation * pos;
		// }
		// for compression test
		// Vector3 oldPos, newPos;
		// // diagonal compression
		// objects.at(0)->getPositions()(Eigen::seq(0, 2)) += 0.01 * Vector3(1, 0, 1);
		// objects.at(0)->getPositions()(Eigen::seq(9, 11)) -= 0.01 * Vector3(1, 0, 1);
		// shearing
		// oldPos = objects.at(0)->getPositions()(Eigen::seq(0, 2));
		// objects.at(0)->getPositions()(Eigen::seq(0, 2)) += 0.05 * Vector3(1, 0, 0);
		// newPos = objects.at(0)->getPositions()(Eigen::seq(0, 2));
		// oldPos = objects.at(0)->getPositions()(Eigen::seq(9, 11));
		// objects.at(0)->getPositions()(Eigen::seq(9, 11)) -= 0.0001 * Vector3(1, 0, 0);
		// newPos = objects.at(0)->getPositions()(Eigen::seq(9, 11));

		/**
		 * Add external forces that define the simulation
		 */
		// for compression test
		// objects.at(0)->getPositions()(Eigen::seq(0, 2)) += 0.1 * Vector3(1, 0, 1);
		// objects.at(0)->getPositions()(Eigen::seq(9, 11)) -= 0.1 * Vector3(1, 0, 1);
		/* for (auto&& object : objects)
		{
		    // object->getPositions()(Eigen::seq(0, 2)) += Vector3({ 0.01, 0, 0 });
		    // object->getPositions()(Eigen::seq(6, 8)) += Vector3({ 0.01, 0, 0 });
		    for (int i = 0; i < 9; i++)
		    {
		        object->getPositions()(Eigen::seq(i * 9 * 3, i * 9 * 3 + 2)) += Vector3({ 0.001, 0, 0 });
		    }
		} */
		for (auto&& object : objects)
		{
			/* Real theta = 0.5;
			Matrix33 rot;

			rot << cos(theta), -sin(theta), 0,
			    sin(theta), cos(theta), 0, 0, 0, 1;
			for (int vtex = 14; vtex <= 28; vtex++)
			{
			    int vid = vtex;
			    argus::Vector3 vec = object->getPositions()(Eigen::seq(3 * vid, 3 * vid + 2));
			    object->getPositions()(Eigen::seq(3 * vid, 3 * vid + 2)) = rot * vec;
			} */
			//
			// // top vertices
			// Real theta = 0.01;
			// Matrix33 rot;
			//
			// // rot << cos(theta), 0, -sin(theta),
			// //     0, 1, 0, sin(theta), 0, cos(theta);
			// rot << 1, 0, 0,
			//     0, cos(theta), -sin(theta),
			//     0, sin(theta), cos(theta);
			// int v1 = 9;
			// int v2 = 0;
			// for (int vtex = 0; vtex <= 9; vtex++)
			// {
			// 	int vid = vtex * 10 + v1;
			// 	argus::Vector3 vec = object->getPositions()(Eigen::seq(3 * vid, 3 * vid + 2));
			// 	object->getPositions()(Eigen::seq(3 * vid, 3 * vid + 2)) = rot * vec;
			// 	object->getVelocities()(Eigen::seq(3 * vid, 3 * vid + 2)) = rot * vec - vec;
			// 	/* Matrix33 newrot;
			// 	newrot << 1, 0, 0,
			// 	    0, cos(-theta), -sin(-theta),
			// 	    0, sin(-theta), cos(-theta);
			//
			// 	vid = vtex * 10 + v2;
			// 	vec = object->getPositions()(Eigen::seq(3 * vid, 3 * vid + 2));
			// 	object->getPositions()(Eigen::seq(3 * vid, 3 * vid + 2)) = (newrot * vec); */
			// }
			// /* for (int vtex = 0; vtex <= 39; vtex++)
			//  * cylinder 16_40
			// {
			// v1 = 15;
			//     int vid = vtex * 16 + v1;
			//     argus::Vector3 vec = object->getPositions()(Eigen::seq(3 * vid, 3 * vid + 2));
			//     object->getPositions()(Eigen::seq(3 * vid, 3 * vid + 2)) = rot * vec;
			// } */
		}
	}

	void WrinklesVBD::applyConstraints(Real timestep, std::vector<std::unique_ptr<DynamicsObject>>& objects)
	{
		for (auto&& object : objects)
		{
			object->applyDynamicConstraints(timestep);
			object->applyTopologicalConstraints();
		}
	}

	/* void WrinklesVBD::step(Real timestep, std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs, nlohmann::json& state, unsigned int frame)
	{
	    // if (frame < 500)
	        addExternalForces(dynamicObjects, timestep);

	    // update amps and positions
	    vbdUpdateAmpPosCombined(timestep, dynamicObjects, obs, m_stats, frame);

	    // init dphis
	    auto t5 = argus::util::captureCurrentTime();
	    initWrinkleMesh(dynamicObjects, false, true);
	    auto t6 = argus::util::captureCurrentTime();
	    m_stats["init"] = m_stats.value("init", 0LL) + std::chrono::duration_cast<std::chrono::milliseconds>(t6 - t5).count();

	    // update dphis
	    auto t7 = argus::util::captureCurrentTime();
	    dynamicObjects[0]->cacheParameters(std::array<bool, 5>({ true, true, true, false, false }));
	    optimize(dynamicObjects, false, true, false);
	    dynamicObjects[0]->clearCache();

	    auto t8 = argus::util::captureCurrentTime();
	    m_stats["optimizeDphis"] = m_stats.value("optimizeDphis", 0LL) + std::chrono::duration_cast<std::chrono::milliseconds>(t8 - t7).count();

	    // apply constraints
	    auto t9 = argus::util::captureCurrentTime();
	    applyConstraints(timestep, dynamicObjects);
	    auto t10 = argus::util::captureCurrentTime();
	    m_stats["constraint"] = m_stats.value("constraint", 0LL) + std::chrono::duration_cast<std::chrono::milliseconds>(t10 - t9).count();
	    // for obstacle sequences
	    int nobs = obs.size();
	    for (int i = 0; i < nobs; i++)
	    {
	        obs.at(i)->updateFrame();
	    }
	    dynamicObjects[0]->clearCache();
	} */

	void WrinklesVBD::step(Real timestep, std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs, nlohmann::json& state, unsigned int frame)
	{
		// if (frame < 40)
		// 	addExternalForces(dynamicObjects, timestep);
		if (m_ampPosUpdate == "separate")
		{ // update amps
			auto t1 = argus::util::captureCurrentTime();
			dynamicObjects[0]->cacheParameters(std::array<bool, 5>({ true, true, false, true, false }));
			vbdUpdateAmps(timestep, dynamicObjects, obs);
			// vbdUpdateAmpsDphis(timestep, dynamicObjects, obs);
			dynamicObjects[0]->clearCache();
			auto t2 = argus::util::captureCurrentTime();
			m_stats["vbdUpdateAmps"] = m_stats.value("vbdUpdateAmps", 0LL) + std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

			// update pos
			auto t3 = argus::util::captureCurrentTime();
			// dynamicObjects[0]->cacheParameters(std::array<bool, 5>({ false, false, true, true, true }));
			updatePos(timestep, dynamicObjects, obs, state, frame);
			// dynamicObjects[0]->clearCache();
			auto t4 = argus::util::captureCurrentTime();
			m_stats["vbdUpdatePos"] = m_stats.value("vbdUpdatePos", 0LL) + std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count();
		}
		else if (m_ampPosUpdate == "combined")
		{
			auto t3 = argus::util::captureCurrentTime();
			vbdUpdateAmpPosCombined(timestep, dynamicObjects, obs, m_stats, frame);
			auto t4 = argus::util::captureCurrentTime();
			m_stats["vbdUpdatePos"] = m_stats.value("vbdUpdatePos", 0LL) + std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count();
		}
		else
		{
			throw "Error in wrinkles-vbd step! Incorrect amp pos update sequence passed!";
		}
		// init dphis
		auto t5 = argus::util::captureCurrentTime();
		initWrinkleMesh(dynamicObjects, false, true);
		auto t6 = argus::util::captureCurrentTime();
		m_stats["init"] = m_stats.value("init", 0LL) + std::chrono::duration_cast<std::chrono::milliseconds>(t6 - t5).count();

		// update dphis
		auto t7 = argus::util::captureCurrentTime();
		dynamicObjects[0]->cacheParameters(std::array<bool, 5>({ true, true, true, false, false }));
		optimize(dynamicObjects, false, true, false);
		dynamicObjects[0]->clearCache();
		auto t8 = argus::util::captureCurrentTime();
		m_stats["optimizeDphis"] = m_stats.value("optimizeDphis", 0LL) + std::chrono::duration_cast<std::chrono::milliseconds>(t8 - t7).count();
		// apply constraints
		auto t9 = argus::util::captureCurrentTime();
		applyConstraints(timestep, dynamicObjects);
		auto t10 = argus::util::captureCurrentTime();
		m_stats["constraint"] = m_stats.value("constraint", 0LL) + std::chrono::duration_cast<std::chrono::milliseconds>(t10 - t9).count();
		// for obstacle sequences
		int nobs = obs.size();
		for (int i = 0; i < nobs; i++)
		{
			obs.at(i)->updateFrame();
		}
		dynamicObjects[0]->clearCache();
	}

	void WrinklesVBD::updatePos(Real timestep, std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs, nlohmann::json& state, unsigned int frame)
	{
		vbdUpdatePos(timestep, dynamicObjects, obs, state, frame);
	}

	void WrinklesVBD::optimizeDphisParallel(std::vector<std::unique_ptr<DynamicsObject>>& objects, const std::vector<int>& faceids)
	{
		if (m_optimizer)
			for (auto&& object : objects)
			{
				unsigned int nfaces = faceids.size();

				for (int fi = 0; fi < nfaces; fi++)
				{
					int faceid = faceids.at(fi);
					std::function<Real(VectorN & X, VectorN & deriv)> func = [&object, faceid](VectorN& x, VectorN& grad)
					{
						if (x.size() != 2)
							throw "Wrong input passed to energy function for dphi optimization!";
						// write energy equation
						// save for restoring
						Vector2 bak = object->getDphis().col(faceid);
						// set
						object->getDphis().col(faceid) = Vector2(x);
						// query energies
						VectorN tderiv;
						Real e = object->getWrinkleShellEnergiesPerFace(faceid, tderiv);
						grad.resize(2);
						grad = tderiv(Eigen::seq(3, 4));
						// reset
						object->getDphis().col(faceid) = bak;
						return e;
					};
					/* // derivative test
					VectorN ip = VectorN::Random(2) * 9998;
					VectorN dx = VectorN::Random(2).normalized();
					Real h = 0.9;
					Real e0, e1, e2;
					VectorN d0, d1, d2, X1, X2;
					e0 = func(ip, d0);
					MatrixNN errors(nfaces, 10);
					errors.setZero();
					for (int i = 0; i < 10; i++)
					{
					    X1 = ip + 0.5 * h * dx;
					    e1 = func(X1, d1);
					    X2 = ip - 0.5 * h * dx;
					    e2 = func(X2, d2);
					    Real error = std::fabs(((e1 - e2) / h) - dx.dot(d0));
					    std::cout << "Error " << error << std::endl;
					    errors(faceid, i) = error;
					    h *= 0.1;
					} */
					VectorN X = object->getDphis().col(faceid);

					try
					{
						m_optimizer->solve(func, X);
						// update dphi
						object->getDphis().col(faceid) = Vector2(X);
					}
					catch (const char* err)
					{
						std::cerr << "Error while optimizing for dphis! " << err << std::endl;
					}
					catch (...)
					{
						std::cerr << "Other error while optimizing for dphis! " << std::endl;
					}
				}
			}
	}

	void WrinklesVBD::optimizeDphisParallel(std::vector<std::unique_ptr<DynamicsObject>>& objects)
	{
		if (m_optimizer)
			for (auto&& object : objects)
			{
				unsigned int nfaces = object->getFaceCount();

#ifdef USE_OMP
#pragma omp parallel for schedule(static)
#endif
				for (int faceid = 0; faceid < nfaces; faceid++)
				{
					std::function<Real(VectorN & X, VectorN & deriv)> func = [&object, faceid](VectorN& x, VectorN& grad)
					{
						if (x.size() != 2)
							throw "Wrong input passed to energy function for dphi optimization!";
						// write energy equation
						// save for restoring
						Vector2 bak = object->getDphis().col(faceid);
						// set
						object->getDphis().col(faceid) = Vector2(x);
						// query energies
						VectorN tderiv;
						Real e = object->getWrinkleShellEnergiesPerFace(faceid, tderiv);
						grad.resize(2);
						grad = tderiv(Eigen::seq(3, 4));
						// reset
						object->getDphis().col(faceid) = bak;
						return e;
					};
					/* // derivative test
					VectorN ip = VectorN::Random(2) * 9998;
					VectorN dx = VectorN::Random(2).normalized();
					Real h = 0.9;
					Real e0, e1, e2;
					VectorN d0, d1, d2, X1, X2;
					e0 = func(ip, d0);
					MatrixNN errors(nfaces, 10);
					errors.setZero();
					for (int i = 0; i < 10; i++)
					{
					    X1 = ip + 0.5 * h * dx;
					    e1 = func(X1, d1);
					    X2 = ip - 0.5 * h * dx;
					    e2 = func(X2, d2);
					    Real error = std::fabs(((e1 - e2) / h) - dx.dot(d0));
					    std::cout << "Error " << error << std::endl;
					    errors(faceid, i) = error;
					    h *= 0.1;
					} */
					VectorN X = object->getDphis().col(faceid);

					try
					{
						m_optimizer->solve(func, X);
						// update dphi
						object->getDphis().col(faceid) = Vector2(X);
					}
					catch (const char* err)
					{
						std::cerr << "Error while optimizing for dphis! " << err << std::endl;
					}
					catch (...)
					{
						std::cerr << "Other error while optimizing for dphis! " << std::endl;
					}
				}
			}
	}
	void WrinklesVBD::optimize(std::vector<std::unique_ptr<DynamicsObject>>& objects, bool amps, bool dphis, bool positions)
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
		/* Real getDphiEnergyPerFace(argus::VectorN X, int faceid, Sheet& sheet, VectorN* deriv = NULL, argus::MatrixNN* hess = NULL)
		    {
		        argus::Vector2 initDphi = sheet.getDphis1VectorPerFace().col(faceid);
		        sheet.getDphis1VectorPerFace().col(faceid) = argus::Vector2(X);
		        argus::Real energy = sheet.getMaterial().getWrinkleShellEnergiesPerFace(sheet, faceid, deriv, hess);
		        sheet.getDphis1VectorPerFace().col(faceid) = initDphi;
		        return energy;
		    } */
		if (dphis)
		{
			// optimizeDphisParallel(objects);
			for (auto&& object : objects)
			{
				object->optimizeParallelEnergies();
			}
			// TODO: resuse it
			/* if (!m_optimizer)
			    throw "Error optimizer is not initialized before calling optimize!\n";
			for (auto&& object : objects)
			{
			    unsigned int vertexCount, faceCount;
			    vertexCount = object->getVertexCount();
			    faceCount = object->getFaceCount();
			    Matrix2N newDphis;
			    newDphis = object->getDphis();
			    Real totale = 0;

			    // parallel solve
			    // #ifdef USE_OMP
			    // #pragma omp parallel for
			    // #endif
			    for (unsigned int i = 0; i < faceCount; i++)
			    {
			        std::function<Real(VectorN&, VectorN&)> func = [&object, i, vertexCount, faceCount](VectorN& X, VectorN& deriv) -> Real
			        {
			            VectorN tempDeriv;
			            Vector2 initDphi = object->getDphis().col(i);
			            object->getDphis().col(i) = Vector2(X);
			            Real e = object->getWrinkleShellEnergiesPerFace(i, tempDeriv, NULL);
			            deriv.resize(2);
			            deriv = tempDeriv(Eigen::seq(3, 4));
			            object->getDphis().col(i) = initDphi;
			            return e;
			        };
			        VectorN X = object->getDphis().col(i);
			        if (m_optimizer)
			        {
			            try
			            {
			                Real e = m_optimizer->solve(func, X);
			                totale += e;
			                object->getDphis().col(i) = Vector2(X);
			                // newDphis.col(i) = Vector2(X);
			            }
			            catch (const char* err)
			            {
			                std::cerr << "Error while optimizing for dphis! " << err << std::endl;
			            }
			            catch (...)
			            {
			                std::cerr << "Other error while optimizing for dphis! " << std::endl;
			            }
			        }
			    }
			    // std::cout << object->getDphis().transpose();
			    // std::cout << totale << std::endl;
			    // exit(0);

			    // total optimization
			    // std::function<Real(VectorN & X, VectorN & deriv)> func = [&object, vertexCount, faceCount](VectorN& X, VectorN& deriv) -> Real
			    // {
			    //     VectorN ip, fderiv;
			    //     ip.resize(4 * vertexCount + 2 * faceCount);
			    //     fderiv.resize(4 * vertexCount + 2 * faceCount);
			    //     ip(Eigen::seq(0, vertexCount - 1)) = object->getAmpVec();
			    //     ip(Eigen::seq(vertexCount, vertexCount + 2 * faceCount - 1)) = X;
			    //     ip(Eigen::seq(vertexCount + 2 * faceCount, 4 * vertexCount + 2 * faceCount - 1)) = object->getPositions().reshaped();
			    //     Real e = object->getEnergies(ip, &fderiv);
			    //     deriv = fderiv(Eigen::seq(vertexCount, vertexCount + 2 * faceCount - 1));
			    //     return e;
			    // };
			    //
			    // VectorN X = object->getDphis().reshaped();
			    // if (m_optimizer)
			    //     m_optimizer->solve(func, X);
			    // else
			    //     throw "Error optimizer is not initialized before calling optimize!\n";
			    // object->getDphis() = X.reshaped(2, faceCount);

			} */
		}

		if (positions)
			throw "Static optimizer not implemented for positions!\n";
	}

	VectorN WrinklesVBD::getAdaptiveAccelerationForVBD(DynamicsObject& object, VectorN& var, bool amp, bool pos)
	{
		unsigned int vertexCount = object.getVertexCount();
		unsigned int faceCount = object.getFaceCount();
		VectorN invInertia = object.getInvInertia();
		VectorN inertia = object.getInertia();
		const SheetInfo& vbdInfo = object.getSheetInfo();
		if (pos && amp)
		{
			MatrixNN massIMat(4 * vertexCount, 4 * vertexCount);
			massIMat.setZero();
			for (int i = 0; i < vertexCount; i++)
			{
				if (inertia(i) != 0)
					massIMat(i, i) = 1 / inertia(i);
				// invInertia(3 * i);
			}
			for (int i = vertexCount; i < 4 * vertexCount; i++)
			{
				massIMat(i, i) = invInertia(i - vertexCount);
			}
			VectorN acc;
			VectorN axt = vbdInfo.axt;
			VectorN aat = vbdInfo.aat;
			VectorN at;
			at.resize(4 * vertexCount);
			at(Eigen::seq(0, vertexCount - 1)) = aat;
			at(Eigen::seq(vertexCount, 4 * vertexCount - 1)) = axt;
			VectorN f, forces;
			f.resize(4 * vertexCount + 2 * faceCount);
			f.setZero();
			VectorN ip;
			ip.resize(4 * vertexCount + 2 * faceCount);
			ip(Eigen::seq(0, vertexCount - 1)) = var(Eigen::seq(0, vertexCount - 1));
			ip(Eigen::seq(vertexCount, vertexCount + 2 * faceCount - 1)) = object.getDphis().reshaped();
			ip(Eigen::seq(vertexCount + 2 * faceCount, 4 * vertexCount + 2 * faceCount - 1)) = var(Eigen::seq(vertexCount, 4 * vertexCount - 1));
			object.getEnergies(ip, &f, NULL, true);
			f = -f;
			forces.resize(4 * vertexCount);
			forces(Eigen::seq(0, vertexCount - 1)) = 2 * f(Eigen::seq(0, vertexCount - 1));
			forces(Eigen::seq(vertexCount, 4 * vertexCount - 1)) = f(Eigen::seq(vertexCount + 2 * faceCount, 4 * vertexCount + 2 * faceCount - 1));
			VectorN aext = massIMat * forces;
			// return non zero acceleration for cases with zero initialial velocities and accelerations
			if (at.norm() < 1e-4)
			{
				return aext;
			}
			Real atilda, atext;
			atext = at.dot(aext.normalized());
			if (atext > aext.norm())
			{
				return aext;
			}
			else if (atext < 0)
			{
				return VectorN::Zero(vertexCount * 4);
			}
			else
			{
				return (atext / aext.norm()) * aext;
			}
		}
		else if (amp)
		{
			MatrixNN massIMat(vertexCount, vertexCount);
			massIMat.setZero();
			for (int i = 0; i < vertexCount; i++)
			{
				massIMat(i, i) = invInertia(3 * i);
			}
			VectorN acc;
			VectorN aat = vbdInfo.aat;
			VectorN at;
			at = aat;
			VectorN f;
			f.resize(4 * vertexCount + 2 * faceCount);
			f.setZero();
			VectorN ip;
			ip.resize(4 * vertexCount + 2 * faceCount);
			ip(Eigen::seq(0, vertexCount - 1)) = var;
			ip(Eigen::seq(vertexCount, vertexCount + 2 * faceCount - 1)) = object.getDphis().reshaped();
			ip(Eigen::seq(vertexCount + 2 * faceCount, 4 * vertexCount + 2 * faceCount - 1)) = object.getPositions().reshaped();
			object.getEnergies(ip, &f, NULL, true);
			f = -f;
			VectorN aext = massIMat * f(Eigen::seq(0, vertexCount - 1));
			// return non zero acceleration for cases with zero initialial velocities and accelerations
			if (at.norm() < 1e-4)
			{
				return aext;
			}
			Real atilda, atext;
			atext = at.dot(aext.normalized());
			if (atext > aext.norm())
			{
				return aext;
			}
			else if (atext < 0)
			{
				return VectorN::Zero(vertexCount * 4);
			}
			else
			{
				return (atext / aext.norm()) * aext;
			}
		}
		else if (pos)
		{
			VectorN acc;
			VectorN at = vbdInfo.axt;
			VectorN aext, f;
			f.resize(4 * vertexCount + 2 * faceCount);
			f.setZero();
			VectorN ip;
			ip.resize(4 * vertexCount + 2 * faceCount);
			ip(Eigen::seq(0, vertexCount - 1)) = object.getAmpVec();
			ip(Eigen::seq(vertexCount, vertexCount + 2 * faceCount - 1)) = object.getDphis().reshaped();
			ip(Eigen::seq(vertexCount + 2 * faceCount, 4 * vertexCount + 2 * faceCount - 1)) = var;
			object.getEnergies(ip, &f, NULL, true);
			f = -f;
			aext = invInertia.asDiagonal() * f(Eigen::seq(vertexCount + 2 * faceCount, vertexCount * 4 + 2 * faceCount - 1));
			// return non zero acceleration for cases with zero initial velocities and accelerations
			if (at.norm() < 1e-4)
				return aext;
			Real atilda, atext;
			atext = at.dot(aext.normalized());
			if (atext > aext.norm())
				return aext;
			else if (atext < 0)
				return VectorN::Zero(vertexCount * 3);
			else
				return (atext / aext.norm()) * aext;
		}
		else
		{
			throw "Calling adaptiveAccelerationForVBD without any dof. Not implemented! Exiting.\n";
		}
	}
	void WrinklesVBD::vbdUpdateAmps(Real timestep, std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs)
	{
		for (auto&& object : dynamicObjects)
		{
			const SheetInfo& vbdInfo = object->getSheetInfo();
			// get mass information
			VectorN& mass = object->getInertia();
			VectorN invInertia = object->getInvInertia();
			// update positions from external forces and inertia
			// adaptive intertia and acceleration
			VectorN xt = object->getAmpVec();
			VectorN vt = object->getAmpVelocities();

			// external forces
			VectorN aext = getAdaptiveAccelerationForVBD(*object, xt, true, false);
			VectorN y = timestep * vt + timestep * timestep * aext;
			y = xt + y;

			object->clampDofs(y);
			VectorN x = y;
			int niters = m_amp_vbd_iterations;
			int ncolors = vbdInfo.coloredSets.rows();
#ifdef ARGUS_CHECKPOINT
			if (argus::simConf.contains("VBD Iterations for Amps") == false)
				argus::simConf["VBD Iterations for Amps"] = niters;
#endif

			for (int iter = 0; iter < niters; iter++)
			{
				// Inside the loop
				for (int color = 0; color < ncolors; color++)
				{
					argus::VectorNi vset;

					vset = vbdInfo.coloredSets.row(color);
					// Eigen::internal::set_is_malloc_allowed(false); // to check if dynamic allocation is happening
#ifdef USE_OMP
#pragma omp parallel for schedule(static)
#endif
					for (int vtex = 0; vtex < vset.size(); vtex++)
					{
						int vid = vset(vtex);
						if (vid == -1)
							continue;
						if (invInertia(3 * vid) == 0)
							continue;
						argus::Real xi, yi;
						xi = x(vid);
						yi = y(vid);
						argus::Real deriv;
						argus::Real hess;
						argus::Real energy = object->getVertexEnergy1(vid, &deriv, &hess);
						argus::Real mvtex = mass(vid);
						argus::Real fi = -(mvtex / (timestep * timestep)) * (xi - yi) - 2 * deriv; // here a factor of 2 due to energy
						argus::Real Hi = (mvtex / (timestep * timestep)) + 2 * hess; // here a factore of 2 due to energy
						// do not change if Hi is singular
						if (std::fabs(Hi) < 1e-6)
							continue;
						if (Hi < 0)
							Hi = -Hi;
						argus::Real dxi
						    = fi / Hi; // a multiplier of 2 as amplitude equation involves double increase in acceleration
						// uses dlib line search
						// function for line search
						auto func = [&object, &yi, &timestep, &mvtex, &xi, &dxi, &x, vid](Real a, Real* deriv = NULL)
						{
							Real abak = object->getAmpVec()(vid);
							object->getAmpVec()(vid) = xi + a * dxi;
							Real e = 0;
							if (deriv)
								e = object->getVertexEnergy1(vid, deriv);
							else
								e = object->getVertexEnergy1(vid, NULL, NULL);
							object->getAmpVec()(vid) = abak;
							e += (mvtex / (2.0 * timestep * timestep)) * (xi + a * dxi - yi) * (xi + a * dxi - yi);
							if (deriv)
							{
								(*deriv) += (mvtex / (timestep * timestep)) * (xi + a * dxi - yi);
							}
							return e;
						};
						// TODO: testing the effect of line serch - causes amp to change suddenly in a time step?
						// argus::Real alpha = dlib::backtracking_line_search(func, func(0), deriv * dxi, 1, 0.5, 5);
						// if (alpha > 1e-4)
						// {
						// 	x(vid) = xi + dxi;
						// }
						if (func(1) < energy + (mvtex / (2.0 * timestep * timestep)) * (xi - yi) * (xi - yi))
						{
							x(vid) = xi + dxi;
						}
					}

					// Eigen::internal::set_is_malloc_allowed(true);
					// update positions
					object->getAmpVec().noalias() = x;
				}

				x.noalias() = object->getAmpVec();
				// clamp
				// object->clampDofs(x);
				// update velocities
				argus::VectorN vnew = (x - xt) / timestep;
				// damping
				vnew *= (1 - m_kdamp_amps);
				object->getAmpVelocities().noalias() = vnew;
				object->getAmpVec().noalias() = xt + vnew * timestep;
			}
		}
	}

	void WrinklesVBD::vbdUpdateAmpsDphis(Real timestep, std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs)
	{
		for (auto&& object : dynamicObjects)
		{
			const std::vector<std::vector<int>> v2f = object->getVertexFaceAdjacencyList();
			const SheetInfo& vbdInfo = object->getSheetInfo();
			// get mass information
			VectorN& mass = object->getInertia();
			VectorN invInertia = object->getInvInertia();
			// update positions from external forces and inertia
			// adaptive intertia and acceleration
			VectorN xt = object->getAmpVec();
			VectorN vt = object->getAmpVelocities();

			// external forces
			VectorN aext = getAdaptiveAccelerationForVBD(*object, xt, true, false);
			VectorN y = timestep * vt + timestep * timestep * aext;
			y = xt + y;

			object->clampDofs(y);
			VectorN x = y;
			int niters = m_amp_vbd_iterations;
			int ncolors = vbdInfo.coloredSets.rows();
#ifdef ARGUS_CHECKPOINT
			if (argus::simConf.contains("VBD Iterations for Amps") == false)
				argus::simConf["VBD Iterations for Amps"] = niters;
#endif

			for (int iter = 0; iter < niters; iter++)
			{
				// Inside the loop
				for (int color = 0; color < ncolors; color++)
				{
					argus::VectorNi vset;

					vset = vbdInfo.coloredSets.row(color);
					// Eigen::internal::set_is_malloc_allowed(false); // to check if dynamic allocation is happening
#ifdef USE_OMP
#pragma omp parallel for schedule(static)
#endif
					for (int vtex = 0; vtex < vset.size(); vtex++)
					{
						int vid = vset(vtex);
						if (vid == -1)
							continue;
						if (invInertia(3 * vid) == 0)
							continue;
						argus::Real xi, yi;
						xi = x(vid);
						yi = y(vid);
						argus::Real deriv;
						argus::Real hess;
						argus::Real energy = object->getVertexEnergy1(vid, &deriv, &hess);
						argus::Real mvtex = mass(vid);
						argus::Real fi = -(mvtex / (timestep * timestep)) * (xi - yi) - 2 * deriv; // here a factor of 2 due to energy
						argus::Real Hi = (mvtex / (timestep * timestep)) + 2 * hess; // here a factore of 2 due to energy
						// do not change if Hi is singular
						if (std::fabs(Hi) < 1e-6)
							continue;
						if (Hi < 0)
							Hi = -Hi;
						argus::Real dxi
						    = fi / Hi; // a multiplier of 2 as amplitude equation involves double increase in acceleration
						// uses dlib line search
						// function for line search
						auto func = [&object, &xi, &dxi, &x, vid](Real a, Real* deriv = NULL)
						{
							object->getAmpVec()(vid) = xi + a * dxi;
							Real e = 0;
							if (deriv)
								e = object->getVertexEnergy1(vid, deriv);
							else
								e = object->getVertexEnergy1(vid, NULL, NULL);
							object->getAmpVec()(vid) = x(vid);
							return e;
						};
						argus::Real alpha = dlib::backtracking_line_search(func, energy, deriv * dxi, 1, 0.5, 5);
						if (alpha > 1e-4)
						{
							x(vid) = xi + dxi;
							object->getAmpVec()(vid) = x(vid);
							optimizeDphisParallel(dynamicObjects, v2f[vid]);
						}
					}

					// Eigen::internal::set_is_malloc_allowed(true);
					// update positions
					object->getAmpVec().noalias() = x;
				}

				x.noalias() = object->getAmpVec();
				// clamp
				object->clampDofs(x);
				// update velocities
				argus::VectorN vnew = (x - xt) / timestep;
				// damping
				vnew *= (1 - m_kdamp_amps);
				object->getAmpVelocities().noalias() = vnew;
				object->getAmpVec().noalias() = xt + vnew * timestep;
			}
		}
	}

	void WrinklesVBD::vbdUpdatePos(Real timestep, std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs, nlohmann::json& state, unsigned int frame)
	{

		// TODO: Remember to save old positions before calling this function to resolve collisions further.
		for (auto&& object : dynamicObjects)
		{

			const SheetInfo& vbdInfo = object->getSheetInfo();
			// get mass information
			argus::VectorN mass = object->getInertia();
			VectorN invInertia = object->getInvInertia();
			// update positions from external forces and inertia
			// adaptive intertia and acceleration
			argus::VectorN xt = object->getPositions();
			argus::VectorN vt = object->getVelocities();
			object->clampDofs(vt);
			// external forces
			argus::VectorN aext = getAdaptiveAccelerationForVBD(*object, xt, false, true);
			argus::VectorN y = timestep * vt + timestep * timestep * aext;
			// TODO: Decide whether to use or not?
			// wme->clampWrinkleDofs(y);
			y = xt + y;
			// get collision free updates
			// wme->getObstacleProjectedPositions(obs, y);
			// y =  getCollisionFreeState(obs, y, xt, timestep, NULL, frame);
			// y= getCollisionFreeState(obs, y, xt, object->getTriangles(), timestep, NULL, frame);
			getCollisionFreeState(obs, y, vt, xt, invInertia, object->getTriangles(), timestep, NULL, frame);
			// Real ccdth = selfCollisionParams.threshold;
			// // self collision substeppng
			// resolveSelfCollisionsSubstepped(
			//     object->getTriangles(),
			//     xt, // Position at the very start of the frame
			//     y, // Predicted position (modified to resolved position)
			//     vt, // Velocity (modified)
			//     selfCollisionParams.self_collision_maxsubsteps,
			//     invInertia,
			//     selfCollisionData,
			//     0.0001, // stiffness
			//     1, // damping
			//     timestep);

			argus::VectorN x = y;
			int nverts = object->getVertexCount();
			Eigen::MatrixXi F = object->getFaces();
			int niters = m_pos_vbd_iterations;
#ifdef ARGUS_CHECKPOINT
			if (argus::simConf.contains("VBD Iterations for Pos") == false)
				argus::simConf["VBD Iterations for Pos"] = niters;
#endif
			std::vector<int> timings;
			std::vector<unsigned int> collverts;
			for (int iter = 0; iter < niters; iter++)
			{
				int ncolors = vbdInfo.coloredSets.rows();
				// Inside the loop
				for (int color = 0; color < ncolors; color++)
				{
					argus::VectorNi set = vbdInfo.coloredSets.row(color);
					// Eigen::internal::set_is_malloc_allowed(false);
					// auto t1 = argus::util::captureCurrentTime();
#ifdef USE_OMP
#pragma omp parallel for schedule(static)
#endif
					for (int vtex = 0; vtex < set.rows(); vtex++)
					{
						int vid = set(vtex);
						if (vid == -1)
							continue;
						if (invInertia(3 * vid) == 0)
							continue;
						Vector3 xi = x(Eigen::seq(3 * vid, 3 * vid + 2));
						Vector3 yi = y(Eigen::seq(3 * vid, 3 * vid + 2));
						Vector3 deriv;
						Matrix33 hess;
						argus::Real energy = object->getVertexEnergy(vid, &deriv, &hess);
						Real mvtex = mass(vid);
						Vector3 fi = -(mvtex / (timestep * timestep)) * (xi - yi) - deriv;
						Matrix33 Hi = (mvtex / (timestep * timestep)) * argus::Matrix33::Identity() + hess;

						// do not change if Hi is singular
						if (std::fabs(Hi.determinant()) < 1e-5)
							continue;
						argus::Vector3 dxi = Hi.inverse() * fi;
						// project for clamping
						// wme->clampVertexWrinkleDofs(dxi, vid);
						argus::Real alpha = 1;
						// alpha = wme->backtrackingLineSearch(func, xi, dxi, alpha, 10);
						auto func = [&object, &xi, &dxi, vid, &x, &yi, &timestep, &mvtex](Real rate)
						{
							object->getPositions()(Eigen::seq(3 * vid, 3 * vid + 2)).noalias() = xi + rate * dxi;
							Real e = object->getVertexEnergy(vid, NULL, NULL);
							object->getPositions()(Eigen::seq(3 * vid, 3 * vid + 2)).noalias() = x(Eigen::seq(3 * vid, 3 * vid + 2));
							e += (mvtex / (2 * timestep * timestep)) * (xi + rate * dxi - yi).norm() * (xi + rate * dxi - yi).norm();
							return e;
						};

						Real f0 = energy + (mvtex / (2 * timestep * timestep)) * (xi - yi).norm() * (xi - yi).norm();
						// TODO: test if func(0) == energy
						Real f1 = func(1);
						alpha = 1;
						if (f1 < f0)
						{
							x(Eigen::seq(3 * vid, 3 * vid + 2)).noalias() = xi + alpha * dxi;
						}
					}

					// auto t2 = argus::util::captureCurrentTime();
					// auto duration = m_stats.value("init", 0LL) + std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
					// timings.push_back(duration);
					// Eigen::internal::set_is_malloc_allowed(true);
					// update positions
					object->getPositions().noalias() = x;
					object->applyDynamicConstraints(timestep);
					object->applyTopologicalConstraints();
				}
				x = object->getPositions();
				// update velocities
				argus::VectorN vnew = (x - xt) / timestep;
				if (selfCollisionParams.enable_self_collisions)
					resolveSelfCollisionsSubstepped(
					    object->getTriangles(),
					    xt, // Position at the very start of the frame
					    x, // Predicted position (modified to resolved position)
					    vnew, // Velocity (modified)
					    selfCollisionParams.self_collision_maxsubsteps,
					    invInertia,
					    selfCollisionData,
					    selfCollisionParams.stiffness, // stiffness
					    selfCollisionParams.damping, // damping
					    timestep);
				vnew
				    = (x - xt) / timestep;
				object->clampDofs(vnew);
				// ether damping
				vnew = vnew * (1 - m_kdamp_pos);
				x = xt + vnew * timestep;

				resolveCollisions(obs, x, xt, vnew, invInertia, object->getTriangles(), timestep, state, frame, collverts);
				// x = xt + vnew * timestep;

				object->getPositions() = x;
				object->getVelocities() = vnew;

				// this function updates positions as well as velocities
			}
			// for (int i = 0; i < timings.size(); i++)
			// {
			// 	std::cout << timings[i] << std::endl;
			// }

			// exit(0);
#ifdef ARGUS_CHECKPOINT
			// object->setCollisionVertices(collverts);
#endif
		}
	}

	void WrinklesVBD::getCollisionFreeState(std::vector<std::unique_ptr<Obstacle>>& obstacles, VectorN& X, VectorN& V, const VectorN& Xold, const VectorN& W, const std::vector<trimesh::triangle_t>& triangles, Real timestep, std::vector<int>* collVerts, unsigned int frame)
	{
		/**
		 * Project updated_pos to collision free states.
		 */
		int nobs = obstacles.size();
		for (int iter = 0; iter < m_collision_iterations; iter++)
			for (int i = 0; i < nobs; i++)
			{
				if (m_sdfCollider == "basic")
					obstacles.at(i)->project(X, V, Xold, timestep, frame, collVerts);
				else if (m_sdfCollider == "quadrature")
					obstacles.at(i)->project(timestep, X, V, Xold, W, triangles, frame, collVerts, &m_config);
				else
					throw "not implemented!";
			}
	}
	void WrinklesVBD::getCollisionFreeState(std::vector<std::unique_ptr<Obstacle>>& obstacles, VectorN& X, const VectorN& Xold, Real timestep, std::vector<int>* collVerts, unsigned int frame)
	{
		/**
		 * Project updated_pos to collision free states.
		 */
		throw "not implemented!";
		// int nobs = obstacles.size();
		// for (int i = 0; i < nobs; i++)
		// {
		// 	obstacles.at(i)->project(X, Xold, timestep, frame, collVerts);
		// }
	}
	void WrinklesVBD::getCollisionFreeState(std::vector<std::unique_ptr<Obstacle>>& obstacles, VectorN& X, const VectorN& Xold, const std::vector<trimesh::triangle_t>& triangles, Real timestep, std::vector<int>* collVerts, unsigned int frame)

	{
		/**
		 * Project updated_pos to collision free states.
		 */
		int nobs = obstacles.size();
		for (int i = 0; i < nobs; i++)
		{
			// X = obstacles.at(i)->project(X, Xold, triangles, timestep, frame, collVerts);
		}
	}
	void WrinklesVBD::resolveCollisions(std::vector<std::unique_ptr<Obstacle>>& obs, VectorN& updated_pos, const VectorN& oldPositions, VectorN& velocities, const VectorN& W, const std::vector<trimesh::triangle_t>& triangles, Real timestep, nlohmann::json& state, unsigned int frame, std::vector<unsigned int>& collVerts)
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
			VectorN inertialPos;
			std::vector<int> collisionVertices;

			// for friction
			if (m_sim_friction)
				inertialPos = updated_pos;
			for (int i = 0; i < coliters; i++)
			{
#ifdef ARGUS_CHECKPOINT
				getCollisionFreeState(obs, updated_pos, velocities, oldPositions, W, triangles, timestep, &collisionVertices, frame);
#else
				getCollisionFreeState(obs, updated_pos, velocities, oldPositions, W, triangles, timestep, NULL, frame);
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
			// old friction pipeline
			/* if (m_sim_friction && std::fabs(frictionK) > 1e-6)
			{

			    std::set<int> collVertsSet(collisionVertices.begin(), collisionVertices.end());
			    std::vector<unsigned int> collVertsVec(collVertsSet.begin(), collVertsSet.end());
			    unsigned int collVertsSize = collVertsVec.size();
#ifdef USE_OMP
#pragma omp parallel for
#endif
			    for (unsigned int i = 0; i < collVertsSize; i++)
			    {
			        unsigned int vid = collVertsVec.at(i);
			        Vector3 pos1, pos2, ncap;
			        pos1 = inertialPos(Eigen::seq(3 * vid, 3 * vid + 2));
			        pos2 = updated_pos(Eigen::seq(3 * vid, 3 * vid + 2));
			        ncap = (pos2 - pos1).normalized();
			        Vector3 vTan = velocities(Eigen::seq(3 * vid, 3 * vid + 2)) - velocities(Eigen::seq(3 * vid, 3 * vid + 2)).dot(ncap) * ncap;
			        velocities(Eigen::seq(3 * vid, 3 * vid + 2)) -= vTan;
			        velocities(Eigen::seq(3 * vid, 3 * vid + 2)) += (1 - frictionK) * vTan;
			        // Vector3 newVel = velocities(Eigen::seq(3 * vid, 3 * vid + 2));
			        // newVel -= (newVel.dot(Vector3 { 0, 0, 1 })) * Vector3({ 0, 0, 1 });
			        // for arm rotating example
			        // velocities(Eigen::seq(3 * vid, 3 * vid + 2))(0) = 0;
			        updated_pos(Eigen::seq(3 * vid, 3 * vid + 2)) = oldPositions(Eigen::seq(3 * vid, 3 * vid + 2)) + velocities(Eigen::seq(3 * vid, 3 * vid + 2)) * timestep;
			    }
			} */
		}
	}

	void WrinklesVBD::resolveCollisions(std::vector<std::unique_ptr<Obstacle>>& obs, VectorN& updated_pos, const VectorN& oldPositions, VectorN& velocities, Real timestep, nlohmann::json& state, unsigned int frame, std::vector<unsigned int>& collVerts)
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
			VectorN inertialPos;
			std::vector<int> collisionVertices;

			// for friction
			if (m_sim_friction)
				inertialPos = updated_pos;
			for (int i = 0; i < coliters; i++)
			{
#ifdef ARGUS_CHECKPOINT
				getCollisionFreeState(obs, updated_pos, oldPositions, timestep, &collisionVertices, frame);
#else
				getCollisionFreeState(obs, updated_pos, oldPositions, timestep, NULL, frame);
#endif
			}
			velocities = (updated_pos - oldPositions) / timestep;
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
			// old friction pipeline
			/* if (m_sim_friction && std::fabs(frictionK) > 1e-6)
			{

			    std::set<int> collVertsSet(collisionVertices.begin(), collisionVertices.end());
			    std::vector<unsigned int> collVertsVec(collVertsSet.begin(), collVertsSet.end());
			    unsigned int collVertsSize = collVertsVec.size();
#ifdef USE_OMP
#pragma omp parallel for
#endif
			    for (unsigned int i = 0; i < collVertsSize; i++)
			    {
			        unsigned int vid = collVertsVec.at(i);
			        Vector3 pos1, pos2, ncap;
			        pos1 = inertialPos(Eigen::seq(3 * vid, 3 * vid + 2));
			        pos2 = updated_pos(Eigen::seq(3 * vid, 3 * vid + 2));
			        ncap = (pos2 - pos1).normalized();
			        Vector3 vTan = velocities(Eigen::seq(3 * vid, 3 * vid + 2)) - velocities(Eigen::seq(3 * vid, 3 * vid + 2)).dot(ncap) * ncap;
			        velocities(Eigen::seq(3 * vid, 3 * vid + 2)) -= vTan;
			        velocities(Eigen::seq(3 * vid, 3 * vid + 2)) += (1 - frictionK) * vTan;
			        // Vector3 newVel = velocities(Eigen::seq(3 * vid, 3 * vid + 2));
			        // newVel -= (newVel.dot(Vector3 { 0, 0, 1 })) * Vector3({ 0, 0, 1 });
			        // for arm rotating example
			        // velocities(Eigen::seq(3 * vid, 3 * vid + 2))(0) = 0;
			        updated_pos(Eigen::seq(3 * vid, 3 * vid + 2)) = oldPositions(Eigen::seq(3 * vid, 3 * vid + 2)) + velocities(Eigen::seq(3 * vid, 3 * vid + 2)) * timestep;
			    }
			} */
		}
	}
	// void estimateAmpOmegaFromStrain(const std::vector<Eigen::Matrix2d>& abars,
	//     const Eigen::MatrixXd& curPos,
	//     const Eigen::MatrixXi& F,
	//     const std::set<int>& clampedVerts,
	//     double amplitudeEstimate,
	//     Eigen::VectorXd& amp,
	//     Eigen::MatrixXd& w,
	//     std::set<int>& tensionFaces)
	//
	void WrinklesVBD::initWrinkleMesh(std::vector<std::unique_ptr<DynamicsObject>>& objects, bool amps, bool dphis)
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

	void WrinklesVBD::vbdUpdate(Real timestep, std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs)
	{
		throw "Not implemented!!";
	}
	void WrinklesVBD::vbdUpdateAmpPosInSameIter(Real timestep, std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs, nlohmann::json& state, unsigned int frame)
	{
	}
	void WrinklesVBD::vbdUpdateAmpPosCombined(Real timestep, std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs, nlohmann::json& state, unsigned int frame)

	{
		// TODO: Remember to save old positions before calling this function to resolve collisions further.
		for (auto&& object : dynamicObjects)
		{
			int nverts = object->getVertexCount();
			int nfaces = object->getFaceCount();
			const SheetInfo& vbdInfo = object->getSheetInfo();
			// get mass information
			argus::VectorN mass = object->getInertia();
			VectorN invInertia = object->getInvInertia();
			// update positions from external forces and inertia
			// adaptive intertia and acceleration
			argus::VectorN xt(4 * nverts);
			xt(Eigen::seq(nverts, 4 * nverts - 1)) = object->getPositions();
			xt(Eigen::seq(0, nverts - 1)) = object->getAmpVec();
			argus::VectorN vt(4 * nverts);
			vt(Eigen::seq(nverts, 4 * nverts - 1)) = object->getVelocities();
			vt(Eigen::seq(0, nverts - 1)) = object->getAmpVelocities();

			object->clampDofs(vt);
			// external forces
			argus::VectorN aext = getAdaptiveAccelerationForVBD(*object, xt, true, true);
			argus::VectorN y = timestep * vt + timestep * timestep * aext;
			// TODO: Decide whether to use or not?
			// wme->clampWrinkleDofs(y);
			y = xt + y;
			// get collision free updates
			// wme->getObstacleProjectedPositions(obs, y);
			// y =  getCollisionFreeState(obs, y, xt, timestep, NULL, frame);
			// y= getCollisionFreeState(obs, y, xt, object->getTriangles(), timestep, NULL, frame);
			VectorN posxt, posy, velxt;
			posy = y(Eigen::seq(nverts, 4 * nverts - 1));
			posxt = xt(Eigen::seq(nverts, 4 * nverts - 1));
			velxt = vt(Eigen::seq(nverts, 4 * nverts - 1));
			getCollisionFreeState(obs, posy, velxt, posxt, invInertia, object->getTriangles(), timestep, NULL, frame);

			y(Eigen::seq(nverts, 4 * nverts - 1)) = posy;

			argus::VectorN x = y;
			Eigen::MatrixXi F = object->getFaces();
			int niters = m_pos_vbd_iterations;
#ifdef ARGUS_CHECKPOINT
			if (argus::simConf.contains("VBD Iterations for Pos") == false)
				argus::simConf["VBD Iterations for Pos"] = niters;
#endif
			std::vector<int> timings;
			std::vector<unsigned int> collverts;
			for (int iter = 0; iter < niters; iter++)
			{
				int ncolors = vbdInfo.coloredSets.rows();
				// Inside the loop
				for (int color = 0; color < ncolors; color++)
				{
					argus::VectorNi set = vbdInfo.coloredSets.row(color);

					// Eigen::internal::set_is_malloc_allowed(false);
					// auto t1 = argus::util::captureCurrentTime();

#ifdef USE_OMP
#pragma omp parallel for schedule(static)
#endif
					for (int vtex = 0; vtex < set.rows(); vtex++)
					{
						int vid = set(vtex);
						if (vid == -1)
							continue;
						if (invInertia(3 * vid) == 0)
							continue;
						Vector4 xi, yi;
						xi(0) = x(vid);
						xi.segment<3>(1) = x(Eigen::seq(nverts + 3 * vid, nverts + 3 * vid + 2));
						yi(0) = y(vid);
						yi.segment<3>(1) = y(Eigen::seq(nverts + 3 * vid, nverts + 3 * vid + 2));
						Vector4 deriv;
						Matrix44 hess;
						argus::Real energy = object->getVertexEnergy4(vid, &deriv, &hess);
						// since amp acceleration is twice as fast as pos acc
						// deriv(0) *= 2; // double force is applied
						// hess.row(0) *= 2;
						// hess.col(0) *= 2;
						// hess(0, 0) /= 2;
						Real mvtex = mass(vid);
						Vector4 fi = -(mvtex / (timestep * timestep)) * (xi - yi) - deriv;
						Matrix44 Hi = (mvtex / (timestep * timestep)) * argus::Matrix44::Identity() + hess;

						// do not change if Hi is singular
						if (std::fabs(Hi.determinant()) < 1e-5)
							continue;
						argus::Vector4 dxi = Hi.inverse() * fi;
						// project for clamping
						// wme->clampVertexWrinkleDofs(dxi, vid);
						argus::Real alpha = 1;
						// alpha = wme->backtrackingLineSearch(func, xi, dxi, alpha, 10);
						auto func = [&object, &xi, &timestep, &yi, &mvtex, &dxi, vid, &x, &nverts](Real rate)
						{
							object->getAmpVec()(vid) = xi(0) + rate * dxi(0);
							object->getPositions()(Eigen::seq(3 * vid, 3 * vid + 2)).noalias() = xi.segment<3>(1) + rate * dxi.segment<3>(1);
							Real e = object->getVertexEnergy4(vid, NULL, NULL);
							object->getPositions()(Eigen::seq(3 * vid, 3 * vid + 2)).noalias() = x(Eigen::seq(nverts + 3 * vid, nverts + 3 * vid + 2));
							object->getAmpVec()(vid) = x(vid);
							e += (mvtex / (2 * timestep * timestep)) * (xi + rate * dxi - yi).transpose() * (xi + rate * dxi - yi);

							return e;
						};

						Real f0 = energy + (mvtex / (2 * timestep * timestep)) * (xi - yi).transpose() * (xi - yi);
						// TODO: test if func(0) == energy
						Real f1 = func(1);
						alpha = 1;
						if (f1 < f0)
						{
							Vector4 update = xi + alpha * dxi;
							x(vid) = update(0);
							x(Eigen::seq(nverts + 3 * vid, nverts + 3 * vid + 2)).noalias() = update.segment<3>(1);
						}
					}

					// auto t2 = argus::util::captureCurrentTime();
					// auto duration = m_stats.value("init", 0LL) + std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
					// timings.push_back(duration);
					// Eigen::internal::set_is_malloc_allowed(true);
					// update positions
					object->getPositions().noalias() = x(Eigen::seq(nverts, 4 * nverts - 1));
					object->getAmpVec().noalias() = x(Eigen::seq(0, nverts - 1));
					object->applyDynamicConstraints(timestep);
					object->applyTopologicalConstraints();
				}
				x(Eigen::seq(nverts, 4 * nverts - 1)) = object->getPositions();
				x(Eigen::seq(0, nverts - 1)) = object->getAmpVec();
				argus::VectorN vnew
				    = (x - xt) / timestep;
				if (selfCollisionParams.enable_self_collisions)
				{
					VectorN posx, velx, posxold;
					posx = x(Eigen::seq(nverts, 4 * nverts - 1));
					posxold = xt(Eigen::seq(nverts, 4 * nverts - 1));
					velx = vnew(Eigen::seq(nverts, 4 * nverts - 1));
					resolveSelfCollisionsSubstepped(
					    object->getTriangles(),
					    posxold, // Position at the very start of the frame
					    posx, // Predicted position (modified to resolved position)
					    velx, // Velocity (modified)
					    selfCollisionParams.self_collision_maxsubsteps,
					    invInertia,
					    selfCollisionData,
					    selfCollisionParams.stiffness, // stiffness
					    selfCollisionParams.damping, // damping
					    timestep);
					x(Eigen::seq(nverts, 4 * nverts - 1)) = posx;
					xt(Eigen::seq(nverts, 4 * nverts - 1)) = posxold;
					vnew(Eigen::seq(nverts, 4 * nverts - 1)) = velx;
				}
				// update velocities

				object->clampDofs(vnew);
				// ether damping
				vnew(Eigen::seq(0, nverts - 1)) *= (1 - m_kdamp_amps);
				vnew(Eigen::seq(nverts, 4 * nverts - 1)) *= (1 - m_kdamp_pos);
				x = xt + vnew * timestep;
				VectorN posxt, posx, vnewx;
				posx = x(Eigen::seq(nverts, 4 * nverts - 1));
				posxt = xt(Eigen::seq(nverts, 4 * nverts - 1));
				vnewx = vnew(Eigen::seq(nverts, 4 * nverts - 1));

				resolveCollisions(obs, posx, posxt, vnewx, invInertia, object->getTriangles(), timestep, state, frame, collverts);
				// x = xt + vnew * timestep;
				x(Eigen::seq(nverts, 4 * nverts - 1)) = posx;
				vnew(Eigen::seq(nverts, 4 * nverts - 1)) = vnewx;
				object->getPositions() = x(Eigen::seq(nverts, 4 * nverts - 1));
				object->getVelocities() = vnew(Eigen::seq(nverts, 4 * nverts - 1));
				object->getAmpVec() = x(Eigen::seq(0, nverts - 1));
				object->getAmpVelocities() = vnew(Eigen::seq(0, nverts - 1));
				// this function updates positions as well as velocities
			}
#ifdef ARGUS_CHECKPOINT
#endif
		}
	}

}
