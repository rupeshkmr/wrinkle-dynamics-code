#include <argus/vbd.hpp>
#include <argus/self_collisions.hpp>
#include <argus/util.hpp>

namespace argus
{

	VBD::~VBD() = default;

	void VBD::init(std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs, nlohmann::json& params)
	{
		/**
		 * Initialize the timestepper
		 * Current config allows global dynamic damping parameters for all sheets.
		 * In future we may want multiple sheets with different parameters.
		 */

		// default values
		m_stats["vbdUpdateAmps"] = 0;
		m_stats["vbdUpdatePos"] = 0;
		m_stats["optimizeDphis"] = 0;
		m_stats["constraint"] = 0;
		m_stats["init"] = 0;

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
		if (params.contains("self_collision_threshold"))
		{
			selfCollisionParams.threshold = params["self_collision_threshold"].get<Real>();
		}
		else
			selfCollisionParams.threshold = 0.01;
		if (params.contains("self_collision_stiffness"))
			selfCollisionParams.stiffness = params["self_collision_stiffness"].get<Real>();
		else
			selfCollisionParams.stiffness = 0.01;
		if (params.contains("self_collision_damping"))
			selfCollisionParams.damping = params["self_collision_damping"].get<Real>();
		else
			selfCollisionParams.damping = 0.01;
		if (params.contains("sdf_collider"))
			m_sdfCollider = params["sdf_collider"].get<std::string>();
		else
			m_sdfCollider = "basic";

		m_kdamp_pos = 0;
		m_sim_friction = false;
		m_pos_vbd_iterations = 10;
		if (params.contains("kdamp_pos"))
			m_kdamp_pos = params["kdamp_pos"].get<Real>();
		if (params.contains("sim_friction"))
		{
			m_sim_friction = params["sim_friction"].get<unsigned int>();
		}
		if (params.contains("pos_vbd_iter"))
			m_pos_vbd_iterations = params["pos_vbd_iter"].get<unsigned int>();

#ifdef ARGUS_CHECKPOINT
		argus::simConf["Friction Enabled"] = m_sim_friction;
#endif
		initMeshes(dynamicObjects, obs);
		initSelfCollision(dynamicObjects[0]->getTriangles(), dynamicObjects[0]->getRestPositions());
		initVBD(dynamicObjects);
	}
	void VBD::initVBD(std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects)
	{
		for (auto&& object : dynamicObjects)
		{
			object->requireVertexColoringInfo();
		}
	}
	void VBD::initMeshes(std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obstacles)
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
	void VBD::initSelfCollision(const std::vector<trimesh::triangle_t>& faces, const VectorN& Xrest)
	{
		// 1. Fill the 15 weights into staticData.flattenedWeights using your helper
		initializeStaticWeights(selfCollisionData);

		// 2. Compute radii per face/probe and determine the grid size (h) using your helper
		// This fills selfCollisionData.flattenedRadii and selfCollisionData.gridCellSize
		precomputeRestRadiiAndGrid(
		    faces,
		    Xrest,
		    selfCollisionData,
		    0.55, // safetyMultiplier
		    0.0001 // minThickness
		);

		// 3. One-time Setup: Initialize the Point Cloud Buffer Metadata
		// This "bakes" the static info into the 128-byte aligned structs
		// so the solver doesn't have to look them up later.
		const int nFaces = faces.size();
		const int K = selfCollisionData.K;
		selfCollisionData.pointCloudBuffer.resize(nFaces * K);

#pragma omp parallel for schedule(static)
		for (int f = 0; f < nFaces; ++f)
		{
			const auto& face = faces[f];
			for (int i = 0; i < K; ++i)
			{
				auto& pt = selfCollisionData.pointCloudBuffer[f * K + i];

				// --- Store Static Data (Never changes during simulation) ---
				pt.radius = selfCollisionData.flattenedRadii[f * K + i];
				pt.faceId = f;
				pt.vIdx[0] = face.v[0];
				pt.vIdx[1] = face.v[1];
				pt.vIdx[2] = face.v[2];

				// Store weights locally in the point for fast repulsion math
				pt.b[0] = selfCollisionData.flattenedWeights[i * 3 + 0];
				pt.b[1] = selfCollisionData.flattenedWeights[i * 3 + 1];
				pt.b[2] = selfCollisionData.flattenedWeights[i * 3 + 2];

				// Initialize position to zero (will be updated in the dynamic Gather phase)
				pt.pos.setZero();

				// Clear padding for memory alignment safety
				for (int p = 0; p < 7; ++p)
					pt.padding[p] = 0.0;
			}
		}
	}
	void VBD::resolveSelfCollisionsSubstepped(
	    const std::vector<trimesh::triangle_t>& faces,
	    const VectorN& xt, // Position at start of frame
	    VectorN& y, // Target position
	    VectorN& v, // Frame velocity
	    int maxSubsteps,
	    const VectorN& invInertia,
	    PrecomputedSelfCollisionData& sd,
	    Real stiffness,
	    Real damping,
	    Real totalDt)
	{
		Real substepDt = totalDt / static_cast<Real>(maxSubsteps);

		// 1. SYNC INITIAL VELOCITY
		// Ensure the velocity we start with matches the total displacement predicted
		v = (y - xt) / totalDt;

		VectorN currentPos = xt;

		for (int i = 0; i < maxSubsteps; ++i)
		{
			// 2. PREDICT SUBSTEP POSITION
			VectorN nextPos = currentPos + v * substepDt;

			// 3. RESOLVE (Solver now uses and updates 'v' correctly for 'substepDt')
			selfcollisionhandler::solveAllSelfCollisionsRepulsions(
			    faces,
			    currentPos, // Xold
			    nextPos, // Xnew
			    v, // V (Correctly scaled to substepDt inside)
			    invInertia,
			    sd,
			    stiffness,
			    damping,
			    substepDt // Passing the tiny DT
			);

			// 4. ADVANCE
			currentPos = nextPos;
		}

		// 5. FINAL SYNC
		y = currentPos;
	}
	void VBD::addExternalForces(std::vector<std::unique_ptr<DynamicsObject>>& objects, Real timestep)
	{
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
		// rotation for cylinder twist
		// std::vector<int> verts = { 4160, 4161, 4162, 4163, 4164, 4165, 4166, 4167, 4168, 4169, 4170, 4171, 4172, 4173, 4174, 4175, 4176, 4177, 4178, 4179, 4180, 4181, 4182, 4183, 4184, 4185, 4186, 4187, 4188, 4189, 4190, 4191, 4192, 4193, 4194, 4195, 4196, 4197, 4198, 4199, 4200, 4201, 4202, 4203, 4204, 4205, 4206, 4207, 4208, 4209, 4210, 4211, 4212, 4213, 4214, 4215, 4216, 4217, 4218, 4219, 4220, 4221, 4222, 4223, 4224 };

		// Real theta = 0.05 * timestep;
		// for (auto&& it : verts)
		// {
		// 	// rotate about z
		// 	Matrix33 rotation;
		// 	rotation << cos(theta), -sin(theta), 0,
		// 	    sin(theta), cos(theta), 0,
		// 	    0, 0, 1;
		// 	int i = it;
		// 	std::cout << i << std::endl;
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
	}

	void VBD::applyConstraints(Real timestep, std::vector<std::unique_ptr<DynamicsObject>>& objects)
	{
		for (auto&& object : objects)
		{
			object->applyDynamicConstraints(timestep);
			object->applyTopologicalConstraints();
		}
	}

	void VBD::step(Real timestep, std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs, nlohmann::json& state, unsigned int frame)
	{
		addExternalForces(dynamicObjects, timestep);
		// dynamic step
		updatePos(timestep, dynamicObjects, obs, state, frame);

		// Apply constraints on all ojects
		applyConstraints(timestep, dynamicObjects);
	}

	void VBD::updatePos(Real timestep, std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs, nlohmann::json& state, unsigned int frame)
	{
		auto t3 = argus::util::captureCurrentTime();

		vbdUpdatePos(timestep, dynamicObjects, obs, state, frame);
		auto t4 = argus::util::captureCurrentTime();
		m_stats["vbdUpdatePos"] = m_stats.value("vbdUpdatePos", 0LL) + std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count();
	}

	VectorN VBD::getAdaptiveAccelerationForVBD(DynamicsObject& object, VectorN& var, bool pos)
	{
		unsigned int vertexCount = object.getVertexCount();
		unsigned int faceCount = object.getFaceCount();
		VectorN invInertia = object.getInvInertia();
		SheetInfo vbdInfo = object.getSheetInfo();
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

	void VBD::vbdUpdatePos(Real timestep, std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs, nlohmann::json& state, unsigned int frame)
	{

		// TODO: Remember to save old positions before calling this function to resolve collisions further.
		for (auto&& object : dynamicObjects)
		{
			SheetInfo vbdInfo = object->getSheetInfo();
			// get mass information
			argus::VectorN mass = object->getInertia();
			VectorN invInertia = object->getInvInertia();
			// update positions from external forces and inertia
			// adaptive intertia and acceleration
			argus::VectorN xt = object->getPositions();
			argus::VectorN vt = object->getVelocities();
			// external forces
			argus::VectorN aext = getAdaptiveAccelerationForVBD(*object, xt, true);
			argus::VectorN y = timestep * vt + timestep * timestep * aext;
			// TODO: Decide whether to use or not?
			// wme->clampWrinkleDofs(y);
			y = xt + y;
			// get collision free updates
			// wme->getObstacleProjectedPositions(obs, y);
			// getCollisionFreeState(obs, y, xt, timestep, NULL, frame);
			getCollisionFreeState(obs, y, vt, xt, invInertia, object->getTriangles(), timestep, NULL, frame);

			argus::VectorN x = y;
			int niters = m_pos_vbd_iterations;
#ifdef ARGUS_CHECKPOINT
			if (argus::simConf.contains("VBD Iterations for Pos") == false)
				argus::simConf["VBD Iterations for Pos"] = niters;
#endif
			std::vector<unsigned int> collverts;
			for (int iter = 0; iter < niters; iter++)
			{
				int ncolors = vbdInfo.coloredSets.rows();
				// Inside the loop
				for (int color = 0; color < ncolors; color++)
				{
					argus::VectorNi set = vbdInfo.coloredSets.row(color);
#ifdef USE_OMP
#pragma omp parallel for
#endif
					for (int vtex = 0; vtex < set.rows(); vtex++)
					{
						int vid = set(vtex);
						if (vid == -1)
							continue;
						if (invInertia(3 * vid) == 0)
							continue;
						Vector3 xi;
						xi = x(Eigen::seq(3 * vid, 3 * vid + 2));
						Vector3 yi;
						yi = y(Eigen::seq(3 * vid, 3 * vid + 2));
						Vector3 deriv;
						Matrix33 hess;
						Real energy = object->getVertexEnergy(vid, &deriv, &hess);
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

					// update positions
					object->getPositions() = x;
					object->applyDynamicConstraints(timestep);
				}
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

				// resolveCollisions(obs, x, xt, vnew, timestep, state, frame, collverts);
				object->getPositions() = x;
				object->getVelocities() = vnew;
				// this function updates positions as well as velocities
				// object->resolveWrinkleSelfCollisions(xt, timestep);
			}
		}
	}
	void VBD::getCollisionFreeState(std::vector<std::unique_ptr<Obstacle>>& obstacles, VectorN& X, VectorN& V, const VectorN& Xold, const VectorN& W, const std::vector<trimesh::triangle_t>& triangles, Real timestep, std::vector<int>* collVerts, unsigned int frame)
	{
		/**
		 * Project updated_pos to collision free states.
		 */
		int nobs = obstacles.size();
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
	void VBD::getCollisionFreeState(std::vector<std::unique_ptr<Obstacle>>& obstacles, VectorN& X, const VectorN& Xold, Real timestep, std::vector<int>* collVerts, unsigned int frame)
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
	void VBD::getCollisionFreeState(std::vector<std::unique_ptr<Obstacle>>& obstacles, VectorN& X, const VectorN& Xold, const std::vector<trimesh::triangle_t>& triangles, Real timestep, std::vector<int>* collVerts, unsigned int frame)

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
	void VBD::resolveCollisions(std::vector<std::unique_ptr<Obstacle>>& obs, VectorN& updated_pos, const VectorN& oldPositions, VectorN& velocities, const VectorN& W, const std::vector<trimesh::triangle_t>& triangles, Real timestep, nlohmann::json& state, unsigned int frame, std::vector<unsigned int>& collVerts)
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

	// 	void VBD::getCollisionFreeState(std::vector<std::unique_ptr<Obstacle>>& obstacles, VectorN& X, const VectorN& Xold, Real timestep, std::vector<int>* collVerts, unsigned int frame)
	// 	{
	// 		/**
	// 		 * Project updated_pos to collision free states.
	// 		 */
	// 		throw "not implemented!";
	// 		// int nobs = obstacles.size();
	// 		// for (int i = 0; i < nobs; i++)
	// 		// {
	// 		// 	obstacles.at(i)->project(X, Xold, timestep, frame, collVerts);
	// 		// }
	// 	}
	// 	void VBD::getCollisionFreeState(std::vector<std::unique_ptr<Obstacle>>& obstacles, VectorN& X, const VectorN& Xold, const std::vector<trimesh::triangle_t>& triangles, Real timestep, std::vector<int>* collVerts, unsigned int frame)

	// 	{
	// 		/**
	// 		 * Project updated_pos to collision free states.
	// 		 */
	// 		int nobs = obstacles.size();
	// 		for (int i = 0; i < nobs; i++)
	// 		{
	// 			X = obstacles.at(i)->project(X, Xold, triangles, timestep, frame, collVerts);
	// 		}
	// 	}
	// 	void VBD::getCollisionFreeState(std::vector<std::unique_ptr<Obstacle>>& obstacles, VectorN& X, VectorN& V, const VectorN& Xold, const VectorN& W, const std::vector<trimesh::triangle_t>& triangles, Real timestep, std::vector<int>* collVerts, unsigned int frame)
	// 	{
	// 		/**
	// 		 * Project updated_pos to collision free states.
	// 		 */
	// 		int nobs = obstacles.size();
	// 		for (int i = 0; i < nobs; i++)
	// 		{
	// 			if (m_sdfCollider == "basic")
	// 				obstacles.at(i)->project(X, V, Xold, timestep, frame, collVerts);
	// 			else if (m_sdfCollider == "quadrature")
	// 				obstacles.at(i)->project(timestep, X, V, Xold, W, triangles, frame, collVerts, &m_config);
	// 			else
	// 				throw "not implemented!";
	// 		}
	// 	}

	// 	void VBD::resolveCollisions(std::vector<std::unique_ptr<Obstacle>>& obs, VectorN& updated_pos, const VectorN& oldPositions, VectorN& velocities, Real timestep, nlohmann::json& state, unsigned int frame, std::vector<unsigned int>& collVerts)
	// 	{
	// 		/**
	// 		 * Project sheet positions to collision free states after resolving collisions from all obstacles.
	// 		 */
	// 		if (obs.size() != 0)
	// 		{
	// 			unsigned int coliters = 10;
	// #ifdef ARGUS_CHECKPOINT
	// 			if (simConf.contains("Obstacle Collision iterations | All objects") == false)
	// 				simConf["Obstacle Collision iterations | All objects"] = coliters;
	// #endif
	// 			// create a copy if normal info is required.
	// 			VectorN inertialPos;
	// 			std::vector<int> collisionVertices;

	// 			// for friction
	// 			if (m_sim_friction)
	// 				inertialPos = updated_pos;
	// 			for (int i = 0; i < coliters; i++)
	// 			{
	// #ifdef ARGUS_CHECKPOINT
	// 				// getCollisionFreeState(obs, updated_pos, &collisionVertices, frame);

	// 				getCollisionFreeState(obs, updated_pos, oldPositions, timestep, &collisionVertices, frame);
	// #else
	// 				getCollisionFreeState(obs, updated_pos, oldPositions, timestep, NULL, frame);
	// #endif
	// 			}
	// 			// update velocities
	// 			velocities = (updated_pos - oldPositions) / timestep;
	// 			// add friciton damping
	// 			// TODO: Right now I am adding global friction coefficient, an elegant solution would be adding friction for individual objects
	// #ifdef ARGUS_CHECKPOINT

	// 			if (collisionVertices.size() == 0)
	// 			{
	// 				state["collision_vertices"] = collisionVertices;
	// 			}
	// 			else
	// 			{
	// 				// store state
	// 				std::set<int> cset(collisionVertices.begin(), collisionVertices.end());
	// 				std::vector<unsigned int> cvec(cset.begin(), cset.end());
	// 				Eigen::VectorXi collverts;
	// 				collverts.resize(cvec.size());
	// 				for (unsigned int i = 0; i < cvec.size(); i++)
	// 					collverts(i) = cvec.at(i);
	// 				state["collision_vertices"] = collverts;
	// 			}
	// #endif
	// 			// old collision resolution pipeline
	// 			/* if (m_sim_friction && std::fabs(frictionK) > 1e-6)
	// 			{

	// 			    std::set<int> collVertsSet(collisionVertices.begin(), collisionVertices.end());
	// 			    std::vector<unsigned int> collVertsVec(collVertsSet.begin(), collVertsSet.end());
	// 			    unsigned int collVertsSize = collVertsVec.size();
	// #ifdef USE_OMP
	// #pragma omp parallel for
	// #endif
	// 			    for (unsigned int i = 0; i < collVertsSize; i++)
	// 			    {
	// 			        unsigned int vid = collVertsVec.at(i);
	// 			        Vector3 pos1, pos2, ncap;
	// 			        pos1 = inertialPos(Eigen::seq(3 * vid, 3 * vid + 2));
	// 			        pos2 = updated_pos(Eigen::seq(3 * vid, 3 * vid + 2));
	// 			        ncap = (pos2 - pos1).normalized();
	// 			        Vector3 vTan = velocities(Eigen::seq(3 * vid, 3 * vid + 2)) - velocities(Eigen::seq(3 * vid, 3 * vid + 2)).dot(ncap) * ncap;
	// 			        velocities(Eigen::seq(3 * vid, 3 * vid + 2)) -= vTan;
	// 			        velocities(Eigen::seq(3 * vid, 3 * vid + 2)) += (1 - frictionK) * vTan;
	// 			        velocities(Eigen::seq(3 * vid, 3 * vid + 2))(0) = 0;
	// 			        updated_pos(Eigen::seq(3 * vid, 3 * vid + 2)) = oldPositions(Eigen::seq(3 * vid, 3 * vid + 2)) + velocities(Eigen::seq(3 * vid, 3 * vid + 2)) * timestep;
	// 			    }
	// 			} */
	// 		}
	// }
}
