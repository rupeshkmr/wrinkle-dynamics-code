#include <argus/wrinkles_solver.hpp>
#include <argus/lbfgsb_wrapper.hpp>
#include <argus/lbfgs_wrapper.hpp>
#include <dlib/optimization.h>
namespace argus
{

	void WrinkleSolver::init(std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, nlohmann::json& params)
	{
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

		initWrinkleMesh(dynamicObjects, true, true);
		// for vbd
		for (auto&& object : dynamicObjects)
		{
			object->requireVertexColoringInfo();
		}
	}

	void WrinkleSolver::initWrinkleMesh(std::vector<std::unique_ptr<DynamicsObject>>& objects, bool amps, bool dphis)
	{
		if (amps == false && dphis == false)
			throw "Error! Calling initWrinkleMesh without chosing a parameter.\n";
		for (auto&& object : objects)
		{
			// TODO: Remove copying here
			VectorN a;
			Matrix2N w;
			object->getInitialWrinkleParameters(0.001, a, w, amps);
			if (amps)
				object->getAmpVec() = a;
			if (dphis)
				object->getDphis() = w;
			std::cout << "initialized wrinkle parameters " << std::endl;
		}
	}

	void WrinkleSolver::optimizeDphisParallel(std::vector<std::unique_ptr<DynamicsObject>>& objects)
	{
		if (m_optimizer)
			for (auto&& object : objects)
			{
				unsigned int nfaces = object->getFaceCount();

#ifdef USE_OMP
#pragma omp parallel for
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

	void WrinkleSolver::optimizeDphisGlobal(std::vector<std::unique_ptr<DynamicsObject>>& objects)
	{
		if (m_optimizer)
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
					ip(Eigen::seq(0, vertexCount - 1)) = object->getAmpVec();
					ip(Eigen::seq(vertexCount, vertexCount + 2 * faceCount - 1)) = X;
					ip(Eigen::seq(vertexCount + 2 * faceCount, 4 * vertexCount + 2 * faceCount - 1)) = object->getPositions().reshaped();
					Real e = object->getEnergies(ip, &fderiv);
					deriv = fderiv(Eigen::seq(vertexCount, vertexCount + 2 * faceCount - 1));
					return e;
				};
				// derivative test
				/* VectorN ip = VectorN::Random(faceCount * 2) * 9998;
				VectorN dx = VectorN::Random(faceCount * 2).normalized();
				Real h = 0.9;
				Real e0, e1, e2;
				VectorN d0, d1, d2, X1, X2;
				e0 = func(ip, d0);
				VectorN errors(10);
				errors.setZero();
				for (int i = 0; i < 10; i++)
				{
				    X1 = ip + 0.5 * h * dx;
				    e1 = func(X1, d1);
				    X2 = ip - 0.5 * h * dx;
				    e2 = func(X2, d2);
				    Real error = std::fabs(((e1 - e2) / h) - dx.dot(d0));
				    std::cout << "Error " << error << std::endl;
				    errors(i) = error;
				    h *= 0.1;
				}
				exit(0); */
				VectorN X = object->getDphis().reshaped();
				try
				{
					m_optimizer->solve(func, X);
					// update dphi
					object->getDphis() = X.reshaped(2, faceCount);
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

	void WrinkleSolver::optimizeAmplitudesGlobal(std::vector<std::unique_ptr<DynamicsObject>>& objects)
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

			/* // deriv test
			VectorN ip, dx, d0, d1, d2, X1, X2;
			Real e0, e1, e2;
			ip = VectorN::Random(vertexCount);
			dx = VectorN::Random(vertexCount).normalized();
			Real h = 0.1;
			e0 = func(ip, d0);
			for (int i = 0; i < 10; i++)
			{
			    X1 = ip + 0.5 * h * dx;
			    e1 = func(X1, d1);
			    X2 = ip - 0.5 * h * dx;
			    e2 = func(X2, d2);
			    Real error = std::fabs(((e1 - e2) / h) - dx.dot(d0));
			    std::cout << "Error " << error << std::endl;
			    h *= 0.1;
			}
			exit(0); */
			VectorN lb, ub;
			ub.resize(vertexCount);
			lb.resize(vertexCount);
			// set the bounds
			for (int i = 0; i < vertexCount; i++)
			{
				ub(i) = std::numeric_limits<argus::Real>::infinity();
				lb(i) = 1e-3;
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

	void WrinkleSolver::optimizeAmplitudesNew(std::vector<std::unique_ptr<DynamicsObject>>& objects)
	{
		for (auto&& object : objects)
		{
			VectorN invInertia = object->getInvInertia();
			VectorN x = object->getAmpVec();
			SheetInfo vbdInfo = object->getSheetInfo();
			int niters = 20;
#ifdef ARGUS_CHECKPOINT
			if (argus::simConf.contains("VBD Iterations for Amps") == false)
				argus::simConf["VBD Iterations for Amps"] = niters;
#endif
			for (int iter = 0; iter < niters; iter++)
			{
				int ncolors = vbdInfo.coloredSets.rows();
				// Inside the loop
				for (int color = 0; color < ncolors; color++)
				{
					argus::VectorNi vset = vbdInfo.coloredSets.row(color);
#ifdef USE_OMP
#pragma omp parallel for
#endif
					for (int vtex = 0; vtex < vset.size(); vtex++)
					{
						int vid = vset(vtex);
						if (vid == -1)
							continue;
						if (invInertia(3 * vid) == 0)
							continue;
						auto ehat = [&object, &vid, this](Real aval, Real* grad, Real* hess) -> Real
						{
							Real abak = object->getAmpVec()(vid);
							VectorN dphisbak = getAdjacentDphis(vid, object);
							VectorN dhda;
							VectorN dphiNew = getNewDphi((grad || hess) ? &dhda : NULL, aval, vid, object);
							setAdjacentDphis(vid, dphiNew, object);
							object->getAmpVec()(vid) = aval;
							std::vector<int> fids;
							VectorN deriv;
							MatrixNN jac;
							Real energy = object->getVertexEnergy1(vid, fids, (grad || hess) ? &deriv : NULL, hess ? &jac : NULL);
							if (grad)
							{
								Real dEda = deriv(0);
								VectorN dEdha = deriv(Eigen::seq(1, Eigen::last));
								if (dEdha.size() != dhda.size())
									throw " Incompatible vector sizes!";
								*grad = dEda + dEdha.dot(dhda);
							}
							if (hess)
							{
								Real d2Eda2 = jac(0, 0);
								int ndofs = dphiNew.size();
								VectorN d2Edadha(ndofs);
								d2Edadha.setZero();
								d2Edadha = jac.col(0)(Eigen::seq(1, ndofs));
								MatrixNN d2Edh2 = jac.block(1, 1, ndofs, ndofs);
								if (d2Edadha.size() != dhda.size())
									throw " Incompatible vector sizes!";
								if (d2Edh2.cols() != dhda.size() || d2Edh2.rows() != dhda.size())
									throw " Incompatible matrix vector product!";
								*hess = d2Eda2 + 2 * d2Edadha.dot(dhda) + dhda.transpose() * d2Edh2 * dhda;
							}
							setAdjacentDphis(vid, dphisbak, object);
							object->getAmpVec()(vid) = abak;
							return energy;
						};
						/* // derivative test
						                        std::cout << "Vid " << vid << std::endl;
						                        Real ip, x1, x2, d0, d1, d2, h0, h1, h2;
						                        Real h = 0.89;
						                        Real ee0, ee1, ee2;
						                        ip = 1.2;
						                        object->getAmpVec()(vid) = ip;
						                        ee0 = object->getVertexEnergy1(vid, &d0, &h0);
						                        for (int i = 0; i < 10; i++)
						                        {

						                            object->getAmpVec()(vid) = ip + 0.5 * h;
						                            ee1 = object->getVertexEnergy1(vid, &d1, &h1);
						                            object->getAmpVec()(vid) = ip - 0.5 * h;
						                            ee2 = object->getVertexEnergy1(vid, &d2, &h2);
						                            Real derror = std::fabs(((ee1 - ee2) / h) - d0);
						                            Real herror = std::fabs(((d1 - d2) / h) - h0);
						                            // std::cout << derror << std::endl;
						                            h *= 0.1;
						                            std::cout << herror << std::endl;
						                        }

						                        std::cout << "end" << vid << std::endl; */
						argus::Real xi;
						xi = x(vid);
						argus::Real deriv;
						argus::Real hess;
						argus::Real energy = ehat(xi, &deriv, &hess);

						argus::Real fi
						    = -deriv; // here a factor of 2 due to energy
						argus::Real Hi = hess; // here a factore of 2 due to energy
						// do not change if Hi is singular
						if (std::fabs(Hi) < 1e-6)
							continue;
						if (Hi < 0)
							Hi *= -1;
						argus::Real dxi
						    = fi / Hi; // a multiplier of 2 as amplitude equation involves double increase in acceleration

						argus::Real alpha = 1;
						energy = object->getVertexEnergy1(vid);
						object->getAmpVec()(vid) = xi + alpha * dxi;
						Real e1;
						while (true)
						{
							object->getAmpVec()(vid) = xi + alpha * dxi;
							e1 = object->getVertexEnergy1(vid, NULL, NULL);
							if (e1 < energy)
								break;
							if (alpha < 1e-5)
							{
								alpha = 0;
								break;
							}
							alpha *= 0.5;
						}
						x(vid) = xi + alpha * dxi;
						if (alpha > 0)
						{
							// update dphis
							VectorN dphiNew = getNewDphi(NULL, x(vid), vid, object);
							setAdjacentDphis(vid, dphiNew, object);
						}
						if (x(vid) < 1e-4)
							x(vid) = 1e-4;
					}
					// update positions
					object->getAmpVec() = x;
				}
				// clamp
				object->clampDofs(x); // update
				object->getAmpVec() = x;
			}
		}
	}

	void WrinkleSolver::optimizeAmplitudesVBD(std::vector<std::unique_ptr<DynamicsObject>>& objects)
	{
		for (auto&& object : objects)
		{
			VectorN invInertia = object->getInvInertia();
			VectorN x = object->getAmpVec();
			SheetInfo vbdInfo = object->getSheetInfo();
			int niters = 100;
#ifdef ARGUS_CHECKPOINT
			if (argus::simConf.contains("VBD Iterations for Amps") == false)
				argus::simConf["VBD Iterations for Amps"] = niters;
#endif
			for (int iter = 0; iter < niters; iter++)
			{
				int ncolors = vbdInfo.coloredSets.rows();
				// Inside the loop
				for (int color = 0; color < ncolors; color++)
				{
					argus::VectorNi vset = vbdInfo.coloredSets.row(color);
#ifdef USE_OMP
#pragma omp parallel for
#endif
					for (int vtex = 0; vtex < vset.size(); vtex++)
					{
						int vid = vset(vtex);
						if (vid == -1)
							continue;
						if (invInertia(3 * vid) == 0)
							continue;
						argus::Real xi;
						xi = x(vid);
						argus::Real deriv;
						argus::Real hess;
						argus::Real energy = object->getVertexEnergy1(vid, &deriv, &hess);
						/* // derivative test
						std::cout << "Vid " << vid << std::endl;
						Real ip, x1, x2, d0, d1, d2, h0, h1, h2;
						Real h = 0.89;
						Real ee0, ee1, ee2;
						ip = 1.2;
						object->getAmpVec()(vid) = ip;
						ee0 = object->getVertexEnergy1(vid, &d0, &h0);
						for (int i = 0; i < 10; i++)
						{

						    object->getAmpVec()(vid) = ip + 0.5 * h;
						    ee1 = object->getVertexEnergy1(vid, &d1, &h1);
						    object->getAmpVec()(vid) = ip - 0.5 * h;
						    ee2 = object->getVertexEnergy1(vid, &d2, &h2);
						    Real derror = std::fabs(((ee1 - ee2) / h) - d0);
						    Real herror = std::fabs(((d1 - d2) / h) - h0);
						    // std::cout << derror << std::endl;
						    h *= 0.1;
						    std::cout << herror << std::endl;
						}

						std::cout << "end" << vid << std::endl; */
						argus::Real fi
						    = -deriv; // here a factor of 2 due to energy
						argus::Real Hi = hess; // here a factore of 2 due to energy
						// do not change if Hi is singular
						if (std::fabs(Hi) < 1e-6)
							continue;
						if (Hi < 0)
							Hi *= -1;
						argus::Real dxi
						    = fi / Hi; // a multiplier of 2 as amplitude equation involves double increase in acceleration

						/* object->getAmpVec()(vid) = xi + alpha * dxi;
						Real e1;
						while (true)
						{
						    object->getAmpVec()(vid) = xi + alpha * dxi;
						    e1 = object->getVertexEnergy1(vid, NULL, NULL);
						    if (e1 < energy)
						        break;
						    if (alpha < 1e-5)
						    {
						        alpha = 0;
						        break;
						    }
						    alpha *= 0.5;
						} */
						auto func = [&object, &xi, &dxi, vid](Real a)
						{
							object->getAmpVec()(vid) = xi + a * dxi;
							Real bak = object->getAmpVec()(vid);
							Real e = object->getVertexEnergy1(vid);
							object->getAmpVec()(vid) = bak;
							return e;
						};
						Real f0 = func(0);
						argus::Real alpha;
						alpha = dlib::backtracking_line_search(func, f0, -fi * dxi, 1, 0.5, 3);
						if (alpha > 1e-4)
						{
							object->getAmpVec()(vid) = xi + alpha * dxi;
							Real f1 = object->getVertexEnergy1(vid);
							if (f1 < f0)
								x(vid) = xi + alpha * dxi;
						}
						if (x(vid) < 1e-4)
							x(vid) = 1e-4;
					}
					// update positions
					object->getAmpVec() = x;
				}
				// clamp
				object->clampDofs(x);
				// update
				object->getAmpVec() = x;
			}
		}
	}
	void WrinkleSolver::optimize(std::vector<std::unique_ptr<DynamicsObject>>& objects, bool amps, bool dphis, bool positions)
	{
		if (!(amps || dphis || positions))
			throw "Error, enable atleast one sheet parameter to be optimized!\n";
		if (amps)
		{
			// optimizeAmplitudesGlobal(objects);
			// optimizeAmplitudesVBD(objects);
			optimizeAmplitudesNew(objects);
		}

		// initWrinkleMesh(objects, false, true);
		if (dphis)
		{
			// optimizeDphisParallel(objects);
			optimizeDphisGlobal(objects);
		}

		if (positions)
			throw "Static optimizer not implemented for positions!\n";
	}

	void WrinkleSolver::solve(std::vector<std::unique_ptr<DynamicsObject>>& dynamicObjects, std::vector<std::unique_ptr<argus::Obstacle>>& obs, nlohmann::json& state)
	{
		// choose a vertex
		int vid = 0;
		int fid = 0;
		int dphiid = 0;
		Real a0, h0, a1, h1;
		a0 = dynamicObjects.at(0)->getAmpVec()(vid);
		h0 = dynamicObjects.at(0)->getDphis().col(fid)(dphiid);
		optimize(dynamicObjects, true, false, false);
		a1 = dynamicObjects.at(0)->getAmpVec()(vid);
		h1 = dynamicObjects.at(0)->getDphis().col(fid)(dphiid);
		// setup energy computations to visualize contours
		// dynamicObjects.at(0)->writeEnergyContour(vid, fid, dphiid, 1, a0, -1, a1, h0, h1);

		a0 = dynamicObjects.at(0)->getAmpVec()(vid);
		h0 = dynamicObjects.at(0)->getDphis().col(fid)(dphiid);
		optimize(dynamicObjects, false, true, false);
		a1 = dynamicObjects.at(0)->getAmpVec()(vid);
		h1 = dynamicObjects.at(0)->getDphis().col(fid)(dphiid);
		// setup energy computations to visualize contours
		// dynamicObjects.at(0)->writeEnergyContour(vid, fid, dphiid, 2, a0, -1, a1, h0, h1);
	}
	VectorN WrinkleSolver::getAdjacentDphis(int vid, std::unique_ptr<DynamicsObject>& dynamicObject)
	{
		VectorN dphis;
		int counter = 0;
		const Eigen::VectorXi& faces = dynamicObject->getSheetInfo().vertexFacesSet.row(vid);
		std::vector<int> adjfaces;
		for (int i = 0; i < faces.size(); i++)
		{
			if (faces(i) != -1)
				adjfaces.push_back(faces(i));
		}
		int nfaces = adjfaces.size();
		dphis.resize(2 * nfaces);
		for (int i = 0; i < nfaces; i++)
		{
			int fid = adjfaces.at(i);
			dphis(Eigen::seq(2 * i, 2 * i + 1)) = dynamicObject->getDphis().col(fid);
		}
		counter++;
		return dphis;
	}

	void WrinkleSolver::setAdjacentDphis(int vid, const VectorN& dphis, std::unique_ptr<DynamicsObject>& dynamicObject)
	{
		const Eigen::VectorXi& faces = dynamicObject->getSheetInfo().vertexFacesSet.row(vid);
		std::vector<int> adjfaces;
		for (int i = 0; i < faces.size(); i++)
		{
			if (faces(i) != -1)
				adjfaces.push_back(faces(i));
		}

		int nfaces = adjfaces.size();
		if (dphis.size() != nfaces * 2)
			throw "Error! Dphis incompatible!";
		for (int i = 0; i < nfaces; i++)
		{
			int fid = adjfaces.at(i);
			dynamicObject->getDphis().col(fid) = dphis(Eigen::seq(2 * i, 2 * i + 1));
		}
	}

	VectorN WrinkleSolver::getNewDphi(VectorN* gradDphis, Real avals, unsigned int vid, std::unique_ptr<DynamicsObject>& dynamicObject)
	{
		/* def omega(self, a):
		    val = self.x[2] - (1 / self.sheet.Eff(self.x)) * (
		        self.sheet.Ef(self.x) + self.sheet.Eaf(self.x) * (a - self.x[1])
		    )
		    # use optimal value and see which converges faster
		    # ip = deepcopy(self.x)
		    # ip[1] = a
		    # ip[2] = val
		    # copy_activedofs = deepcopy(self.sheet.activeDofs)
		    # self.sheet.activeDofs = [2]
		    # res = optimize.minimize(self.sheet.E, ip, jac=self.sheet.gradE, method="L-BFGS-B")
		    # self.sheet.activeDofs = deepcopy(copy_activedofs)
		    return val

		def gradOmega(self, a):
		    return -self.sheet.Eaf(self.x) / self.sheet.Eff(self.x) */
		VectorN oldDphis = getAdjacentDphis(vid, dynamicObject);
		Real aold = dynamicObject->getAmpVec()(vid);
		std::vector<int> fids;
		VectorN deriv;
		MatrixNN hess;
		dynamicObject->getVertexEnergy1(vid, fids, &deriv, &hess);
		int nadjFaces = dynamicObject->getSheetInfo().n_adjacentFaces(vid);
		MatrixNN Eff;
		VectorN Ef, Eaf, dphi, dphinew, graddphinew;
		Ef = deriv(Eigen::seq(1, 2 * nadjFaces));
		Eaf = hess.col(0)(Eigen::seq(1, Eigen::last));
		Eff = hess.block(1, 1, 2 * nadjFaces, 2 * nadjFaces);
		dphi = oldDphis;
		// compatiblity check
		int Effrows, Effcols, Efsize, Eafsize, dphisize;
		dphisize = dphi.size();
		Effrows = Eff.rows();
		Effcols = Eff.cols();
		Efsize = Ef.size();
		Eafsize = Eaf.size();
		if (Effcols != Efsize)
			throw "Eff incompatible with Ef!";
		if (Effcols != Eafsize)
			throw "Eff incompatible with Eaf!";
		if (Effrows != dphisize)
			throw "Eff incompatible with dphi!";
		dphinew = dphi - Eff.inverse() * (Ef + Eaf * (avals - aold));
		if (gradDphis)
			*gradDphis = -Eff.inverse() * Eaf;
		return dphinew;
		// val = self.x[2] - (1 / self.sheet.Eff(self.x)) * (
		//     self.sheet.Ef(self.x) + self.sheet.Eaf(self.x) * (a - self.x[1])
		// )
	}
}
