#include <argus/dphi_penalty.hpp>

namespace argus
{
#ifdef ARGUS_CHECKPOINT
	void DphiPenalty::saveInfo()
	{
		nlohmann::json data;
		data["Energy"] = "Slow dphi penalty";
		data["Description"] = "Penalize dphis that become too small";
		data["Is external"]
		    = isExternal();
		data["penalty threshold"] = m_threshold;
		simConf["Energies"].push_back(data);
	}
#endif

	// constructor
	DphiPenalty::DphiPenalty(argus::Real d)
	{
		isExternal() = true;
		m_threshold = d;
	}

	void DphiPenalty::getEnergy(Sheet& sheet, Real& energy, VectorN* deriv, std::vector<Eigen::Triplet<Real>>* hess) const
	{
		energy = 0;
		int nfaces = sheet.getTriangles().size();
		int nvertices = sheet.getPositionsN3().cols();
		std::vector<VectorN> derivatives;
		std::vector<Real> energies(nfaces, 0.0);
		if (deriv)
		{
			derivatives.resize(nfaces);
			for (int i = 0; i < nfaces; i++)
			{
				derivatives[i].resize(3 + 2 + 3 * 3 + 3 * 3);
				derivatives[i].setZero();
			}
		}
#ifdef USE_OMP
#pragma omp parallel for
#endif
		for (int i = 0; i < nfaces; i++)
		{
			// get energies and derivatives from all faces
			argus::Real tempe;
			if (deriv)
				getWrinkleShellEnergyPerFace(sheet, tempe, i, &derivatives.at(i), NULL);
			else
				getWrinkleShellEnergyPerFace(sheet, tempe, i, NULL, NULL);
			energies.at(i) = tempe;
			/*if(deriv)
			{
			    std::cout << "received derivatives " << std::endl;
			    std::cout << derivatives.at(i) << std::endl;
			    exit(0);
			} */
		}

		// consolidate energies
		for (int i = 0; i < nfaces; i++)
		{
			energy += energies.at(i);
		}
		if (deriv)
		{
			deriv->resize(nvertices + 2 * nfaces + 3 * nvertices);
			deriv->setZero();
			// consolidate derivatives
			for (int i = 0; i < nfaces; i++)
			{
				deriv->coeffRef(nvertices + 2 * i) += derivatives.at(i)[3];
				deriv->coeffRef(nvertices + 2 * i + 1) += derivatives.at(i)[4];
			}
		}
		/* if(deriv)
		{
		    std::cout << "Returning final energy " << energy << std::endl;
		    std::cout << "Returning final derivatives " << *deriv << std::endl;
		    exit(0);
		} */
	}

	void DphiPenalty::getWrinkleShellEnergy(Sheet& sheet, Real& energy, VectorN* deriv, std::vector<Eigen::Triplet<Real>>* hess, std::vector<bool>* params) const
	{
		if (deriv)
			getEnergy(sheet, energy, deriv, NULL);
		else
			getEnergy(sheet, energy, NULL, NULL);
	}

	void DphiPenalty::getWrinkleShellEnergyPerFace(Sheet& sheet, Real& energy, int faceId, VectorN* deriv, MatrixNN* hess) const
	{
		energy = 0;
		VectorN derivative;
		if (deriv)
		{
			derivative.resize(2);
			derivative.setZero();
		}
		// get per face penalty energy
		argus::Vector2 dphi = sheet.getDphis1VectorPerFace().col(faceId);
		const argus::Matrix32& M = sheet.getDphiToGradPhis().at(faceId);
		Real delta = sheet.getMinimumEdgeLengths().at(faceId) * m_threshold;
		argus::Vector3 gradphi = M * dphi;
		if (gradphi.norm() < delta)
		{
			// std::cout << "Dphis too low " << std::endl;
			argus::Real penalty = gradphi.norm() - delta;
			energy = energy + penalty * penalty; // std::pow(penalty, 2);
			if (deriv)
			{
				argus::Vector3 grad_i = M * argus::Vector2({ 1, 0 });
				argus::Vector3 grad_j = M * argus::Vector2({ 0, 1 });
				// handle zero grad phi
				if (gradphi.norm() == 0)
					gradphi = argus::Vector3::Ones() * 1e-4;
				derivative(0) = 2 * penalty * (gradphi / gradphi.norm()).dot(grad_i);
				derivative(1) = 2 * penalty * (gradphi / gradphi.norm()).dot(grad_j);
				// std::cout << "Derivative returning " << derivative.transpose() << std::endl;
			}
		}

		// assuming that this function is used only for dphi per face optimization
		if (deriv)
		{
			deriv->resize(3 + 2 + 3 * 3 + 3 * 3);
			deriv->setZero();
			(*deriv)(Eigen::seq(3, 4)) += derivative;
		}
	}

	Real DphiPenalty::getEnergyPerFace(Sheet& sheet, Eigen::Vector<Real, 23>* forces, int faceId) const
	{
		Real energy = 0;
		VectorN derivative;
		if (forces)
		{
			derivative.resize(2);
			derivative.setZero();
		}
		// get per face penalty energy
		argus::Vector2 dphi = sheet.getDphis1VectorPerFace().col(faceId);
		const argus::Matrix32& M = sheet.getDphiToGradPhis().at(faceId);
		Real delta = sheet.getMinimumEdgeLengths().at(faceId) * m_threshold;
		argus::Vector3 gradphi = M * dphi;
		if (gradphi.norm() < delta)
		{
			// std::cout << "Dphis too low " << std::endl;
			argus::Real penalty = gradphi.norm() - delta;
			energy = energy + penalty * penalty; // std::pow(penalty, 2);
			if (forces)
			{
				argus::Vector3 grad_i = M * argus::Vector2({ 1, 0 });
				argus::Vector3 grad_j = M * argus::Vector2({ 0, 1 });
				// handle zero grad phi
				if (gradphi.norm() == 0)
					gradphi = argus::Vector3::Ones() * 1e-4;
				derivative(0) = 2 * penalty * (gradphi / gradphi.norm()).dot(grad_i);
				derivative(1) = 2 * penalty * (gradphi / gradphi.norm()).dot(grad_j);
				// std::cout << "Derivative returning " << derivative.transpose() << std::endl;
			}
		}

		// assuming that this function is used only for dphi per face optimization
		if (forces)
		{
			forces->setZero();
			(*forces)(Eigen::seq(3, 4)) += derivative;
		}
		return energy;
	}

	void DphiPenalty::addForces(Sheet& sheet, Matrix3N& forces) const
	{
		throw "Not Implemented!";
	}
	void DphiPenalty::addForces(Sheet& sheet, VectorN& forces) const
	{
		throw "Not Implemented!";
	}
};

// 	void Sheet::applySmallDphiConstraints()
// 	{
// #ifdef USE_OMP
// #pragma omp parallel for
// #endif
// 		for (int faceId = 0; faceId < m_FaceCount; faceId++)
// 		{
// 			argus::Vector2 dphi = m_dPhis1VectorPerFace.col(faceId);
// 			const argus::Matrix32& M = m_dphi_to_gradphis.at(faceId);
// 			Real delta = m_minimum_edge_lengths.at(faceId) * PI / 2.0;
// 			argus::Vector3 gradphi = M * dphi;
// 			if (gradphi.norm() < delta)
// 			{
// 				m_dPhis1VectorPerFace.col(faceId) *= 1e-4;
// 			}
// 		}
// 	}