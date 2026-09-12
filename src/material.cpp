#include <argus/material.hpp>
#include <argus/sheet.hpp>

namespace argus
{
	Material::Material()
	{
#ifdef ARGUS_DEBUG
		std::cout << "Warning! Add material parameters before startig simulation." << std::endl;
#endif
	}
	Material::Material(const Real density, const std::vector<std::shared_ptr<SheetForces>> forces)
	    : m_Density(density)
	    , m_SheetForces(forces)
	{
#ifdef ARGUS_CHECKPOINT
		simConf["density"] = density;
		simConf["Energies"] = nlohmann::json::array();
		for (auto& force : forces)
		{
			force->saveInfo();
		}
#endif
	}

	Material::Material(const Real density, const Real Y, Real mu, const Real thickness, const std::vector<std::shared_ptr<SheetForces>> forces)
	    : m_Density(density)
	    , m_thickness(thickness)
	    , m_SheetForces(forces)
	    , m_YoungsModulus(Y)
	{
		m_lameAlpha = Y * mu / (1 - mu * mu);
		m_lameBeta = Y / (2 * (1 + mu));
#ifdef ARGUS_CHECKPOINT
		simConf["density"] = density;
		simConf["lame_alpha"] = m_lameAlpha;
		simConf["lame_beta"] = m_lameBeta;
		simConf["Y"] = Y;
		simConf["mu"] = mu;
		simConf["thickness"] = thickness;
		simConf["Energies"] = nlohmann::json::array();
		for (auto& force : forces)
		{
			force->saveInfo();
		}
#endif
	}

	Material::Material(const Material& material)
	    : m_Density(material.m_Density)
	    , m_lameAlpha(material.m_lameAlpha)
	    , m_lameBeta(material.m_lameBeta)
	    , m_thickness(material.m_thickness)
	    , m_SheetForces(material.m_SheetForces)
	{
#ifdef ARGUS_CHECKPOINT
		simConf["density"] = m_Density;
		simConf["lame_alpha"] = m_lameAlpha;
		simConf["lame_beta"] = m_lameBeta;
		simConf["thickness"] = m_thickness;
		simConf["Energies"] = nlohmann::json::array();
#endif
	}

	void Material::setParams(const Real density, const Real Y, const Real mu, const Real thickness)
	{
		m_Density = density;
		m_lameAlpha = Y * mu / (1 - mu * mu);
		m_lameBeta = Y / (2 * (1 + mu));
		m_thickness = thickness;
#ifdef ARGUS_CHECKPOINT
		simConf["density"] = density;
		simConf["lame_alpha"] = m_lameAlpha;
		simConf["lame_beta"] = m_lameBeta;
		simConf["Y"] = Y;
		simConf["mu"] = mu;
		simConf["thickness"] = thickness;
#endif
	}

	void Material::clearSheetForces()
	{
		m_SheetForces.clear();
	}
	void Material::addSheetForces(std::vector<std::shared_ptr<SheetForces>> forces)
	{
		m_SheetForces = forces;
#ifdef ARGUS_CHECKPOINT
		for (auto& force : forces)
		{
			force->saveInfo();
		}
#endif
	}

	void Material::registerRestShapeData(Sheet& sheet)
	{
		for (auto& force : m_SheetForces)
		{
			force->registerRestShapeData(sheet);
		}
	}

	Real Material::getEnergies(Sheet& sheet, VectorN& dofs, VectorN* deriv, std::vector<Eigen::Triplet<Real>>* hess, bool vbd)
	{
		int nV = sheet.getVertexCount();
		int nF = sheet.getFaceCount();
		if (dofs.rows() != 4 * nV + 2 * nF)
			throw "Error! Calling getEnergies over dofs with insufficient size.\n";
		VectorN bak_dofs = sheet.getDofs();
		sheet.getAmplitudes() = dofs(Eigen::seq(0, nV - 1));
		sheet.getPositionsN3() = dofs(Eigen::seq(nV + 2 * nF, 4 * nV + 2 * nF - 1)).reshaped(3, nV);
		sheet.getDphis1VectorPerFace() = dofs(Eigen::seq(nV, nV + 2 * nF - 1)).reshaped(2, nF);
		Real energy = 0;
		if (deriv)
		{
			if (deriv->size() == 0)
				throw "Memory not allocated to derivative variable!";
			deriv->setZero();
		}
		for (auto&& force : m_SheetForces)
		{
			// if vbd iteration is running, then the force must be external
			if (vbd && !force->isExternal())
				continue;

			Real forceEnergy;
			if (deriv && hess)
			{
				std::vector<Eigen::Triplet<Real>> tempHess;
				VectorN tempDeriv;
				force->getEnergy(sheet, forceEnergy, &tempDeriv, &tempHess);
				(*deriv) += tempDeriv;
				hess->insert(hess->end(), tempHess.begin(), tempHess.end());
			}
			else if (deriv && !hess)
			{
				VectorN tempDeriv;
				force->getEnergy(sheet, forceEnergy, &tempDeriv, NULL);
				(*deriv) += tempDeriv;
			}
			else
				force->getEnergy(sheet, forceEnergy, NULL);
			energy += forceEnergy;
		}
		// reset sheet
		sheet.getAmplitudes() = bak_dofs(Eigen::seq(0, nV - 1));
		sheet.getPositionsN3() = bak_dofs(Eigen::seq(nV + 2 * nF, 4 * nV + 2 * nF - 1)).reshaped(3, nV);
		sheet.getDphis1VectorPerFace() = bak_dofs(Eigen::seq(nV, nV + 2 * nF - 1)).reshaped(2, nF);
		return energy;
	}

	Real Material::getEnergies(Sheet& sheet, VectorN* deriv, std::vector<Eigen::Triplet<Real>>* hess, bool vbd)
	{
		Real energy = 0;
		if (deriv)
		{
			if (deriv->size() == 0)
				throw "Memory not allocated to derivative variable!";
			deriv->setZero();
		}
		for (auto&& force : m_SheetForces)
		{
			// if vbd iteration is running, then the force must be external
			if (vbd && !force->isExternal())
				continue;

			Real forceEnergy;
			if (deriv && hess)
			{
				std::vector<Eigen::Triplet<Real>> tempHess;
				VectorN tempDeriv;
				force->getEnergy(sheet, forceEnergy, &tempDeriv, &tempHess);
				(*deriv) += tempDeriv;
				hess->insert(hess->end(), tempHess.begin(), tempHess.end());
			}
			else if (deriv && !hess)
			{
				VectorN tempDeriv;
				force->getEnergy(sheet, forceEnergy, &tempDeriv, NULL);
				(*deriv) += tempDeriv;
			}
			else
				force->getEnergy(sheet, forceEnergy, NULL);
			energy += forceEnergy;
		}
		return energy;
	}
	Real Material::getVertexEnergy1(Sheet& sheet, int vid, std::vector<int>& fids, VectorN* deriv, MatrixNN* hess)
	{
		Real energy = 0;
		int nadjacentFaces = sheet.getSheetInfo().n_adjacentFaces(vid);
		if (deriv)
		{
			deriv->resize(2 * nadjacentFaces + 1);
			deriv->setZero();
		}

		if (hess)
		{
			hess->resize(2 * nadjacentFaces + 1, 2 * nadjacentFaces + 1);
			hess->setZero();
		}
		/* if (deriv && deriv->size() == 0)
		    throw "Error while calling getVertexEnergy1 derivatives not declared!";
		if (hess && hess->size() == 0)
		    throw "Error while calling getVertexEnergy1 hess not declared!"; */
		for (auto&& force : m_SheetForces)
		{
			// if vbd iteration is running, then the force must be external
			if (force->isExternal())
				continue;
			try
			{
				VectorN tempDeriv;
				MatrixNN tempHess;
				energy += force->getVertexEnergy1(sheet, vid, fids, deriv ? &tempDeriv : NULL, hess ? &tempHess : NULL);
				if (deriv)
				{
					(*deriv) += tempDeriv;
				}
				if (hess)
					(*hess) += tempHess;
			}
			catch (const char* err)
			{
				std::cout << "Caught an exception of character const in Material::getVertexEnergy.\n"
				          << err << "\nExiting" << std::endl;
				exit(0);
			}
			catch (...)
			{
				std::cerr << "Caught an unknown exception in Material::getVertexEnergy. Exiting" << std::endl;
				exit(0);
			}
		}
		return energy;
	}
	Real Material::getVertexEnergy1(Sheet& sheet, int vid, Real* deriv, Real* hess)
	{
		Real energy = 0;
		if (deriv)
			*deriv = 0;
		if (hess)
			*hess = 0;
		for (auto&& force : m_SheetForces)
		{
			// if vbd iteration is running, then the force must be external
			if (force->isExternal())
				continue;
			try
			{
				Real tempDeriv;
				Real tempHess;
				energy += force->getVertexEnergy1(sheet, vid, deriv ? &tempDeriv : NULL, hess ? &tempHess : NULL);
				if (deriv)
					(*deriv) += tempDeriv;
				if (hess)
					(*hess) += tempHess;
			}
			catch (const char* err)
			{
				std::cout << "Caught an exception of character const in Material::getVertexEnergy.\n"
				          << err << "\nExiting" << std::endl;
				exit(0);
			}
			catch (...)
			{
				std::cerr << "Caught an unknown exception in Material::getVertexEnergy. Exiting" << std::endl;
				exit(0);
			}
		}
		return energy;
	}
	Real Material::getVertexEnergy1New(Sheet& sheet, int vid, Real* deriv, Real* hess)
	{
		Real energy = 0;
		if (deriv)
			*deriv = 0;
		if (hess)
			*hess = 0;
		for (auto&& force : m_SheetForces)
		{
			// if vbd iteration is running, then the force must be external
			if (force->isExternal())
				continue;
			try
			{
				Real tempDeriv;
				Real tempHess;
				energy += force->getVertexEnergy1New(sheet, vid, deriv ? &tempDeriv : NULL, hess ? &tempHess : NULL);
				if (deriv)
					(*deriv) += tempDeriv;
				if (hess)
					(*hess) += tempHess;
			}
			catch (const char* err)
			{
				std::cout << "Caught an exception of character const in Material::getVertexEnergy.\n"
				          << err << "\nExiting" << std::endl;
				exit(0);
			}
			catch (...)
			{
				std::cerr << "Caught an unknown exception in Material::getVertexEnergy. Exiting" << std::endl;
				exit(0);
			}
		}
		return energy;
	}
	Real Material::getVertexEnergy4(Sheet& sheet, int vid, Vector4* deriv, Matrix44* hess)
	{
		Real energy = 0;
		if (deriv)
			deriv->setZero();
		if (hess)
			hess->setZero();
		for (auto&& force : m_SheetForces)
		{
			// if vbd iteration is running, then the force must be external
			if (force->isExternal())
				continue;
			try
			{
				Vector4 tempDeriv;
				Matrix44 tempHess;
				energy += force->getVertexEnergy4(sheet, vid, deriv ? &tempDeriv : NULL, hess ? &tempHess : NULL);
				if (deriv)
					(*deriv) += tempDeriv;
				if (hess)
					(*hess) += tempHess;
			}
			catch (const char* err)
			{
				std::cout << "Caught an exception of character const in Material::getVertexEnergy.\n"
				          << err << "\nExiting" << std::endl;
				exit(0);
			}
			catch (...)
			{
				std::cerr << "Caught an unknown exception in Material::getVertexEnergy. Exiting" << std::endl;
				exit(0);
			}
		}
		return energy;
	}
	Real Material::getVertexEnergy4New(Sheet& sheet, int vid, Vector4* deriv, Matrix44* hess)
	{
		Real energy = 0;
		if (deriv)
			deriv->setZero();
		if (hess)
			hess->setZero();
		for (auto&& force : m_SheetForces)
		{
			// if vbd iteration is running, then the force must be external
			if (force->isExternal())
				continue;
			// try
			// {
			Vector4 tempDeriv;
			Matrix44 tempHess;
			energy += force->getVertexEnergy4New(sheet, vid, deriv ? &tempDeriv : NULL, hess ? &tempHess : NULL);
			if (deriv)
				(*deriv) += tempDeriv;
			if (hess)
				(*hess) += tempHess;
			// }
			// catch (const char* err)
			// {
			// 	std::cout << "Caught an exception of character const in Material::getVertexEnergy4New.\n"
			// 	          << err << "\nExiting" << std::endl;
			// 	exit(0);
			// }
			// catch (...)
			// {
			// 	std::cerr << "Caught an unknown exception in Material::getVertexEnergy4New. Exiting" << std::endl;
			// 	exit(0);
			// }
		}
		return energy;
	}

	Real Material::getVertexEnergy(Sheet& sheet, int vid, Vector3* deriv, Matrix33* hess)
	{
		Real energy = 0;
		if (deriv)
			deriv->setZero();
		if (hess)
			hess->setZero();
		for (auto&& force : m_SheetForces)
		{
			// if vbd iteration is running, then the force must be external
			if (force->isExternal())
				continue;
			try
			{
				Vector3 tempDeriv;
				Matrix33 tempHess;

				energy += force->getVertexEnergy(sheet, vid, deriv ? &tempDeriv : NULL, hess ? &tempHess : NULL);
				if (deriv)
					(*deriv) += tempDeriv;
				if (hess)
					(*hess) += tempHess;
			}
			catch (const char* err)
			{
				std::cout << "Caught an exception of character const in Material::getVertexEnergy.\n"
				          << err << "\nExiting" << std::endl;
				exit(0);
			}
			catch (...)
			{
				std::cerr << "Caught an unknown exception in Material::getVertexEnergy. Exiting" << std::endl;
				exit(0);
			}
		}
		return energy;
	}
	Real Material::getVertexEnergyNew(Sheet& sheet, int vid, Vector3* deriv, Matrix33* hess)
	{
		Real energy = 0;
		if (deriv)
			deriv->setZero();
		if (hess)
			hess->setZero();
		for (auto&& force : m_SheetForces)
		{
			// if vbd iteration is running, then the force must be external
			if (force->isExternal())
				continue;
			try
			{
				Vector3 tempDeriv;
				Matrix33 tempHess;

				energy += force->getVertexEnergyNew(sheet, vid, deriv ? &tempDeriv : NULL, hess ? &tempHess : NULL);
				if (deriv)
					(*deriv) += tempDeriv;
				if (hess)
					(*hess) += tempHess;
			}
			catch (const char* err)
			{
				std::cout << "Caught an exception of character const in Material::getVertexEnergy.\n"
				          << err << "\nExiting" << std::endl;
				exit(0);
			}
			catch (...)
			{
				std::cerr << "Caught an unknown exception in Material::getVertexEnergy. Exiting" << std::endl;
				exit(0);
			}
		}
		return energy;
	}
	Real Material::getWrinkleShellEnergiesPerFace(Sheet& sheet, int faceId, VectorN* deriv, MatrixNN* hess)
	{
		Real energy = 0;
		if (deriv)
		{
			deriv->resize(3 + 2 + 3 * 3 + 3 * 3);
			deriv->setZero();
		}
		if (hess)
		{
			hess->resize(23, 23);
			hess->setZero();
		}

		for (auto&& force : m_SheetForces)
		{
			// std::cout << "Sheet Force for face " << faceId << std::endl;
			Real forceEnergy;
			if (hess && deriv)
			{
				MatrixNN tempHess;
				VectorN tempDeriv;
				force->getWrinkleShellEnergyPerFace(sheet, forceEnergy, faceId, &tempDeriv, &tempHess);
				(*deriv) += tempDeriv;
				(*hess) += tempHess;
			}
			else if (deriv && !hess)
			{
				// std::cout << "Getting deriv " << std::endl;
				VectorN tempDeriv;
				force->getWrinkleShellEnergyPerFace(sheet, forceEnergy, faceId, &tempDeriv, NULL);
				(*deriv) += tempDeriv;
			}
			else
				force->getWrinkleShellEnergyPerFace(sheet, forceEnergy, faceId, NULL);
			energy += forceEnergy;
		}
		return energy;
	}
	// void Material::optimizeDphis1VectorPerFace(Sheet& sheet, std::vector<Real>& energies,int max_iter)
	// {
	// 	LBFGSpp::LBFGSParam<argus::Real> param;
	// 	param.epsilon = 1e-6;
	// 	param.max_iterations = max_iter;
	// 	int nfaces = sheet.getTriangles().size();
	//
	// }
	// TODO: Remove
	/* 	Real Material::getWrinkleShellEnergies(Sheet& sheet, std::vector<bool> params, VectorN* deriv, std::vector<Eigen::Triplet<Real> >* hess)
	    {
	        Real energy = 0;
	        if(deriv)
	        {
	            if(deriv->size() == 0)
	                throw "Memory not allocated to derivative variable!";
	        }

	        for (auto&& force : m_SheetForces)
	        {
	            Real forceEnergy;
	            if(deriv && hess)
	            {
	                std::vector<Eigen::Triplet<Real> > tempHess;
	                VectorN tempDeriv;
	                force->getEnergy(sheet, forceEnergy, &tempDeriv, &tempHess);
	                (*deriv) += tempDeriv;
	                hess->insert(hess->end(), tempHess.begin(), tempHess.end());
	            }
	            else if(deriv && !hess)
	            {
	                VectorN tempDeriv;
	                force->getEnergy(sheet, forceEnergy, &tempDeriv, NULL);
	                (*deriv) += tempDeriv;
	            }
	            else
	                force->getEnergy(sheet, forceEnergy, NULL);
	            energy += forceEnergy;
	        }
	        return energy;
	    } */
	Matrix3N& Material::getForces(Sheet& sheet,
	    Matrix3N& forces) const
	{
		forces.setZero();

		for (auto&& force : m_SheetForces)
		{
			force->addForces(sheet, forces);
		}

		return forces;
	};
	Eigen::Vector<Real, 23>& Material::getForcesPerFace(Sheet& sheet, Eigen::Vector<Real, 23>& forces, int faceId) const
	{
		forces.setZero();
		for (auto&& force : m_SheetForces)
		{
			force->addForcesPerFace(sheet, forces, faceId);
		}
		return forces;
	}
	Real Material::getEnergyPerFace(Sheet& sheet, Eigen::Vector<Real, 23>* forces, int faceId) const
	{
		Real energy = 0;
		if (forces)
		{
			forces->setZero();
		}
		for (auto&& force : m_SheetForces)
		{
			energy += force->getEnergyPerFace(sheet, forces, faceId);
		}
		return energy;
	}
	VectorN& Material::getForces(Sheet& sheet,
	    VectorN& forces) const
	{
		forces.setZero();

		for (auto&& force : m_SheetForces)
		{
			force->addForces(sheet, forces);
		}

		return forces;
	};
	Real Material::getEnergy(Sheet& sheet) const
	{
		Real energy = 0;
		for (auto&& force : m_SheetForces)
		{
			energy += force->getEnergy(sheet);
		}

		return energy;
	};
#ifdef BOGUS
	void Material::getForcesJacobian(Sheet& sheet, BlockSparseMatrix& forcesPositionDeriv, BlockSparseMatrix& forcesVelocityDeriv)
	{
		forcesPositionDeriv.setBlocksToZero();
		forcesVelocityDeriv.setBlocksToZero();

		for (auto&& force : m_SheetForces)
		{
			force->addForcesJacobian(sheet, forcesPositionDeriv, forcesVelocityDeriv);
		}
	}
#endif

	std::set<std::pair<size_t, size_t>> Material::getJacobianIndices(Sheet& sheet)
	{
		std::set<std::pair<size_t, size_t>> indices;
		for (auto&& force : m_SheetForces)
		{
			force->addJacobianIndices(sheet, indices);
		}
		return indices;
	}

	const void Material::printDebInfo(Sheet& sheet, int faceId, int quadId) const
	{
		for (auto&& force : m_SheetForces)
		{
			force->printDebInfo(sheet, faceId, quadId);
		}
	}

}
