#ifndef ARGUS_DPHI_PENALTY_HPP
#define ARGUS_DPHI_PENALTY_HPP
#include "common.hpp"
#include "sheet.hpp"
#include "geometry_functions.hpp"
#include "../../vendor/MeshLib/MeshConnectivity.h"
#include "../../vendor/MeshLib/MeshGeometry.h"
#include "../../vendor/SecondFundamentalForm/SecondFundamentalFormDiscretization.h"
#include "../../vendor/SecondFundamentalForm/MidedgeAverageFormulation.h"
#include "../vendor/MeshLib/IntrinsicGeometry.h"
namespace argus
{
	class DphiPenalty : public SheetForces
	{
		/**
		 * Penalize slow wrinkle frequencies.
		 */
		Real m_threshold;

	public:
		DphiPenalty(argus::Real d);
#ifdef ARGUS_CHECKPOINT
		void saveInfo();
#endif
		void getEnergy(Sheet& sheet, Real& energy, VectorN* deriv = NULL, std::vector<Eigen::Triplet<Real>>* hess = NULL) const override;
		void getWrinkleShellEnergy(Sheet& sheet, Real& energy, VectorN* deriv = NULL, std::vector<Eigen::Triplet<Real>>* hess = NULL, std::vector<bool>* params = NULL) const;
		void getWrinkleShellEnergyPerFace(Sheet& sheet, Real& energy, int faceId, VectorN* deriv = NULL, MatrixNN* hess = NULL) const override;

		void addForces(Sheet& sheet,
		    Matrix3N& forces) const override;
		void addForces(Sheet& sheet, VectorN& forces) const override;

		void addJacobianIndices(Sheet& sheet, std::set<std::pair<size_t, size_t>>& indices) const override { }
		// This matrix is 3N x 3N
		Real getEnergyPerFace(Sheet& sheet, Eigen::Vector<Real, 23>* forces, int faceId) const override;

#ifdef BOGUS
		void addForcesJacobian(Sheet& sheet,
		    BlockSparseMatrix& forcesPositionDeriv,
		    BlockSparseMatrix& forcesVelocityDeriv) const override { }
#endif
	};
}
#endif
