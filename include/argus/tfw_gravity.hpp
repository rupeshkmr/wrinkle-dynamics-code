#ifndef ARGUS_TFW_GRAVITY_HPP
#define ARGUS_TFW_GRAVITY_HPP
#include "common.hpp"
#include "sheet.hpp"
namespace argus
{
	class TFWGravity : public SheetForces
	{

		bool amp_enabled;
		Vector3 g;

	public:
		TFWGravity(Vector3 gravitationalCoefficient, bool amp_e = false)
		{
			g = gravitationalCoefficient;
			amp_enabled = amp_e;
			isExternal() = true;
		}
#ifdef ARGUS_CHECKPOINT
		void saveInfo() override;
#endif
		void getEnergy(Sheet& sheet, Real& energy, VectorN* deriv = NULL, std::vector<Eigen::Triplet<Real>>* hess = NULL) const override;
		void getWrinkleShellEnergyPerFace(Sheet& sheet, Real& energy, int faceId, VectorN* deriv = NULL, MatrixNN* hess = NULL) const override;
		void addForces(Sheet& sheet,
		    Matrix3N& forces) const override;
		// Real getEnergy(Sheet& sheet) const override;
		void addForcesPerFace(Sheet& sheet, Eigen::Vector<Real, 23>& forces, int faceId) const override;
		Real getEnergyPerFace(Sheet& sheet, Eigen::Vector<Real, 23>* forces, int faceId) const override { return 0.0; }
		Real getEnergy(Sheet& sheet) const override;
		void addForces(Sheet& sheet, VectorN& forces) const override;

		void addJacobianIndices(Sheet& sheet, std::set<std::pair<size_t, size_t>>& indices) const override { }
		// This matrix is 3N x 3N
#ifdef BOGUS
		void addForcesJacobian(Sheet& sheet,
		    BlockSparseMatrix& forcesPositionDeriv,
		    BlockSparseMatrix& forcesVelocityDeriv) const override { }
#endif
	};
}
#endif
