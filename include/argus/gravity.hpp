#ifndef GRAVITY_HPP
#define GRAVITY_HPP

#include "common.hpp"
#include "sheet_forces.hpp"

namespace argus
{
	class Gravity : public SheetForces
	{
		Vector3 m_Gravity;

	public:
#ifdef ARGUS_CHECKPOINT
		void saveInfo() override;
#endif

		Gravity(const Vector3& Gravity)
		    : m_Gravity(Gravity) { };
		void getEnergy(Sheet& sheet, Real& energy, VectorN* deriv = NULL, std::vector<Eigen::Triplet<Real>>* hess = NULL) const override; // returns gravitational potential
		void addForces(Sheet& sheet,
		    Matrix3N& forces) const override;
		Real getEnergy(Sheet& sheet) const override;

		void addJacobianIndices(Sheet& sheet, std::set<std::pair<size_t, size_t>>& indices) const override { }
		// This matrix is 3N x 3N

#ifdef BOGUS
		void addForcesJacobian(Sheet& sheet,
		    BlockSparseMatrix& forcesPositionDeriv,
		    BlockSparseMatrix& forcesVelocityDeriv) const override { }
#endif
	};
}

#endif /* GRAVITY_HPP */
