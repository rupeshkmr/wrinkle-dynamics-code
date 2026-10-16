#ifndef SHEET_FORCES_HPP
#define SHEET_FORCES_HPP

#include "common.hpp"

#include <set>

namespace argus
{
	class Sheet;
	class SheetForces
	{
	protected:
		bool external;

	public:
		bool& isExternal() { return external; };
		// for cloth models with poisitons are dofs
#ifdef ARGUS_CHECKPOINT
		virtual void saveInfo() { throw "Not Implemented! Add this function to ensure proper checkpointing."; }
#endif
		virtual void getEnergy(Sheet& sheet, Real& energy, VectorN* deriv = NULL, std::vector<Eigen::Triplet<Real>>* hess = NULL) const { throw "Not implemented"; }
		virtual Real getVertexEnergy(Sheet& sheet, int vid, Vector3* deriv = NULL, Matrix33* hess = NULL) const { throw "Not implemented"; }
		virtual Real getVertexEnergyNew(Sheet& sheet, int vid, Vector3* deriv = NULL, Matrix33* hess = NULL) const { throw "Not implemented"; }
		virtual Real getVertexEnergy4(Sheet& sheet, int vid, Vector4* deriv = NULL, Matrix44* hess = NULL) const { throw "Not implemented"; }
		virtual Real getVertexEnergy4New(Sheet& sheet, int vid, Vector4* deriv = NULL, Matrix44* hess = NULL) const { throw "Not implemented"; }
		// function to get forces on amplitudes only
		virtual Real getVertexEnergy1(Sheet& sheet, int vid, Real* deriv = NULL, Real* hess = NULL) const { throw "Not implemented"; }
		virtual Real getVertexEnergy1New(Sheet& sheet, int vid, Real* deriv = NULL, Real* hess = NULL) const { throw "Not implemented"; }

		// 0th index for deriv and hess are for amps
		virtual Real getVertexEnergy1(Sheet& sheet, int vid, std::vector<int>& fids, VectorN* deriv = NULL, MatrixNN* hess = NULL) const { throw "Not implemented"; }
		// for cloth models with wrinkle parameters
		// TODO: Remove
		virtual void getWrinkleShellEnergy(Sheet& sheet, Real& energy, std::vector<bool> params, VectorN* deriv = NULL, std::vector<Eigen::Triplet<Real>>* hess = NULL) const { throw "Not implemented"; }
		virtual void getWrinkleShellEnergyPerFace(Sheet& sheet, Real& energy, int faceId, VectorN* deriv = NULL, MatrixNN* hess = NULL) const { throw "Not implemented"; }
		// rest shape information initialization
		virtual void registerRestShapeData(Sheet& sheet) { }
		// get sheet forces
		virtual void addForces(Sheet& sheet, Matrix3N& forces) const = 0;
		virtual void addForces(Sheet& sheet, VectorN& forces) const { throw "Not Implemented"; }
		virtual Real getEnergyPerFace(Sheet& sheet, Eigen::Vector<Real, 23>* forces, int faceId) const { throw "Not Implemented"; }
		virtual void addForcesPerFace(Sheet& sheet, Eigen::Vector<Real, 23>& forces, int faceId) const { throw "Not Implemented"; }
		virtual Real getEnergy(Sheet& sheet) const { throw "Not Implemented! Add this function to ensure proper checkpointing."; }
		// Jacobian matrices are symmetric so we need only one pair of indices
		virtual void addJacobianIndices(Sheet& sheet, std::set<std::pair<size_t, size_t>>& indices) const { throw "Not implemented"; }
// This matrix is 3N x 3N
#ifdef BOGUS
		virtual void addForcesJacobian(Sheet& sheet,
		    BlockSparseMatrix& forcesPositionDeriv,
		    BlockSparseMatrix& forcesVelocityDeriv) const { throw "Not implemented"; }
#endif
		virtual const void printDebInfo(Sheet& sheet, int faceId, int quadId) const { throw "Not Implemented!"; }
	};
}

#endif /* SHEET_FORCES_HPP */
