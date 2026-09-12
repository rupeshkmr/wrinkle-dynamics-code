#ifndef ARGUS_TFW_BENDING_HPP
#define ARGUS_TFW_BENDING_HPP
#include "common.hpp"
#include "sheet.hpp"
namespace argus
{
	// get energy amps per vertex, dphis 1 vector per face
	Real getBendingEnergyDensityFromQuadAmpsPerVtexDphis1VectorPerFace(int faceId, int quadId, Sheet& sheet, VectorN* deriv, MatrixNN* hess);
	Real getBendingEnergyAmpsPerVtexDphis1VectorPerFace(Sheet& sheet, VectorN* deriv, std::vector<Eigen::Triplet<Real>>* hess);
	class TFWBending : public SheetForces
	{
		bool isWrinkleEnergy;

	public:
		TFWBending()
		{
			isWrinkleEnergy = true;
			isExternal() = false;
		}

#ifdef ARGUS_CHECKPOINT
		void saveInfo() override;
#endif
		void getEnergy(Sheet& sheet, Real& energy, VectorN* deriv = NULL, std::vector<Eigen::Triplet<Real>>* hess = NULL) const override;
		void getWrinkleShellEnergy(Sheet& sheet, Real& energy, std::vector<bool> params, VectorN* deriv = NULL, std::vector<Eigen::Triplet<Real>>* hess = NULL) const override;
		void getWrinkleShellEnergyPerFace(Sheet& sheet, Real& energy, int faceId, VectorN* deriv = NULL, MatrixNN* hess = NULL) const override;
		void addForces(Sheet& sheet,
		    Matrix3N& forces) const override;

		Real getEnergyDensityFromQuadPerVtex(int vid, int faceId, int quadId, Sheet& sheet, Vector3* deriv, Matrix33* hess) const;
		Real getVertexEnergy(Sheet& sheet, int vid, Vector3* deriv, Matrix33* hess) const override;
		Real getEnergyDensityFromQuadPerVtex4(int vid, int faceId, int quadId, Sheet& sheet, Vector4* deriv, Matrix44* hess) const;
		Real getVertexEnergy4(Sheet& sheet, int vid, Vector4* deriv, Matrix44* hess) const override;
		Real getVertexEnergy4New(Sheet& sheet, int vid, Vector4* deriv, Matrix44* hess) const override;
		Real getEnergyDensityFromQuadPerVtex1(int vid, int faceId, int quadId, Sheet& sheet, Real* deriv, Real* hess) const;
		Real getVertexEnergy1(Sheet& sheet, int vid, Real* deriv, Real* hess) const override;
		Real getVertexEnergy1(Sheet& sheet, int vid, std::vector<int>& fids, VectorN* deriv, MatrixNN* hess) const override;
		Real getEnergyDensityFromQuadPerVtex1(int vid, int faceId, int quadId, Sheet& sheet, VectorN* deriv, MatrixNN* hess) const;
		Real getVertexEnergy1New(Sheet& sheet, int vid, Real* deriv, Real* hess) const override;
		Real getVertexEnergyNew(Sheet& sheet, int vid, Vector3* deriv, Matrix33* hess) const override;

		// Real getEnergyDensityFromQuadPerVtex(int vid, int faceId, int quadId, Sheet& sheet, Vector3& deriv, Matrix33& hess) const;
		// Real getVBDEnergy(Sheet& sheet, int vid, Vector3* deriv, Matrix33* hess);
		// void addJacobianIndices(Sheet& sheet, std::set<std::pair<size_t, size_t>>& indices) const override { }
		// // This matrix is 3N x 3N
		// void addForcesJacobian(Sheet& sheet,
		//     BlockSparseMatrix& forcesPositionDeriv,
		//     BlockSparseMatrix& forcesVelocityDeriv) const override { }
		//
		void addForcesPerFace(Sheet& sheet, Eigen::Vector<Real, 23>& forces, int faceId) const override;
		Real getEnergyPerFace(Sheet& sheet, Eigen::Vector<Real, 23>* forces, int faceId) const override;
		void addForces(Sheet& sheet, VectorN& forces) const override;
	};
}
#endif
