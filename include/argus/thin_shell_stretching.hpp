#ifndef ARGUS_THIN_SHELL_STRETCHING_HPP
#define ARGUS_THIN_SHELL_STRETCHING_HPP
#include "common.hpp"
#include "sheet.hpp"
/* #include "geometry_functions.hpp"
#include "../../vendor/MeshLib/MeshConnectivity.h"
#include "../../vendor/MeshLib/MeshGeometry.h"
#include "../../vendor/SecondFundamentalForm/SecondFundamentalFormDiscretization.h"
#include "../../vendor/SecondFundamentalForm/MidedgeAverageFormulation.h"
#include "../vendor/MeshLib/IntrinsicGeometry.h" */
namespace argus
{
	class ThinShellStretching : public SheetForces
	{

	public:
		ThinShellStretching()
		{
			isExternal() = false;
		}
#ifdef ARGUS_CHECKPOINT
		void saveInfo() override;
#endif
		void getEnergy(Sheet& sheet, Real& energy, VectorN* deriv = NULL, std::vector<Eigen::Triplet<Real>>* hess = NULL) const override;
		void getWrinkleShellEnergy(Sheet& sheet, Real& energy, VectorN* deriv = NULL, std::vector<Eigen::Triplet<Real>>* hess = NULL) const;

		void getWrinkleShellEnergyPerFace(Sheet& sheet, Real& energy, int faceId, VectorN* deriv, MatrixNN* hess) const override;
		void addForces(Sheet& sheet,
		    Matrix3N& forces) const override;

		void addJacobianIndices(Sheet& sheet, std::set<std::pair<size_t, size_t>>& indices) const override { }
		// This matrix is 3N x 3N

#ifdef BOGUS
		void addForcesJacobian(Sheet& sheet,
		    BlockSparseMatrix& forcesPositionDeriv,
		    BlockSparseMatrix& forcesVelocityDeriv) const override { }
#endif
		// amps per vertex dphis 1 vector per face
		Real getStretchingEnergyDensityFromQuadAmpsPerVtexDphis1VectorPerFace(int faceId, int quadId, Sheet& sheet, VectorN* deriv, MatrixNN* hess) const;
		Real getStretchingEnergyAmpsPerVtexDphis1VectorPerFace(Sheet& sheet, VectorN* deriv, std::vector<Eigen::Triplet<Real>>* hess) const;
		Real getEnergyDensityFromQuadPerVtex(int vid, int faceId, int quadId, Sheet& sheet, Vector3* deriv, Matrix33* hess) const;
		Real getEnergyDensityFromQuadPerVtex4(int vid, int faceId, int quadId, Sheet& sheet, Vector4* deriv, Matrix44* hess) const;
		Real getEnergyDensityFromQuadPerVtex1(int vid, int faceId, int quadId, Sheet& sheet, Real* deriv, Real* hess) const;
		Real getVertexEnergy(Sheet& sheet, int vid, Vector3* deriv, Matrix33* hess) const override;
		Real getVertexEnergy4(Sheet& sheet, int vid, Vector4* deriv, Matrix44* hess) const override;
		Real getVertexEnergy1(Sheet& sheet, int vid, Real* deriv, Real* hess) const override;
		Real getEnergy(Sheet& sheet) const override;
	};

}
#endif
