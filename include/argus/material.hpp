#ifndef MATERIAL_HPP
#define MATERIAL_HPP

// #include "../lbfgspp/LBFGS.h"
// #include "../lbfgspp/LBFGSB.h"
#include "common.hpp"

#include "sheet_forces.hpp"

#include <memory>
#include <set>

namespace argus
{
	class Material
	{
		Real m_Density;
		Real m_lameAlpha;
		Real m_lameBeta;
		Real m_thickness;
		Real m_YoungsModulus;
		std::vector<std::shared_ptr<SheetForces>> m_SheetForces;

	public:
		Material();
		Material(const Real density, const std::vector<std::shared_ptr<SheetForces>> forces = {});
		Material(const Real density, const Real Y = 0.0, const Real mu = 0.0, const Real thickness = 0.1, const std::vector<std::shared_ptr<SheetForces>> forces = {});
		Material(const Material& material);

		void setParams(const Real density, const Real Y = 0.0, const Real mu = 0.0, const Real thickness = 0.0);
		Real getEnergies(Sheet& sheet, VectorN* deriv = NULL, std::vector<Eigen::Triplet<Real>>* hess = NULL, bool vbd = false);
		Real getEnergies(Sheet& sheet, VectorN& dofs, VectorN* deriv = NULL, std::vector<Eigen::Triplet<Real>>* hess = NULL, bool vbd = false);
		Real getVertexEnergy(Sheet& sheet, int vid, Vector3* deriv = NULL, Matrix33* hess = NULL);
		Real getVertexEnergyNew(Sheet& sheet, int vid, Vector3* deriv = NULL, Matrix33* hess = NULL);
		Real getVertexEnergy4New(Sheet& sheet, int vid, Vector4* deriv = NULL, Matrix44* hess = NULL);

		Real getVertexEnergy4(Sheet& sheet, int vid, Vector4* deriv = NULL, Matrix44* hess = NULL);
		Real getVertexEnergy1(Sheet& sheet, int vid, Real* deriv = NULL, Real* hess = NULL);
		Real getVertexEnergy1New(Sheet& sheet, int vid, Real* deriv = NULL, Real* hess = NULL);
		Real getEnergy(Sheet& sheet) const;
		Real getVertexEnergy1(Sheet& sheet, int vid, std::vector<int>& fids, VectorN* deriv = NULL, MatrixNN* hess = NULL);
		// todo delete moved to sheet class
		// void optimizeDphis1VectorPerFace(Sheet& sheet, std::vector<Real>& energies, int max_iter= 1000);
		// TODO: Remove
		// Real getWrinkleShellEnergies(Sheet& sheet, std::vector<bool> params, VectorN* deriv=NULL, std::vector<Eigen::Triplet<Real> >* hess=NULL);
		Real getWrinkleShellEnergiesPerFace(Sheet& sheet, int faceId, VectorN* deriv = NULL, MatrixNN* hess = NULL);
		Matrix3N& getForces(Sheet& sheet, Matrix3N& forces) const;
		Eigen::Vector<Real, 23>& getForcesPerFace(Sheet& sheet, Eigen::Vector<Real, 23>& forces, int faceId) const;
		Real getEnergyPerFace(Sheet& sheet, Eigen::Vector<Real, 23>* forces, int faceId) const;

		VectorN& getForces(Sheet& sheet, VectorN& forces) const;
#ifdef BOGUS
		void getForcesJacobian(Sheet& sheet, BlockSparseMatrix& forcesPositionDeriv, BlockSparseMatrix& forcesVelocityDeriv);
#endif
		void registerRestShapeData(Sheet& sheet);
		std::set<std::pair<size_t, size_t>> getJacobianIndices(Sheet& sheet);
		void clearSheetForces();
		void addSheetForces(std::vector<std::shared_ptr<SheetForces>> forces);
		const Real getDensity() const { return m_Density; }
		const Real getLameAlpha() { return m_lameAlpha; }
		const Real getLameBeta() { return m_lameBeta; }
		const Real getThickness() { return m_thickness; }
		const Real getY() { return m_YoungsModulus; }
		const void printDebInfo(Sheet& sheet, int faceId, int quadId) const;
	};
}

#endif /* MATERIAL_HPP */
