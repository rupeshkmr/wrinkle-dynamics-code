#ifndef GEOMETRY_FUNCTIONS_HPP
#define GEOMETRY_FUNCTIONS_HPP

#include "common.hpp"
#include "sheet.hpp"
namespace argus
{

	inline Real doubleContraction(const Matrix22& M1, const Matrix22& M2)
	{
		return M1.cwiseProduct(M2.transpose()).sum();

		// return (M1 * M2).trace();
	}
	inline Real cotan(const Vector3& v0, const Vector3& v1, const Vector3& v2)
	{
		Real e0 = (v2 - v1).norm();
		Real e1 = (v2 - v0).norm();
		Real e2 = (v0 - v1).norm();
		Real angle0 = acos((e1 * e1 + e2 * e2 - e0 * e0) / (2 * e1 * e2));
		Real cot = 1.0 / tan(angle0);

		return cot;
	}

	inline Real d2Edqiqj(const Matrix22& T1, const Matrix22& T2, const Matrix22& T3, const Matrix22& T4)
	{
		return (T1.cwiseProduct(T2)).sum() + (T3.cwiseProduct(T4)).sum();
	}

	inline Real SVNormSquared(Real alpha, Real beta, const Matrix22& M, Matrix22* deriv = NULL)
	{
		Real Mtrace = M.trace();
		if (deriv)
		{
			// (*deriv) = alpha * M.trace() * Matrix22::Identity() + 2 * beta * M;
			(*deriv) = 2 * beta * M;
			(*deriv).diagonal().noalias() += alpha * Mtrace * Vector2::Ones();
		}
		return (alpha / 2.0) * Mtrace * Mtrace + beta * (M * M).trace();
	}
	// void loadActiveVector3(Eigen::Vector<Real,3>& ip, Eigen::Vector<AReal15,3> *op);
	// Real d2Edqiqj(Matrix22 T1, Matrix22 T2, Matrix22 T3, Matrix22 T4);
	// void firstFundamentalForm(
	//     Eigen::Vector3<AReal15>& q0,
	//     Eigen::Vector3<AReal15>& q1,
	//     Eigen::Vector3<AReal15>& q2, Eigen::Matrix2<AReal15>* I);
	Matrix22 firstFundamentalForm(
	    Matrix33& vertexPositions,
	    bool requireDerivative,
	    bool requireHessian,
	    MatrixNN* derivative = NULL, // F(face, i)
	    std::vector<Eigen::Matrix<Real, 9, 9>>* hessian = NULL);
	Eigen::Matrix2d firstFundamentalForm(
	    const MeshConnectivity& mesh,
	    const Eigen::MatrixXd& curPos,
	    int face,
	    Eigen::Matrix<double, 4, 9>* derivative, // F(face, i)
	    std::vector<Eigen::Matrix<double, 9, 9>>* hessian);
	Matrix22 firstFundamentalFormNew(
	    const Matrix33& vertexPositions,
	    bool requireDerivative = false,
	    bool requireHessian = false,
	    Eigen::Matrix<Real, 4, 9>* derivative = NULL, // F(face, i)
	    std::array<Eigen::Matrix<Real, 9, 9>, 4>* hessian = NULL);
	Matrix22 firstFundamentalFormNew(
	    const Matrix33& vertexPositions,
	    Eigen::Matrix<Real, 4, 9>* derivative = NULL,
	    std::array<Eigen::Matrix<Real, 9, 9>, 4>* hessian = NULL);
	// Real SVNormSquared(Real alpha, Real beta, Matrix22& M, Matrix22* deriv = NULL);
	// Real SVNormSquared(Real alpha, Real beta, Eigen::Matrix<AReal4, 2, 2> M);

	Matrix22 secondFundamentalForm(
	    int face,
	    MatrixNN* derivative,
	    std::vector<MatrixNN>* hessian);

	// Real cotan(const Vector3 v0, const Vector3 v1, const Vector3 v2);
}

#endif
