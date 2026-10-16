#include <argus/geometry_functions.hpp>
#include <argus/common.hpp>
namespace argus
{

	// Real SVNormSquared(Real alpha, Real beta, Eigen::Matrix<AReal4, 2, 2> M)
	// {
	// 	return ((alpha / 2.0) * M.trace() * M.trace() + beta * (M * M).trace()).val;
	// }

	// void loadActiveVector3(Eigen::Vector<Real,3>& ip, Eigen::Vector<AReal15,3> *op)
	// {

	// }
	// void firstFundamentalForm(
	//     Eigen::Vector3<AReal15>& q0,
	//     Eigen::Vector3<AReal15>& q1,
	//     Eigen::Vector3<AReal15>& q2, Eigen::Matrix2<AReal15>* I)
	// {
	// 	// (*I)(0,0) = (q1-q0).dot(q1-q0);
	// 	(*I) << (q1 - q0).dot(q1 - q0), (q1 - q0).dot(q2 - q0),
	// 	    (q2 - q0).dot(q1 - q0), (q2 - q0).dot(q2 - q0);
	// }

	Matrix22 firstFundamentalForm(
	    Matrix33& vertexPositions,
	    bool requireDerivative,
	    bool requireHessian,
	    MatrixNN* derivative, // F(face, i)
	    std::vector<Eigen::Matrix<Real, 9, 9>>* hessian)
	{
		Vector3 q0 = vertexPositions.col(0);
		Vector3 q1 = vertexPositions.col(1);
		Vector3 q2 = vertexPositions.col(2);
		Matrix22 result;
		result << (q1 - q0).dot(q1 - q0), (q1 - q0).dot(q2 - q0),
		    (q2 - q0).dot(q1 - q0), (q2 - q0).dot(q2 - q0);

		if (requireDerivative)
		{
			// std::cout<< "derivative" << std::endl;
			derivative->resize(4, 9);
			derivative->setZero();
			derivative->block<1, 3>(0, 3) += 2.0 * (q1 - q0).transpose();
			derivative->block<1, 3>(0, 0) -= 2.0 * (q1 - q0).transpose();
			derivative->block<1, 3>(1, 6) += (q1 - q0).transpose();
			derivative->block<1, 3>(1, 3) += (q2 - q0).transpose();
			derivative->block<1, 3>(1, 0) += -(q1 - q0).transpose() - (q2 - q0).transpose();
			derivative->block<1, 3>(2, 6) += (q1 - q0).transpose();
			derivative->block<1, 3>(2, 3) += (q2 - q0).transpose();
			derivative->block<1, 3>(2, 0) += -(q1 - q0).transpose() - (q2 - q0).transpose();
			derivative->block<1, 3>(3, 6) += 2.0 * (q2 - q0).transpose();
			derivative->block<1, 3>(3, 0) -= 2.0 * (q2 - q0).transpose();
		}

		if (requireHessian)
		{
			hessian->resize(4);
			for (int i = 0; i < 4; i++)
			{
				(*hessian)[i].resize(9, 9);
				(*hessian)[i].setZero();
			}
			// std::cout << "hessian" << std::endl;
			Matrix33 I = Matrix33::Identity();
			(*hessian)[0].block<3, 3>(0, 0) += 2.0 * I;
			(*hessian)[0].block<3, 3>(3, 3) += 2.0 * I;
			(*hessian)[0].block<3, 3>(0, 3) -= 2.0 * I;
			(*hessian)[0].block<3, 3>(3, 0) -= 2.0 * I;

			(*hessian)[1].block<3, 3>(3, 6) += I;
			(*hessian)[1].block<3, 3>(6, 3) += I;
			(*hessian)[1].block<3, 3>(0, 3) -= I;
			(*hessian)[1].block<3, 3>(0, 6) -= I;
			(*hessian)[1].block<3, 3>(3, 0) -= I;
			(*hessian)[1].block<3, 3>(6, 0) -= I;
			(*hessian)[1].block<3, 3>(0, 0) += 2.0 * I;

			(*hessian)[2].block<3, 3>(3, 6) += I;
			(*hessian)[2].block<3, 3>(6, 3) += I;
			(*hessian)[2].block<3, 3>(0, 3) -= I;
			(*hessian)[2].block<3, 3>(0, 6) -= I;
			(*hessian)[2].block<3, 3>(3, 0) -= I;
			(*hessian)[2].block<3, 3>(6, 0) -= I;
			(*hessian)[2].block<3, 3>(0, 0) += 2.0 * I;

			(*hessian)[3].block<3, 3>(0, 0) += 2.0 * I;
			(*hessian)[3].block<3, 3>(6, 6) += 2.0 * I;
			(*hessian)[3].block<3, 3>(0, 6) -= 2.0 * I;
			(*hessian)[3].block<3, 3>(6, 0) -= 2.0 * I;
		}

		return result;
	}
	Matrix22 firstFundamentalFormNew(
	    const Matrix33& vertexPositions,
	    Eigen::Matrix<Real, 4, 9>* derivative,
	    std::array<Eigen::Matrix<Real, 9, 9>, 4>* hessian)
	{
		Vector3 q0 = vertexPositions.col(0);
		Vector3 q1 = vertexPositions.col(1);
		Vector3 q2 = vertexPositions.col(2);
		Matrix22 result;
		result << (q1 - q0).dot(q1 - q0), (q1 - q0).dot(q2 - q0),
		    (q2 - q0).dot(q1 - q0), (q2 - q0).dot(q2 - q0);

		if (derivative)
		{
			// std::cout<< "derivative" << std::endl;
			derivative->setZero();
			derivative->block<1, 3>(0, 3).noalias() += 2.0 * (q1 - q0).transpose();
			derivative->block<1, 3>(0, 0).noalias() -= 2.0 * (q1 - q0).transpose();
			derivative->block<1, 3>(1, 6).noalias() += (q1 - q0).transpose();
			derivative->block<1, 3>(1, 3).noalias() += (q2 - q0).transpose();
			derivative->block<1, 3>(1, 0).noalias() += -(q1 - q0).transpose() - (q2 - q0).transpose();
			derivative->block<1, 3>(2, 6).noalias() += (q1 - q0).transpose();
			derivative->block<1, 3>(2, 3).noalias() += (q2 - q0).transpose();
			derivative->block<1, 3>(2, 0).noalias() += -(q1 - q0).transpose() - (q2 - q0).transpose();
			derivative->block<1, 3>(3, 6).noalias() += 2.0 * (q2 - q0).transpose();
			derivative->block<1, 3>(3, 0).noalias() -= 2.0 * (q2 - q0).transpose();
		}

		if (hessian)
		{
			for (int i = 0; i < 4; i++)
			{
				(*hessian)[i].setZero();
			}
			// std::cout << "hessian" << std::endl;
			(*hessian)[0].block<3, 3>(0, 0).diagonal().noalias() += 2.0 * Vector3::Ones();
			(*hessian)[0].block<3, 3>(3, 3).diagonal().noalias() += 2.0 * Vector3::Ones();
			(*hessian)[0].block<3, 3>(0, 3).diagonal().noalias() -= 2.0 * Vector3::Ones();
			(*hessian)[0].block<3, 3>(3, 0).diagonal().noalias() -= 2.0 * Vector3::Ones();

			(*hessian)[1].block<3, 3>(3, 6).diagonal().noalias() += Vector3::Ones();
			(*hessian)[1].block<3, 3>(6, 3).diagonal().noalias() += Vector3::Ones();
			(*hessian)[1].block<3, 3>(0, 3).diagonal().noalias() -= Vector3::Ones();
			(*hessian)[1].block<3, 3>(0, 6).diagonal().noalias() -= Vector3::Ones();
			(*hessian)[1].block<3, 3>(3, 0).diagonal().noalias() -= Vector3::Ones();
			(*hessian)[1].block<3, 3>(6, 0).diagonal().noalias() -= Vector3::Ones();
			(*hessian)[1].block<3, 3>(0, 0).diagonal().noalias() += 2.0 * Vector3::Ones();

			(*hessian)[2].block<3, 3>(3, 6).diagonal().noalias() += Vector3::Ones();
			(*hessian)[2].block<3, 3>(6, 3).diagonal().noalias() += Vector3::Ones();
			(*hessian)[2].block<3, 3>(0, 3).diagonal().noalias() -= Vector3::Ones();
			(*hessian)[2].block<3, 3>(0, 6).diagonal().noalias() -= Vector3::Ones();
			(*hessian)[2].block<3, 3>(3, 0).diagonal().noalias() -= Vector3::Ones();
			(*hessian)[2].block<3, 3>(6, 0).diagonal().noalias() -= Vector3::Ones();
			(*hessian)[2].block<3, 3>(0, 0).diagonal().noalias() += 2.0 * Vector3::Ones();

			(*hessian)[3].block<3, 3>(0, 0).diagonal().noalias() += 2.0 * Vector3::Ones();
			(*hessian)[3].block<3, 3>(6, 6).diagonal().noalias() += 2.0 * Vector3::Ones();
			(*hessian)[3].block<3, 3>(0, 6).diagonal().noalias() -= 2.0 * Vector3::Ones();
			(*hessian)[3].block<3, 3>(6, 0).diagonal().noalias() -= 2.0 * Vector3::Ones();
		}

		return result;
	}

	Matrix22 firstFundamentalFormNew(
	    const Matrix33& vertexPositions,
	    bool requireDerivative,
	    bool requireHessian,
	    Eigen::Matrix<Real, 4, 9>* derivative, // F(face, i)
	    std::array<Eigen::Matrix<Real, 9, 9>, 4>* hessian)
	{
		Vector3 q0 = vertexPositions.col(0);
		Vector3 q1 = vertexPositions.col(1);
		Vector3 q2 = vertexPositions.col(2);
		Matrix22 result;
		result << (q1 - q0).dot(q1 - q0), (q1 - q0).dot(q2 - q0),
		    (q2 - q0).dot(q1 - q0), (q2 - q0).dot(q2 - q0);

		if (requireDerivative)
		{
			// std::cout<< "derivative" << std::endl;
			derivative->setZero();
			derivative->block<1, 3>(0, 3).noalias() += 2.0 * (q1 - q0).transpose();
			derivative->block<1, 3>(0, 0).noalias() -= 2.0 * (q1 - q0).transpose();
			derivative->block<1, 3>(1, 6).noalias() += (q1 - q0).transpose();
			derivative->block<1, 3>(1, 3).noalias() += (q2 - q0).transpose();
			derivative->block<1, 3>(1, 0).noalias() += -(q1 - q0).transpose() - (q2 - q0).transpose();
			derivative->block<1, 3>(2, 6).noalias() += (q1 - q0).transpose();
			derivative->block<1, 3>(2, 3).noalias() += (q2 - q0).transpose();
			derivative->block<1, 3>(2, 0).noalias() += -(q1 - q0).transpose() - (q2 - q0).transpose();
			derivative->block<1, 3>(3, 6).noalias() += 2.0 * (q2 - q0).transpose();
			derivative->block<1, 3>(3, 0).noalias() -= 2.0 * (q2 - q0).transpose();
		}

		if (requireHessian)
		{
			for (int i = 0; i < 4; i++)
			{
				(*hessian)[i].setZero();
			}
			// std::cout << "hessian" << std::endl;
			(*hessian)[0].block<3, 3>(0, 0).diagonal().noalias() += 2.0 * Vector3::Ones();
			(*hessian)[0].block<3, 3>(3, 3).diagonal().noalias() += 2.0 * Vector3::Ones();
			(*hessian)[0].block<3, 3>(0, 3).diagonal().noalias() -= 2.0 * Vector3::Ones();
			(*hessian)[0].block<3, 3>(3, 0).diagonal().noalias() -= 2.0 * Vector3::Ones();

			(*hessian)[1].block<3, 3>(3, 6).diagonal().noalias() += Vector3::Ones();
			(*hessian)[1].block<3, 3>(6, 3).diagonal().noalias() += Vector3::Ones();
			(*hessian)[1].block<3, 3>(0, 3).diagonal().noalias() -= Vector3::Ones();
			(*hessian)[1].block<3, 3>(0, 6).diagonal().noalias() -= Vector3::Ones();
			(*hessian)[1].block<3, 3>(3, 0).diagonal().noalias() -= Vector3::Ones();
			(*hessian)[1].block<3, 3>(6, 0).diagonal().noalias() -= Vector3::Ones();
			(*hessian)[1].block<3, 3>(0, 0).diagonal().noalias() += 2.0 * Vector3::Ones();

			(*hessian)[2].block<3, 3>(3, 6).diagonal().noalias() += Vector3::Ones();
			(*hessian)[2].block<3, 3>(6, 3).diagonal().noalias() += Vector3::Ones();
			(*hessian)[2].block<3, 3>(0, 3).diagonal().noalias() -= Vector3::Ones();
			(*hessian)[2].block<3, 3>(0, 6).diagonal().noalias() -= Vector3::Ones();
			(*hessian)[2].block<3, 3>(3, 0).diagonal().noalias() -= Vector3::Ones();
			(*hessian)[2].block<3, 3>(6, 0).diagonal().noalias() -= Vector3::Ones();
			(*hessian)[2].block<3, 3>(0, 0).diagonal().noalias() += 2.0 * Vector3::Ones();

			(*hessian)[3].block<3, 3>(0, 0).diagonal().noalias() += 2.0 * Vector3::Ones();
			(*hessian)[3].block<3, 3>(6, 6).diagonal().noalias() += 2.0 * Vector3::Ones();
			(*hessian)[3].block<3, 3>(0, 6).diagonal().noalias() -= 2.0 * Vector3::Ones();
			(*hessian)[3].block<3, 3>(6, 0).diagonal().noalias() -= 2.0 * Vector3::Ones();

			// Matrix33 I = Matrix33::Identity();
			// (*hessian)[0].block<3, 3>(0, 0).noalias() += 2.0 * I;
			// (*hessian)[0].block<3, 3>(3, 3).noalias() += 2.0 * I;
			// (*hessian)[0].block<3, 3>(0, 3) -= 2.0 * I;
			// (*hessian)[0].block<3, 3>(3, 0) -= 2.0 * I;
			//
			// (*hessian)[1].block<3, 3>(3, 6) +=Vector3::Ones();
			// (*hessian)[1].block<3, 3>(6, 3) +=Vector3::Ones();
			// (*hessian)[1].block<3, 3>(0, 3) -=Vector3::Ones();
			// (*hessian)[1].block<3, 3>(0, 6) -=Vector3::Ones();
			// (*hessian)[1].block<3, 3>(3, 0) -=Vector3::Ones();
			// (*hessian)[1].block<3, 3>(6, 0) -=Vector3::Ones();
			// (*hessian)[1].block<3, 3>(0, 0) += 2.0 * I;
			//
			// (*hessian)[2].block<3, 3>(3, 6) +=Vector3::Ones();
			// (*hessian)[2].block<3, 3>(6, 3) +=Vector3::Ones();
			// (*hessian)[2].block<3, 3>(0, 3) -=Vector3::Ones();
			// (*hessian)[2].block<3, 3>(0, 6) -=Vector3::Ones();
			// (*hessian)[2].block<3, 3>(3, 0) -=Vector3::Ones();
			// (*hessian)[2].block<3, 3>(6, 0) -=Vector3::Ones();
			// (*hessian)[2].block<3, 3>(0, 0) += 2.0 * I;
			//
			// (*hessian)[3].block<3, 3>(0, 0) += 2.0 * I;
			// (*hessian)[3].block<3, 3>(6, 6) += 2.0 * I;
			// (*hessian)[3].block<3, 3>(0, 6) -= 2.0 * I;
			// (*hessian)[3].block<3, 3>(6, 0) -= 2.0 * I;
		}

		return result;
	}

	Eigen::Matrix2d firstFundamentalForm(
	    const MeshConnectivity& mesh,
	    const Eigen::MatrixXd& curPos,
	    int face,
	    Eigen::Matrix<double, 4, 9>* derivative, // F(face, i)
	    std::vector<Eigen::Matrix<double, 9, 9>>* hessian)
	{
		Eigen::Vector3d q0 = curPos.row(mesh.faceVertex(face, 0));
		Eigen::Vector3d q1 = curPos.row(mesh.faceVertex(face, 1));
		Eigen::Vector3d q2 = curPos.row(mesh.faceVertex(face, 2));
		Eigen::Matrix2d result;
		result << (q1 - q0).dot(q1 - q0), (q1 - q0).dot(q2 - q0),
		    (q2 - q0).dot(q1 - q0), (q2 - q0).dot(q2 - q0);

		if (derivative)
		{
			derivative->setZero();
			derivative->block<1, 3>(0, 3) += 2.0 * (q1 - q0).transpose();
			derivative->block<1, 3>(0, 0) -= 2.0 * (q1 - q0).transpose();
			derivative->block<1, 3>(1, 6) += (q1 - q0).transpose();
			derivative->block<1, 3>(1, 3) += (q2 - q0).transpose();
			derivative->block<1, 3>(1, 0) += -(q1 - q0).transpose() - (q2 - q0).transpose();
			derivative->block<1, 3>(2, 6) += (q1 - q0).transpose();
			derivative->block<1, 3>(2, 3) += (q2 - q0).transpose();
			derivative->block<1, 3>(2, 0) += -(q1 - q0).transpose() - (q2 - q0).transpose();
			derivative->block<1, 3>(3, 6) += 2.0 * (q2 - q0).transpose();
			derivative->block<1, 3>(3, 0) -= 2.0 * (q2 - q0).transpose();
		}

		if (hessian)
		{
			hessian->resize(4);
			for (int i = 0; i < 4; i++)
			{
				(*hessian)[i].setZero();
			}
			Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
			(*hessian)[0].block<3, 3>(0, 0) += 2.0 * I;
			(*hessian)[0].block<3, 3>(3, 3) += 2.0 * I;
			(*hessian)[0].block<3, 3>(0, 3) -= 2.0 * I;
			(*hessian)[0].block<3, 3>(3, 0) -= 2.0 * I;

			(*hessian)[1].block<3, 3>(3, 6) += I;
			(*hessian)[1].block<3, 3>(6, 3) += I;
			(*hessian)[1].block<3, 3>(0, 3) -= I;
			(*hessian)[1].block<3, 3>(0, 6) -= I;
			(*hessian)[1].block<3, 3>(3, 0) -= I;
			(*hessian)[1].block<3, 3>(6, 0) -= I;
			(*hessian)[1].block<3, 3>(0, 0) += 2.0 * I;

			(*hessian)[2].block<3, 3>(3, 6) += I;
			(*hessian)[2].block<3, 3>(6, 3) += I;
			(*hessian)[2].block<3, 3>(0, 3) -= I;
			(*hessian)[2].block<3, 3>(0, 6) -= I;
			(*hessian)[2].block<3, 3>(3, 0) -= I;
			(*hessian)[2].block<3, 3>(6, 0) -= I;
			(*hessian)[2].block<3, 3>(0, 0) += 2.0 * I;

			(*hessian)[3].block<3, 3>(0, 0) += 2.0 * I;
			(*hessian)[3].block<3, 3>(6, 6) += 2.0 * I;
			(*hessian)[3].block<3, 3>(0, 6) -= 2.0 * I;
			(*hessian)[3].block<3, 3>(6, 0) -= 2.0 * I;
		}

		return result;
	}
}
