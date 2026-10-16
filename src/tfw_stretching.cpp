#include <argus/tfw_stretching.hpp>
#include <argus/util.hpp>
#include <argus/geometry_functions.hpp>
#include "../../vendor/SecondFundamentalForm/SecondFundamentalFormDiscretization.h"
#include "../../vendor/SecondFundamentalForm/MidedgeAverageFormulation.h"
#include "halfedge/trimesh.h"
#include "halfedge/trimesh_types.h"
#define idxAt(i, j) 3 * i + j
#define mat(A, i, j) A(0, idxAt(i, j)), A(1, idxAt(i, j)), A(2, idxAt(i, j)), A(3, idxAt(i, j))
#define hessMat(A, i, j) A[0](i, j), A[1](i, j), A[2](i, j), A[3](i, j)
namespace argus
{
#ifdef ARGUS_CHECKPOINT
	void TFWStretch::saveInfo()
	{
		nlohmann::json data;
		data["Energy"] = "TFW Stretching";
		data["Is external"] = isExternal();
		simConf["Energies"].push_back(data);
	}
#endif
	// function for per vertex derivatives and hessians amps
	Real TFWStretch::getEnergyDensityFromQuadPerVtex1(int vid, int faceId, int quadId, Sheet& sheet, VectorN* deriv, MatrixNN* hess) const
	{
		Real energy = 0;
		const std::vector<trimesh::triangle_t>& faces = sheet.getTriangles();
		const std::vector<trimesh::edge_t>& edges = sheet.getEdges();
		const Matrix3N& positions = sheet.getPositionsN3();
		const Eigen::MatrixXi& edgeOppositeVertices = sheet.getEdgeOppositeVertices();
		const std::vector<trimesh::triangle_t>& vertexOppositeEdges = sheet.getVertexOppositeEdgesPerFace();
		int nfaces = faces.size();
		int nedges = edges.size();
		int nvertices = positions.cols();
		Matrix33 facePositions;
		int vidxi = faces[faceId].v[0];
		int vidxj = faces[faceId].v[1];
		int vidxk = faces[faceId].v[2];
		// facePositions.col(0) = positions.block<3, 1>(0, faces[faceId].v[0]);
		// facePositions.col(1) = positions.block<3, 1>(0, faces[faceId].v[1]);
		// facePositions.col(2) = positions.block<3, 1>(0, faces[faceId].v[2]);

		facePositions.col(0) = positions.col(faces[faceId].v[0]);
		facePositions.col(1) = positions.col(faces[faceId].v[1]);
		facePositions.col(2) = positions.col(faces[faceId].v[2]);
		MatrixNN gradI, gradII; // gradI(i) contains dI(i)/dX
		std::vector<Eigen::Matrix<Real, 9, 9>> hessI;
		Matrix22 I;
		I = firstFundamentalForm(facePositions, true, true, &gradI, &hessI); // correct till here
		std::vector<Eigen::MatrixXd> hessII;
		std::shared_ptr<SecondFundamentalFormDiscretization> sff;
		sff = std::make_shared<MidedgeAverageFormulation>();
		Matrix22 II = sff->secondFundamentalForm(sheet.getBaseMesh(), sheet.getPositionsN3().transpose(), sheet.getEdgeDofs(), faceId, &gradII, &hessII);
		Matrix22 Ibar = sheet.getRestShape().getIbars()[faceId];
		Matrix22 Ibarinv = Ibar.inverse();
		Matrix22 IIbar = sheet.getRestShape().getIIbars()[faceId];
		Vector3 gradA; // datilda/da
		Vector2 da;
		Matrix23 gradDA; // dda/da
		Real u = sheet.getQuadPoints()[quadId].u;
		Real v = sheet.getQuadPoints()[quadId].v;
		Vector3 ampUVWeights;
		ampUVWeights << 1 - u - v, u, v;
		Vector2 dphi = sheet.computeDphi1VectorPerFace(faceId, NULL);
		Real lameAlpha = sheet.getLameAlpha();
		Real lameBeta = sheet.getLameBeta();
		Vector3 amps;
		for (int i = 0; i < 3; i++)
		{
			// vertices(i) = _state.baseMesh.faceVertex(faceId, i);
			int vid = sheet.getTriangles()[faceId].v[i];
			amps(i) = sheet.getAmplitudes()(vid);
		}
		auto func = [vid, u, v, I, II, Ibar, IIbar, &sheet](Eigen::Vector<Real, 5> input, VectorN* d, MatrixNN* h)
		{
			Real lA, lB;
			lA = sheet.getLameAlpha();
			lB = sheet.getLameBeta();
			AReal3 energy = 0;
			Vector3 finput;
			finput(0) = input(vid);
			finput(1) = input(3);
			finput(2) = input(4);
			Eigen::Vector<AReal3, 3> X;

			X = AReal3::make_active(finput);
			Eigen::Vector<AReal3, 2> da;
			da(0) = X(1) - X(0);
			da(1) = X(2) - X(0);
			AReal3 a;
			if (vid == 0)
				a = u * input(1) + v * input(2) + (1 - u - v) * X(0);
			else if (vid == 1)
				a = u * X(0) + v * input(2) + (1 - u - v) * input(0);
			else if (vid == 2)
				a = u * input(1) + v * X(0) + (1 - u - v) * input(0);
			else
				throw "error in tfw_stretching";
			// dphis
			Eigen::Vector<AReal3, 2> dphi = X(Eigen::seq(1, 2));
			Eigen::Matrix<AReal3, 2, 2> M1, M3, M4;
			// first term
			M1 = Ibar.inverse() * (I - Ibar + 0.5 * da * da.transpose() + 0.5 * a * a * dphi * dphi.transpose());
			AReal3 e1 = 0.5 * lA * (M1.trace() * M1.trace()) + lB * (M1 * M1).trace();
			energy += e1;

			// second term
			if (dphi.norm() > 1e-4)
			{
				Eigen::Vector<AReal3, 2> w = Ibar.inverse() * dphi;
				AReal3 num = w.norm() * w.norm() * II.trace() - w.transpose() * II * w;
				AReal3 den = w.norm() * w.norm() * Ibar.trace() - w.transpose() * Ibar * w;
				AReal3 e2 = 4.0 * a * a * lB * ((lA + lB) / (lA + 2.0 * lB)) * num * num / (den * den);
				energy += e2;
			}

			// third term
			M3 = Ibar.inverse() * (dphi * da.transpose() + da * dphi.transpose());
			AReal3 e3 = a * a * (1 / 32.0) * (0.5 * lA * (M3.trace() * M3.trace()) + lB * (M3 * M3).trace());
			energy += e3;

			// forth term
			M4 = Ibar.inverse() * da * da.transpose();
			AReal3 e4 = (1 / 8.0) * (0.5 * lA * (M4.trace() * M4.trace()) + lB * (M4 * M4).trace());
			energy += e4;
			if (d)
			{
				(*d) = energy.grad;
			}
			if (h)
				(*h) = energy.Hess;
			return TinyAD::to_passive(energy);
		};
		if (deriv)
		{
			deriv->resize(3);
			deriv->setZero();
		}
		if (hess)
		{
			hess->resize(3, 3);
			hess->setZero();
		}
		Eigen::Vector<Real, 5> input;
		input(Eigen::seq(0, 2)) = amps;
		input(Eigen::seq(3, 4)) = sheet.getDphis1VectorPerFace().col(faceId);
		energy = func(input, deriv, hess);
		Real area
		    = 0.5 * sqrt(Ibar.determinant());
		if (deriv)
			*deriv = area * (*deriv);
		if (hess)
		{
			(*hess) = area * (*hess);
			// (*hess) = util::lowRankApprox(*hess);
		}

		return area * energy;
	}

	Real TFWStretch::getVertexEnergy1(Sheet& sheet, int vid, std::vector<int>& fids, VectorN* deriv, MatrixNN* hess) const
	{
		// need to compute deriv and hessian wrt the vertex at index vid
		// just for testing
		Real energy = 0;
		const Eigen::VectorXi& faces = sheet.sheetInfo.vertexFacesSet.row(vid);
		unsigned int nfaces = 0;
		for (int i = 0; i < faces.rows(); i++)
		{
			if (faces(i) == -1)
				break;
			nfaces++;
		}
		// create buffers to store parallel results
		std::vector<Real> energies(nfaces, 0);
		std::vector<VectorN> derivatives(nfaces);

		if (deriv)
		{
			deriv->resize(nfaces * 2 + 1);
			deriv->setZero();
		}
		if (hess)
		{
			hess->resize(nfaces * 2 + 1, nfaces * 2 + 1);
			hess->setZero();
		}

		std::vector<MatrixNN> hessians(nfaces);
		for (unsigned int i = 0; i < nfaces; i++)
		{
			derivatives[i].resize(2 * nfaces + 1);
			derivatives[i].setZero();
			hessians[i].resize(2 * nfaces + 1, 2 * nfaces + 1);
			hessians[i].setZero();
		}

		// #ifdef USE_OMP
		// #pragma omp parallel for
		// #endif
		for (int i = 0; i < nfaces; i++)
		{
			int faceId = faces(i);
			int localVid = -1;
			for (int ii = 0; ii < 3; ii++)
			{
				int vidii = sheet.getTriangles().at(faceId).v[ii];
				if (vidii == vid)
				{
					localVid = ii;
					break;
				}
			}
			if (localVid == -1)
			{
				std::cout << "Vertex " << localVid << " not found! Face " << faceId << ". Error in tfw_stretching.cpp\n";
				exit(0);
			}
			for (int quad = 0; quad < sheet.getQuadPoints().size(); quad++)
			{
				VectorN tempDeriv;
				MatrixNN tempHess;
				energies[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * TFWStretch::getEnergyDensityFromQuadPerVtex1(localVid, faceId, quad, sheet, deriv ? &tempDeriv : NULL, hess ? &tempHess : NULL) / 4.0;
				if (deriv)
				{
					tempDeriv = (sheet.getQuadPoints()[quad].weight * sheet.getThickness() * tempDeriv / 4.0).eval();
					derivatives[i](0) += tempDeriv(0);
					derivatives[i](1 + 2 * i) += tempDeriv(1);
					derivatives[i](1 + 2 * i + 1) += tempDeriv(2);
				}
				if (hess)
				{
					tempHess = (sheet.getQuadPoints()[quad].weight * sheet.getThickness() * tempHess / 4.0).eval();
					hessians[i](0, 0) += tempHess(0, 0);
					hessians[i](0, 2 * i + 1) += tempHess(0, 1);
					hessians[i](0, 2 * i + 2) += tempHess(0, 2);
					hessians[i](2 * i + 1, 0) += tempHess(1, 0);
					hessians[i](2 * i + 2, 0) += tempHess(2, 0);
					hessians[i](1 + 2 * i, 1 + 2 * i) += tempHess(1, 1);
					hessians[i](1 + 2 * i, 1 + 2 * i + 1) += tempHess(1, 2);
					hessians[i](1 + 2 * i + 1, 1 + 2 * i) += tempHess(2, 1);
					hessians[i](1 + 2 * i + 1, 1 + 2 * i + 1) += tempHess(2, 2);
				}
			}
		}

		// consolidate
		for (int i = 0; i < nfaces; i++)
		{
			if (std::find(fids.begin(), fids.end(), faces(i)) == fids.end())
				fids.push_back(faces(i));
			energy += energies[i];
			if (deriv)
				(*deriv) += derivatives[i];
			if (hess)
				(*hess) += hessians[i];
		}
		return energy;
	}
	// function for per vertex derivatives and hessians amps
	Real TFWStretch::getEnergyDensityFromQuadPerVtex1(int vid, int faceId, int quadId, Sheet& sheet, Real* deriv, Real* hess) const
	{
		Real energy = 0;
		const std::vector<trimesh::triangle_t>& faces = sheet.getTriangles();
		const std::vector<trimesh::edge_t>& edges = sheet.getEdges();
		const Matrix3N& positions = sheet.getPositionsN3();
		const Eigen::MatrixXi& edgeOppositeVertices = sheet.getEdgeOppositeVertices();
		const std::vector<trimesh::triangle_t>& vertexOppositeEdges = sheet.getVertexOppositeEdgesPerFace();
		int nfaces = faces.size();
		int nedges = edges.size();
		int nvertices = positions.cols();

		Matrix33 facePositions;
		int vidxi = faces[faceId].v[0];
		int vidxj = faces[faceId].v[1];
		int vidxk = faces[faceId].v[2];
		// facePositions.col(0) = positions.block<3, 1>(0, faces[faceId].v[0]);
		// facePositions.col(1) = positions.block<3, 1>(0, faces[faceId].v[1]);
		// facePositions.col(2) = positions.block<3, 1>(0, faces[faceId].v[2]);

		facePositions.col(0) = positions.col(faces[faceId].v[0]);
		facePositions.col(1) = positions.col(faces[faceId].v[1]);
		facePositions.col(2) = positions.col(faces[faceId].v[2]);
		Matrix22 I = firstFundamentalFormNew(facePositions, false, false, NULL, NULL); // correct till here
		Matrix22 Ibar = sheet.getRestShape().getIbars()[faceId];
		Matrix22 Ibarinv = Ibar.inverse();

		Vector3 gradA; // datilda/da
		Vector2 da;
		// Matrix23 gradDA; // dda/da
		Real u = sheet.getQuadPoints()[quadId].u;
		Real v = sheet.getQuadPoints()[quadId].v;
		Vector3 ampUVWeights;
		ampUVWeights << 1 - u - v, u, v;
		Real a = sheet.computeAmplitudesFromQuad(faceId, quadId, &da, &gradA, NULL, nullptr, nullptr);
		// Matrix22 gradDphi;
		Vector2 dphi = sheet.computeDphi1VectorPerFace(faceId, NULL);
		// std::array<Matrix22,2> gradDphiDphiTensor;
		// std::vector<Matrix22> hessDphiDphiT;
		Matrix22 dphidphiTensor;
		dphidphiTensor = sheet.computeDphiDphi1VectorPerFaceTensor(faceId, quadId, NULL, NULL);
		std::array<Matrix22, 3> gradDaDaTensor;
		std::array<Matrix33, 4> hessDaDaTensor;
		Matrix22 dadaTensor;
		dadaTensor = sheet.computeDaDaTensorNew(faceId, quadId, &gradDaDaTensor, &hessDaDaTensor);
		std::array<Matrix22, 5> gradDaDphiTensor;
		// std::vector<Eigen::Matrix<Real, 5, 5>> hessDaDphiTensor, hessDphiDaTensor;
		Matrix22 dadphiTensor = sheet.computeDaDphi1VectorPerFaceTensorNew(faceId, quadId, &gradDaDphiTensor, NULL);
		std::array<Matrix22, 5> gradDphiDaTensor;
		Matrix22 dphidaTensor = sheet.computeDphi1VectorPerFaceDaTensorNew(faceId, quadId, &gradDphiDaTensor, NULL);
		Real lameAlpha = sheet.getLameAlpha();
		Real lameBeta = sheet.getLameBeta();

		if (deriv)
			*deriv = 0;
		if (hess)
			*hess = 0;
		// First term
		// Correct numerically
		Matrix22 M0 = Ibarinv * (I + 0.5 * a * a * dphidphiTensor + 0.5 * dadaTensor - Ibar);
		Matrix22 dMSVdM;
		Real MSV;
		MSV = SVNormSquared(lameAlpha, lameBeta, M0, &dMSVdM);
		energy += MSV;
		Matrix22 dMda;
		if (deriv || hess)
			dMda = Ibarinv * (a * gradA(vid) * dphidphiTensor + 0.5 * gradDaDaTensor[vid]);

		if (deriv)
		{
			// wrt amp
			(*deriv) += doubleContraction(dMSVdM, dMda);
		}

		if (hess)
		{
			// amps wrt amps
			Matrix22 hessDaDa;
			hessDaDa << hessMat(hessDaDaTensor, vid, vid);
			Matrix22 d2Mdada = Ibarinv * (gradA(vid) * gradA(vid) * dphidphiTensor + 0.5 * hessDaDa);
			// Matrix22 dMda = Ibarinv * (a * gradA(vid) * dphidphiTensor + 0.5 * gradDaDaTensor.at(vid));
			Matrix22 d2MSVdadM = lameAlpha * (dMda).trace() * Matrix22::Identity() + 2 * lameBeta * dMda;
			(*hess) += doubleContraction(d2MSVdadM, dMda) + doubleContraction(dMSVdM, d2Mdada);
		}
		// Second term
		if (dphi.norm() > 1e-4)
		{
			std::shared_ptr<SecondFundamentalFormDiscretization> sff;
			sff = std::make_shared<MidedgeAverageFormulation>();
			Matrix22 II = sff->secondFundamentalFormNew(sheet.getBaseMesh(), sheet.getPositionsN3(), sheet.getEdgeDofs(), faceId, NULL, NULL);
			Matrix22 IIbar = sheet.getRestShape().getIIbars()[faceId];

			Vector2 w = Ibarinv * dphi;
			Real num = w.norm() * w.norm() * II.trace() - w.transpose() * II * w;
			Real den = w.norm() * w.norm() * Ibar.trace() - w.transpose() * Ibar * w;
			Real c = 4 * lameBeta * ((lameAlpha + lameBeta) / (lameAlpha + 2 * lameBeta)) / (den * den);
			energy += a * a * c * num * num;
			Matrix22 dnumdII;
			if (deriv || hess)
				dnumdII = w.norm() * w.norm() * Matrix22::Identity() - w * w.transpose();
			if (deriv)
			{
				// wrt amps
				(*deriv) += 2 * a * gradA(vid) * c * num * num;
			}
			if (hess)
			{
				// amps vs amps
				(*hess) += 2 * gradA(vid) * gradA(vid) * c * num * num;
			}
		}

		// Third term
		// Numerically correct
		Matrix22 M2 = (Ibarinv * (dphidaTensor + dadphiTensor));
		Matrix22 dM2SVdM;
		Real M2SVM = SVNormSquared(sheet.getLameAlpha(), sheet.getLameBeta(), M2, &dM2SVdM);
		energy += (a * a / 32.0) * M2SVM;
		if (deriv)
		{
			// wrt amp
			Matrix22 DdphidaDa, DdadphiDa;
			DdphidaDa = gradDphiDaTensor[vid];
			DdadphiDa = gradDaDphiTensor[vid];
			// Matrix22 dMda = (gradA(vid) / std::sqrt(32.0)) * (Ibarinv * (dphidaTensor + dadphiTensor)) + (a / std::sqrt(32.0)) * (Matrix22 );
			// (*deriv)(0) += doubleContraction(dM2SVdM, dMda);
			(*deriv) += (a * gradA(vid) / 16.0) * M2SVM + (a * a / 32.0) * doubleContraction(dM2SVdM, Ibarinv * (DdphidaDa + DdadphiDa));
		}
		if (hess)
		{
			Matrix22 DdphidaDa, DdadphiDa;
			DdphidaDa = gradDphiDaTensor[vid];
			DdadphiDa = gradDaDphiTensor[vid];
			Matrix22 dMda, d2MSVdadM;
			dMda = Ibarinv * (gradDphiDaTensor[vid] + gradDaDphiTensor[vid]);
			d2MSVdadM = lameAlpha * dMda.trace() * Matrix22::Identity() + 2 * lameBeta * dMda;
			// wrt amps
			(*hess)
			    += (2 * gradA(vid) * gradA(vid) / 32.0) * M2SVM + 2 * (2 * gradA(vid) * a / 32.0) * doubleContraction(dM2SVdM, Ibarinv * (DdphidaDa + DdadphiDa))
			    + (a * a / 32.0) * doubleContraction(d2MSVdadM, dMda);
		}

		// Fourth term
		// Correct numerically
		Matrix22 M3 = Ibarinv * dadaTensor;
		Matrix22 dM3SVdM;
		energy += (1 / 8.0) * SVNormSquared(sheet.getLameAlpha(), sheet.getLameBeta(), M3, &dM3SVdM);
		if (deriv)
		{
			// wrt amps
			Matrix22 dMda = Ibarinv * gradDaDaTensor[vid];
			(*deriv) += (1 / 8.0) * doubleContraction(dM3SVdM, dMda);
		}
		if (hess)
		{
			// wrt amps
			Matrix22 dMda = Ibarinv * gradDaDaTensor[vid];
			Matrix22 d2MSVdadM = lameAlpha * dMda.trace() * Matrix22::Identity() + 2 * lameBeta * dMda;
			Matrix22 hessDaDa;
			hessDaDa << hessMat(hessDaDaTensor, vid, vid);
			(*hess) += (1 / 8.0) * (doubleContraction(d2MSVdadM, dMda) + doubleContraction(dM3SVdM, Ibarinv * hessDaDa));
		}

		Real area
		    = 0.5 * sqrt(Ibar.determinant());
		if (deriv)
			*deriv = area * (*deriv);
		if (hess)
			(*hess) = area * (*hess);

		return area * energy;
	}

	Real TFWStretch::getVertexEnergy1New(Sheet& sheet, int vid, Real* deriv, Real* hess) const
	{

		// need to compute deriv and hessian wrt the vertex at index vid
		// just for testing
		Real energy = 0;
		const Eigen::Vector<int, 32>& faces = sheet.sheetInfo.vertexFacesSet.row(vid);
		unsigned int nfaces = 0;
		for (int i = 0; i < faces.rows(); i++)
		{
			if (faces(i) == -1)
				break;
			nfaces++;
		}
		// create buffers to store parallel results
		if (deriv)
			*deriv = 0;
		if (hess)
			*hess = 0;

		for (int i = 0; i < nfaces; i++)
		{
			int faceId = faces(i);
			int localVid = -1;
			for (int ii = 0; ii < 3; ii++)
			{
				int vidii = sheet.getTriangles().at(faceId).v[ii];
				if (vidii == vid)
				{
					localVid = ii;
					break;
				}
			}
			if (localVid == -1)
			{
				std::cout << "Vertex " << localVid << " not found! Face " << faceId << ". Error in tfw_stretching.cpp\n";
				exit(0);
			}
			for (int quad = 0; quad < sheet.getQuadPoints().size(); quad++)
			{
				Real tempDeriv, tempHess;
				energy += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * TFWStretch::getEnergyDensityFromQuadPerVtex1(localVid, faceId, quad, sheet, deriv ? &tempDeriv : NULL, hess ? &tempHess : NULL) / 4.0;
				if (deriv)
					(*deriv) += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * tempDeriv / 4.0;
				if (hess)
					(*hess) += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * tempHess / 4.0;
			}
		}
		// if (*hess < 0)
		// (*hess) = -1 * (*hess);
		// consolidate
		return energy;
	}
	Real TFWStretch::getVertexEnergy1(Sheet& sheet, int vid, Real* deriv, Real* hess) const
	{
		// need to compute deriv and hessian wrt the vertex at index vid
		// just for testing
		Real totalEnergy = 0;
		const Eigen::Vector<int, 32>& facesVec = sheet.sheetInfo.vertexFacesSet.row(vid);
		unsigned int nAdjFaces = 0;
		for (int i = 0; i < facesVec.rows(); i++)
		{
			if (facesVec(i) == -1)
				break;
			nAdjFaces++;
		}
		// create buffers to store parallel results
		if (deriv)
			*deriv = 0;
		if (hess)
			*hess = 0;

		for (int i = 0; i < nAdjFaces; i++)
		{
			int faceId = facesVec(i);
			int localVid = -1;
			for (int ii = 0; ii < 3; ii++)
			{
				int vidii = sheet.getTriangles().at(faceId).v[ii];
				if (vidii == vid)
				{
					localVid = ii;
					break;
				}
			}
			if (localVid == -1)
			{
				std::cout << "Vertex " << localVid << " not found! Face " << faceId << ". Error in tfw_stretching.cpp\n";
				exit(0);
			}
			const std::vector<trimesh::triangle_t>& faces = sheet.getTriangles();
			const std::vector<trimesh::edge_t>& edges = sheet.getEdges();
			const Matrix3N& positions = sheet.getPositionsN3();
			const Eigen::MatrixXi& edgeOppositeVertices = sheet.getEdgeOppositeVertices();
			const std::vector<trimesh::triangle_t>& vertexOppositeEdges = sheet.getVertexOppositeEdgesPerFace();
			int nfaces = faces.size();
			int nedges = edges.size();
			int nvertices = positions.cols();
			Matrix22 I;

			if (sheet.cachedIsValid())
			{
				I.noalias() = sheet.getCachedIs().second[faceId];
			}
			else
			{
				Matrix33 facePositions;
				int vidxi = faces[faceId].v[0];
				int vidxj = faces[faceId].v[1];
				int vidxk = faces[faceId].v[2];
				// facePositions.col(0) = positions.block<3, 1>(0, faces[faceId].v[0]);
				// facePositions.col(1) = positions.block<3, 1>(0, faces[faceId].v[1]);
				// facePositions.col(2) = positions.block<3, 1>(0, faces[faceId].v[2]);
				facePositions.col(0) = positions.col(faces[faceId].v[0]);
				facePositions.col(1) = positions.col(faces[faceId].v[1]);
				facePositions.col(2) = positions.col(faces[faceId].v[2]);
				I.noalias() = firstFundamentalFormNew(facePositions, false, false, NULL, NULL); // correct till here
			}
			Matrix22 Ibar = sheet.getRestShape().getIbars()[faceId];
			Matrix22 Ibarinv = Ibar.inverse();
			Matrix22 dphidphiTensor;
			if (sheet.cachedDphiDphiTValid())
			{
				dphidphiTensor = sheet.getCachedDphiDpiT().second[faceId];
			}
			else
			{
				dphidphiTensor = sheet.computeDphiDphi1VectorPerFaceTensor(faceId, 0, NULL, NULL);
			}

			Vector3 gradA; // datilda/da
			Vector2 da = sheet.computeDaFromFace(faceId);
			// Matrix23 gradDA; // dda/da

			// Matrix22 gradDphi;
			Vector2 dphi = sheet.computeDphi1VectorPerFace(faceId, NULL);
			// std::array<Matrix22,2> gradDphiDphiTensor;
			// std::vector<Matrix22> hessDphiDphiT;

			std::array<Matrix22, 3> gradDaDaTensor;
			std::array<Matrix33, 4> hessDaDaTensor;
			Matrix22 dadaTensor;
			std::array<Matrix22, 5> gradDaDphiTensor;
			std::array<Matrix22, 5> gradDphiDaTensor;
			Real lameAlpha = sheet.getLameAlpha();
			Real lameBeta = sheet.getLameBeta();
			dadaTensor = sheet.computeDaDaTensorNew(faceId, 0, &gradDaDaTensor, &hessDaDaTensor);

			Matrix22 dadphiTensor = sheet.computeDaDphi1VectorPerFaceTensorNew(faceId, 0, &gradDaDphiTensor, NULL);
			Matrix22 dphidaTensor = sheet.computeDphi1VectorPerFaceDaTensorNew(faceId, 0, &gradDphiDaTensor, NULL);
			Matrix22 II;
			Matrix22 IIbar = sheet.getRestShape().getIIbars()[faceId];

			// Second term
			Vector2 w;
			Real num, den, c2;
			if (dphi.norm() > 0)
			{
				if (sheet.cachedIIsValid())
				{
					II.noalias() = sheet.getCachedIIs().second[faceId];
				}
				else
				{
					// std::shared_ptr<SecondFundamentalFormDiscretization> sff;
					// sff = std::make_shared<MidedgeAverageFormulation>();
					exit(0);
					II.noalias() = sffobject->secondFundamentalFormNew(sheet.getBaseMesh(), sheet.getPositionsN3(), sheet.getEdgeDofs(), faceId, NULL, NULL);
				}

				w.noalias() = Ibarinv * dphi;
				num = w.norm() * w.norm() * II.trace() - w.transpose() * II * w;
				den = w.norm() * w.norm() * Ibar.trace() - w.transpose() * Ibar * w;
				c2 = 4 * lameBeta * ((lameAlpha + lameBeta) / (lameAlpha + 2 * lameBeta)) / (den * den);
			}
			// Third term
			Matrix22 M2 = (Ibarinv * (dphidaTensor + dadphiTensor));
			Matrix22 dM2SVdM;
			Real M2SVM = SVNormSquared(sheet.getLameAlpha(), sheet.getLameBeta(), M2, &dM2SVdM);

			// Fourth term
			// Correct numerically
			Matrix22 M3 = Ibarinv * dadaTensor;
			Matrix22 dM3SVdM;
			Real de4, he4;
			Real e4 = (1 / 8.0) * SVNormSquared(sheet.getLameAlpha(), sheet.getLameBeta(), M3, &dM3SVdM);
			if (deriv)
			{
				// wrt amps
				Matrix22 dMda = Ibarinv * gradDaDaTensor[localVid];
				de4 = (1 / 8.0) * doubleContraction(dM3SVdM, dMda);
			}
			if (hess)
			{
				// wrt amps
				Matrix22 dMda = Ibarinv * gradDaDaTensor[localVid];
				Matrix22 d2MSVdadM = lameAlpha * dMda.trace() * Matrix22::Identity() + 2 * lameBeta * dMda;
				Matrix22 hessDaDa;
				hessDaDa << hessMat(hessDaDaTensor, localVid, localVid);
				he4 = (1 / 8.0) * (doubleContraction(d2MSVdadM, dMda) + doubleContraction(dM3SVdM, Ibarinv * hessDaDa));
			}
			for (int quadId = 0; quadId < sheet.getQuadPoints().size(); quadId++)
			{
				Real energy = 0;
				Real tempDeriv = 0;
				Real tempHess = 0;
				Real a = sheet.computeAmplitudesFromQuad(faceId, quadId, NULL, &gradA, NULL, nullptr, nullptr);

				// First term
				// Correct numerically
				Matrix22 M0 = Ibarinv * (I + 0.5 * a * a * dphidphiTensor + 0.5 * dadaTensor - Ibar);
				Matrix22 dMSVdM;
				Real MSV;
				MSV = SVNormSquared(lameAlpha, lameBeta, M0, &dMSVdM);
				energy += MSV;
				Matrix22 dMda;
				if (deriv || hess)
					dMda = Ibarinv * (a * gradA(localVid) * dphidphiTensor + 0.5 * gradDaDaTensor[localVid]);

				if (deriv)
				{
					// wrt amp
					tempDeriv += doubleContraction(dMSVdM, dMda);
				}

				if (hess)
				{
					// amps wrt amps
					Matrix22 hessDaDa;
					hessDaDa << hessMat(hessDaDaTensor, localVid, localVid);
					Matrix22 d2Mdada = Ibarinv * (gradA(localVid) * gradA(localVid) * dphidphiTensor + 0.5 * hessDaDa);
					// Matrix22 dMda = Ibarinv * (a * gradA(localVid) * dphidphiTensor + 0.5 * gradDaDaTensor.at(localVid));
					Matrix22 d2MSVdadM = lameAlpha * (dMda).trace() * Matrix22::Identity() + 2 * lameBeta * dMda;
					tempHess += doubleContraction(d2MSVdadM, dMda) + doubleContraction(dMSVdM, d2Mdada);
				}
				// Second term
				if (dphi.norm() > 1e-4)
				{

					energy += a * a * c2 * num * num;
					Matrix22 dnumdII;
					if (deriv || hess)
						dnumdII = w.norm() * w.norm() * Matrix22::Identity() - w * w.transpose();
					if (deriv)
					{
						// wrt amps
						tempDeriv += 2 * a * gradA(localVid) * c2 * num * num;
					}
					if (hess)
					{
						// amps vs amps
						tempHess += 2 * gradA(localVid) * gradA(localVid) * c2 * num * num;
					}
				}

				// Third term
				// Numerically correct

				energy += (a * a / 32.0) * M2SVM;
				if (deriv)
				{
					// wrt amp
					// Matrix22 DdphidaDa, DdadphiDa;
					// DdphidaDa = gradDphiDaTensor[localVid];
					// DdadphiDa = gradDaDphiTensor[localVid];
					// Matrix22 dMda = (gradA(localVid) / std::sqrt(32.0)) * (Ibarinv * (dphidaTensor + dadphiTensor)) + (a / std::sqrt(32.0)) * (Matrix22 );
					// (*deriv)(0) += doubleContraction(dM2SVdM, dMda);
					tempDeriv += (a * gradA(localVid) / 16.0) * M2SVM + (a * a / 32.0) * doubleContraction(dM2SVdM, Ibarinv * (gradDphiDaTensor[localVid] + gradDaDphiTensor[localVid]));
				}
				if (hess)
				{
					// Matrix22 DdphidaDa, DdadphiDa;
					// DdphidaDa = gradDphiDaTensor[localVid];
					// DdadphiDa = gradDaDphiTensor[localVid];
					Matrix22 dMda, d2MSVdadM;
					dMda = Ibarinv * (gradDphiDaTensor[localVid] + gradDaDphiTensor[localVid]);
					// d2MSVdadM = lameAlpha * dMda.trace() * Matrix22::Identity() + 2 * lameBeta * dMda;
					// wrt amps
					tempHess
					    += (2 * gradA(localVid) * gradA(localVid) / 32.0) * M2SVM + 2 * (2 * gradA(localVid) * a / 32.0) * doubleContraction(dM2SVdM, Ibarinv * (gradDphiDaTensor[localVid] + gradDaDphiTensor[localVid]))
					    + (a * a / 32.0) * doubleContraction(lameAlpha * dMda.trace() * Matrix22::Identity() + 2 * lameBeta * dMda, dMda);
				}
				// 4th term
				energy += e4;
				if (deriv)
					tempDeriv += de4;
				if (hess)
					tempHess += he4;
				Real area = 0.5 * sqrt(Ibar.determinant());
				if (deriv)
					tempDeriv = area * tempDeriv;
				if (hess)
					tempHess = area * tempHess;

				totalEnergy += sheet.getQuadPoints()[quadId].weight * sheet.getThickness() * energy * area / 4.0;
				if (deriv)
					(*deriv) += sheet.getQuadPoints()[quadId].weight * sheet.getThickness() * tempDeriv / 4.0;
				if (hess)
					(*hess) += sheet.getQuadPoints()[quadId].weight * sheet.getThickness() * tempHess / 4.0;
			}
		}
		// if (*hess < 0)
		// (*hess) = -1 * (*hess);
		// consolidate
		return totalEnergy;
	}

	// 	Real TFWStretch::getVertexEnergy1(Sheet& sheet, int vid, Real* deriv, Real* hess) const
	// 	{
	// 		// need to compute deriv and hessian wrt the vertex at index vid
	// 		// just for testing
	// 		Real energy = 0;
	// 		const Eigen::VectorXi& faces = sheet.vbdInfo.vertexFacesSet.row(vid);
	// 		unsigned int nfaces = 0;
	// 		for (int i = 0; i < faces.rows(); i++)
	// 		{
	// 			if (faces(i) == -1)
	// 				break;
	// 			nfaces++;
	// 		}
	// 		// create buffers to store parallel results
	// 		std::vector<Real> energies(nfaces, 0);
	// 		std::vector<Real> derivatives(nfaces, 0);
	// 		if (deriv)
	// 			*deriv = 0;
	// 		if (hess)
	// 			*hess = 0;
	//
	// 		std::vector<Real> hessians(nfaces, 0);
	//
	// #ifdef USE_OMP
	// #pragma omp parallel for
	// #endif
	// 		for (int i = 0; i < nfaces; i++)
	// 		{
	// 			int faceId = faces(i);
	// 			int localVid = -1;
	// 			for (int ii = 0; ii < 3; ii++)
	// 			{
	// 				int vidii = sheet.getTriangles().at(faceId).v[ii];
	// 				if (vidii == vid)
	// 				{
	// 					localVid = ii;
	// 					break;
	// 				}
	// 			}
	// 			if (localVid == -1)
	// 			{
	// 				std::cout << "Vertex " << localVid << " not found! Face " << faceId << ". Error in tfw_stretching.cpp\n";
	// 				exit(0);
	// 			}
	// 			for (int quad = 0; quad < sheet.getQuadPoints().size(); quad++)
	// 			{
	// 				Real tempDeriv, tempHess;
	// 				energies[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * TFWStretch::getEnergyDensityFromQuadPerVtex1(localVid, faceId, quad, sheet, deriv ? &tempDeriv : NULL, hess ? &tempHess : NULL) / 4.0;
	// 				if (deriv)
	// 					derivatives[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * tempDeriv / 4.0;
	// 				if (hess)
	// 					hessians[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * tempHess / 4.0;
	// 			}
	// 		}
	//
	// 		// consolidate
	// 		for (int i = 0; i < nfaces; i++)
	// 		{
	// 			energy += energies[i];
	// 			if (deriv)
	// 				(*deriv) += derivatives[i];
	// 			if (hess)
	// 				(*hess) += hessians[i];
	// 		}
	// 		return energy;
	// 	}

	// function for per vertex derivatives and hessians amps + pos
	Real TFWStretch::getEnergyDensityFromQuadPerVtex4(int vid, int faceId, int quadId, Sheet& sheet, Vector4* deriv, Matrix44* hess) const
	{
		Real energy = 0;
		const std::vector<trimesh::triangle_t>& faces = sheet.getTriangles();
		const std::vector<trimesh::edge_t>& edges = sheet.getEdges();
		const Matrix3N& positions = sheet.getPositionsN3();
		const Eigen::MatrixXi& edgeOppositeVertices = sheet.getEdgeOppositeVertices();
		const std::vector<trimesh::triangle_t>& vertexOppositeEdges = sheet.getVertexOppositeEdgesPerFace();
		int nfaces = faces.size();
		int nedges = edges.size();
		int nvertices = positions.cols();
		Matrix33 facePositions;
		int vidxi = faces[faceId].v[0];
		int vidxj = faces[faceId].v[1];
		int vidxk = faces[faceId].v[2];
		// facePositions.col(0) = positions.block<3, 1>(0, faces[faceId].v[0]);
		// facePositions.col(1) = positions.block<3, 1>(0, faces[faceId].v[1]);
		// facePositions.col(2) = positions.block<3, 1>(0, faces[faceId].v[2]);

		facePositions.col(0) = positions.col(faces[faceId].v[0]);
		facePositions.col(1) = positions.col(faces[faceId].v[1]);
		facePositions.col(2) = positions.col(faces[faceId].v[2]);
		MatrixNN gradI, gradII; // gradI(i) contains dI(i)/dX
		std::vector<Eigen::Matrix<Real, 9, 9>> hessI;
		Matrix22 I;
		I = firstFundamentalForm(facePositions, true, true, &gradI, &hessI); // correct till here
		std::vector<Eigen::MatrixXd> hessII;
		std::shared_ptr<SecondFundamentalFormDiscretization> sff;
		sff = std::make_shared<MidedgeAverageFormulation>();
		Matrix22 II = sff->secondFundamentalForm(sheet.getBaseMesh(), sheet.getPositionsN3().transpose(), sheet.getEdgeDofs(), faceId, &gradII, &hessII);
		Matrix22 Ibar = sheet.getRestShape().getIbars()[faceId];
		Matrix22 Ibarinv = Ibar.inverse();
		Matrix22 IIbar = sheet.getRestShape().getIIbars()[faceId];
		Vector3 gradA; // datilda/da
		Vector2 da;
		Matrix23 gradDA; // dda/da
		Real u = sheet.getQuadPoints()[quadId].u;
		Real v = sheet.getQuadPoints()[quadId].v;
		Vector3 ampUVWeights;
		ampUVWeights << 1 - u - v, u, v;
		Real a = sheet.computeAmplitudesFromQuad(faceId, quadId, &da, &gradA, &gradDA, nullptr, nullptr);
		Matrix22 gradDphi;
		Vector2 dphi = sheet.computeDphi1VectorPerFace(faceId, &gradDphi);
		std::vector<Matrix22> gradDphiDphiTensor;
		std::vector<Matrix22> hessDphiDphiT;
		Matrix22 dphidphiTensor;
		dphidphiTensor = sheet.computeDphiDphi1VectorPerFaceTensor(faceId, quadId, &gradDphiDphiTensor, &hessDphiDphiT);
		std::vector<Matrix22> gradDaDaTensor;
		std::vector<Matrix33> hessDaDaTensor;
		Matrix22 dadaTensor;
		dadaTensor = sheet.computeDaDaTensor(faceId, quadId, &gradDaDaTensor, &hessDaDaTensor);
		std::vector<Matrix22> gradDaDphiTensor;
		std::vector<Eigen::Matrix<Real, 5, 5>> hessDaDphiTensor, hessDphiDaTensor;
		Matrix22 dadphiTensor = sheet.computeDaDphi1VectorPerFaceTensor(faceId, quadId, &gradDaDphiTensor, &hessDaDphiTensor);
		std::vector<Matrix22> gradDphiDaTensor;
		Matrix22 dphidaTensor = sheet.computeDphi1VectorPerFaceDaTensor(faceId, quadId, &gradDphiDaTensor, &hessDphiDaTensor);
		Real lameAlpha = sheet.getLameAlpha();
		Real lameBeta = sheet.getLameBeta();

		if (deriv)
			deriv->setZero();
		if (hess)
			hess->setZero();
		// First term
		// Correct numerically
		Matrix22 M0 = Ibarinv * (I + 0.5 * a * a * dphidphiTensor + 0.5 * dadaTensor - Ibar);
		Matrix22 dMSVdM;
		Real MSV;
		MSV = SVNormSquared(lameAlpha, lameBeta, M0, &dMSVdM);
		energy += MSV;
		Matrix22 dMda;
		if (deriv || hess)
			dMda = Ibarinv * (a * gradA(vid) * dphidphiTensor + 0.5 * gradDaDaTensor.at(vid));

		if (deriv)
		{
			// wrt amp
			(*deriv)(0) += doubleContraction(dMSVdM, dMda);
			// Derivatives wrt face positions
			for (int i = 0; i < 3; i++)
			{
				Matrix22 dIdxij;
				dIdxij << mat(gradI, vid, i);
				Matrix22 dM0dxij = Ibarinv * dIdxij;
				(*deriv)(i + 1) += doubleContraction(lameAlpha * M0.trace() * Matrix22::Identity() + 2 * lameBeta * M0, dM0dxij);
			}
		}

		if (hess)
		{
			// amps wrt amps
			Matrix22 hessDaDa;
			hessDaDa << hessMat(hessDaDaTensor, vid, vid);
			Matrix22 d2Mdada = Ibarinv * (gradA(vid) * gradA(vid) * dphidphiTensor + 0.5 * hessDaDa);
			// Matrix22 dMda = Ibarinv * (a * gradA(vid) * dphidphiTensor + 0.5 * gradDaDaTensor.at(vid));
			Matrix22 d2MSVdadM = lameAlpha * (dMda).trace() * Matrix22::Identity() + 2 * lameBeta * dMda;
			(*hess)(0, 0) += doubleContraction(d2MSVdadM, dMda) + doubleContraction(dMSVdM, d2Mdada);

			// pos vs pos
			for (unsigned int j = 0; j < 3; j++)
			{
				Matrix22 dIdxj;
				dIdxj << mat(gradI, vid, j);
				// pos vs amps
				Matrix22 dMdxj = Ibarinv * dIdxj;
				Real crossHess = doubleContraction(lameAlpha * dMdxj.trace() * Matrix22::Identity() + 2 * lameBeta * dMdxj, dMda);
				(*hess)(0, j + 1) += crossHess;
				(*hess)(j + 1, 0) += crossHess;
				// pos wrt pos
				for (unsigned int k = 0; k <= j; k++)
				{
					Matrix22 d2MSVdxjdM, dMdxk, d2Mdxjdxk;
					// d2MSVdxijdM
					d2MSVdxjdM = lameAlpha * doubleContraction(Ibarinv * dIdxj, Matrix22::Identity()) * Matrix22::Identity() + 2 * lameBeta * Ibarinv * dIdxj;
					// dMdxik
					Matrix22 dIdxk;
					dIdxk << mat(gradI, vid, k);
					dMdxk = Ibarinv * dIdxk;
					// d2Mdxijdxik
					Matrix22 hessIjk;
					hessIjk << hessMat(hessI, idxAt(vid, j), idxAt(vid, k));
					d2Mdxjdxk = Ibarinv * hessIjk;
					Real d2Edxjdxk = doubleContraction(d2MSVdxjdM, dMdxk)
					    + doubleContraction(dMSVdM, d2Mdxjdxk);
					(*hess)(j + 1, k + 1) += d2Edxjdxk;
					if (j != k)
						(*hess)(k + 1, j + 1) += d2Edxjdxk;
				}
			}
		}
		// Second term
		if (dphi.norm() > 1e-4)
		{
			Vector2 w = Ibarinv * dphi;
			Real num = w.norm() * w.norm() * II.trace() - w.transpose() * II * w;
			Real den = w.norm() * w.norm() * Ibar.trace() - w.transpose() * Ibar * w;
			Real c = 4 * lameBeta * ((lameAlpha + lameBeta) / (lameAlpha + 2 * lameBeta)) / (den * den);
			energy += a * a * c * num * num;
			Matrix22 dnumdII;
			if (deriv || hess)
				dnumdII = w.norm() * w.norm() * Matrix22::Identity() - w * w.transpose();
			if (deriv)
			{
				// wrt amps
				(*deriv)(0) += 2 * a * gradA(vid) * c * num * num;
				// wrt pos
				for (unsigned int j = 0; j < 3; j++)
				{
					Matrix22 dIIdxj;
					dIIdxj << mat(gradII, vid, j);
					(*deriv)(j + 1) += a * a * c * 2 * num * doubleContraction(dnumdII, dIIdxj);
				}
			}
			if (hess)
			{
				// amps vs amps
				(*hess)(0, 0) += 2 * gradA(vid) * gradA(vid) * c * num * num;
				// amps vs pos
				for (unsigned int j = 0; j < 3; j++)
				{
					Matrix22 dIIdxj;
					dIIdxj << mat(gradII, vid, j);

					Real h0j = 2 * gradA(vid) * a * c * 2 * num * doubleContraction(dnumdII, dIIdxj);
					(*hess)(0, j + 1) += h0j;
					(*hess)(j + 1, 0) += h0j;
				}
				// pos vs pos
				for (unsigned int j = 0; j < 3; j++)
				{
					Matrix22 dIIdxj;
					dIIdxj << mat(gradII, vid, j);
					for (unsigned int k = 0; k <= j; k++)
					{
						Real hjk = 0;
						Matrix22 dIIdxk;
						dIIdxk << mat(gradII, vid, k);
						Matrix22 d2IIdxjdxk;
						d2IIdxjdxk << hessMat(hessII, idxAt(vid, j), idxAt(vid, k));
						hjk = c * a * a * 2 * (doubleContraction(dnumdII, dIIdxj) * doubleContraction(dnumdII, dIIdxk) + num * doubleContraction(dnumdII, d2IIdxjdxk));
						(*hess)(j + 1, k + 1) += hjk;
						if (j != k)
							(*hess)(k + 1, j + 1) += hjk;
					}
				}
			}
		}

		// Third term
		// Numerically correct
		Matrix22 M2 = (Ibarinv * (dphidaTensor + dadphiTensor));
		Matrix22 dM2SVdM;
		Real M2SVM = SVNormSquared(sheet.getLameAlpha(), sheet.getLameBeta(), M2, &dM2SVdM);
		energy += (a * a / 32.0) * M2SVM;
		if (deriv)
		{
			// wrt amp
			Matrix22 DdphidaDa, DdadphiDa;
			DdphidaDa = gradDphiDaTensor.at(vid);
			DdadphiDa = gradDaDphiTensor.at(vid);
			// Matrix22 dMda = (gradA(vid) / std::sqrt(32.0)) * (Ibarinv * (dphidaTensor + dadphiTensor)) + (a / std::sqrt(32.0)) * (Matrix22 );
			// (*deriv)(0) += doubleContraction(dM2SVdM, dMda);
			(*deriv)(0) += (a * gradA(vid) / 16.0) * M2SVM + (a * a / 32.0) * doubleContraction(dM2SVdM, Ibarinv * (DdphidaDa + DdadphiDa));
		}
		if (hess)
		{
			Matrix22 DdphidaDa, DdadphiDa;
			DdphidaDa = gradDphiDaTensor.at(vid);
			DdadphiDa = gradDaDphiTensor.at(vid);
			Matrix22 dMda, d2MSVdadM;
			dMda = Ibarinv * (gradDphiDaTensor.at(vid) + gradDaDphiTensor.at(vid));
			d2MSVdadM = lameAlpha * dMda.trace() * Matrix22::Identity() + 2 * lameBeta * dMda;
			// wrt amps
			(*hess)(0, 0)
			    += (2 * gradA(vid) * gradA(vid) / 32.0) * M2SVM + 2 * (2 * gradA(vid) * a / 32.0) * doubleContraction(dM2SVdM, Ibarinv * (DdphidaDa + DdadphiDa))
			    + (a * a / 32.0) * doubleContraction(d2MSVdadM, dMda);
		}

		// Fourth term
		// Correct numerically
		Matrix22 M3 = Ibarinv * dadaTensor;
		Matrix22 dM3SVdM;
		energy += (1 / 8.0) * SVNormSquared(sheet.getLameAlpha(), sheet.getLameBeta(), M3, &dM3SVdM);
		if (deriv)
		{
			// wrt amps
			Matrix22 dMda = Ibarinv * gradDaDaTensor.at(vid);
			(*deriv)(0) += (1 / 8.0) * doubleContraction(dM3SVdM, dMda);
		}
		if (hess)
		{
			// wrt amps
			Matrix22 dMda = Ibarinv * gradDaDaTensor.at(vid);
			Matrix22 d2MSVdadM = lameAlpha * dMda.trace() * Matrix22::Identity() + 2 * lameBeta * dMda;
			Matrix22 hessDaDa;
			hessDaDa << hessMat(hessDaDaTensor, vid, vid);
			(*hess)(0, 0) += (1 / 8.0) * (doubleContraction(d2MSVdadM, dMda) + doubleContraction(dM3SVdM, Ibarinv * hessDaDa));
		}

		// project
		// if (hess)
		// 	*hess = argus::util::hessianProjection((*hess), true);
		Real area
		    = 0.5 * sqrt(Ibar.determinant());
		if (deriv)
			*deriv = (area * (*deriv)).eval();
		if (hess)
			(*hess) = (area * (*hess)).eval();
		return area * energy;
	}

	Real TFWStretch::getVertexEnergy4New(Sheet& sheet, int vid, Vector4* deriv, Matrix44* hess) const
	{
		// need to compute deriv and hessian wrt the vertex at index vid
		// just for testing
		Real energy = 0;
		const Eigen::VectorXi& faces = sheet.sheetInfo.vertexFacesSet.row(vid);
		unsigned int nfaces = 0;
		for (int i = 0; i < faces.rows(); i++)
		{
			if (faces(i) == -1)
				break;
			nfaces++;
		}
		// create buffers to store parallel results
		std::vector<Real> energies(nfaces, 0);
		std::vector<Vector4> derivatives(nfaces);
		if (deriv)
			deriv->setZero();
		if (hess)
			hess->setZero();
		if (deriv)
			for (int i = 0; i < nfaces; i++)
			{
				derivatives[i].setZero();
			}
		std::vector<Matrix44> hessians(nfaces);
		if (hess)
			for (int i = 0; i < nfaces; i++)
			{
				hessians[i].setZero();
			}

#ifdef USE_OMP
#pragma omp parallel for
#endif
		for (int i = 0; i < nfaces; i++)
		{
			int faceId = faces(i);
			int localVid = -1;
			for (int ii = 0; ii < 3; ii++)
			{
				int vidii = sheet.getTriangles().at(faceId).v[ii];
				if (vidii == vid)
				{
					localVid = ii;
					break;
				}
			}
			if (localVid == -1)
			{
				std::cout << "Vertex " << localVid << " not found! Face " << faceId << ". Error in tfw_stretching.cpp\n";
				exit(0);
			}
			for (int quad = 0; quad < sheet.getQuadPoints().size(); quad++)
			{
				Vector4 tempDeriv1, tempDeriv2, tempDeriv;
				Matrix44 tempHess, tempHess1, tempHess2;
				/* Vector3 td1, td2, td;
				Matrix33 th;
				// test hessians
				Vector3 dx = Vector3::Random().normalized();
				// dx(0) = 0;
				Real h = 0.1;
				sheet.getPositionsN3().col(vid) = argus::Vector3::Random() * 99;
				Vector3 X = sheet.getPositionsN3().col(vid);
				sheet.getAmplitudes()(vid) = 99;
				Real ainit = sheet.getAmplitudes()(vid);
				TFWStretch::getEnergyDensityFromQuadPerVtex(localVid, faceId, quad, sheet, &td, &th);
				// std::cout << "Temphess\n"
				//           << tempHess << std::endl;
				for (int ii = 0; ii < 10; ii++)
				{
				    Vector3 hdx = X + 0.5 * h * dx;
				    sheet.getPositionsN3().col(vid) = hdx;
				    // sheet.getAmplitudes()(vid) = ainit + 0.5 * h * dx(0);
				    Real e1 = TFWStretch::getEnergyDensityFromQuadPerVtex(localVid, faceId, quad, sheet, &td1, &th);
				    hdx = X - 0.5 * h * dx;
				    sheet.getPositionsN3().col(vid) = hdx;
				    // sheet.getAmplitudes()(vid) = ainit - 0.5 * h * dx(0);
				    Real e2 = TFWStretch::getEnergyDensityFromQuadPerVtex(localVid, faceId, quad, sheet, &td2, &th);
				    // std::cout << (tempDeriv1 - tempDeriv2).transpose()/h << " - " << (tempHess*dx).transpose() << std::endl;
				    // std::cout << "ErrorHessians " << (((tempDeriv1 - tempDeriv2) / h) - tempHess * dx).norm() << std::endl;
				    std::cout << "Error deriv " << std::fabs(((e1 - e2) / h) - dx.dot(td)) << std::endl;
				    h /= 10.0;
				}
				exit(0); */
				energies[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * TFWStretch::getEnergyDensityFromQuadPerVtex4(localVid, faceId, quad, sheet, deriv ? &tempDeriv : NULL, hess ? &tempHess : NULL) / 4.0;
				if (deriv)
					derivatives[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * tempDeriv / 4.0;
				if (hess)
					hessians[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * tempHess / 4.0;
			}
		}

		// consolidate
		for (int i = 0; i < nfaces; i++)
		{
			energy += energies[i];
			if (deriv)
				(*deriv) += derivatives[i];
			if (hess)
				(*hess) += hessians[i];
		}
		return energy;
	}

	Real TFWStretch::getVertexEnergy4(Sheet& sheet, int vid, Vector4* deriv, Matrix44* hess) const
	{
		// need to compute deriv and hessian wrt the vertex at index vid
		// just for testing
		Real totalEnergy = 0;
		Real lameAlpha = sheet.getLameAlpha();
		Real lameBeta = sheet.getLameBeta();

		const Eigen::Vector<int, 32>& faces = sheet.sheetInfo.vertexFacesSet.row(vid);
		unsigned int nfaces = 0;
		for (int i = 0; i < faces.rows(); i++)
		{
			if (faces(i) == -1)
				break;
			nfaces++;
		}
		if (deriv)
			deriv->setZero();
		if (hess)
			hess->setZero();

		for (int i = 0; i < nfaces; i++)
		{
			int faceId = faces(i);
			int localVid = -1;
			for (int ii = 0; ii < 3; ii++)
			{
				int vidii = sheet.getTriangles().at(faceId).v[ii];
				if (vidii == vid)
				{
					localVid = ii;
					break;
				}
			}
			if (localVid == -1)
			{
				std::cout << "Vertex " << localVid << " not found! Face " << faceId << ". Error in tfw_stretching.cpp\n";
				exit(0);
			}

			const std::vector<trimesh::triangle_t>& faces = sheet.getTriangles();
			const std::vector<trimesh::edge_t>& edges = sheet.getEdges();
			const Matrix3N& positions = sheet.getPositionsN3();
			const Eigen::MatrixXi& edgeOppositeVertices = sheet.getEdgeOppositeVertices();
			const std::vector<trimesh::triangle_t>& vertexOppositeEdges = sheet.getVertexOppositeEdgesPerFace();

			Matrix33 facePositions;
			int vidxi = faces[faceId].v[0];
			int vidxj = faces[faceId].v[1];
			int vidxk = faces[faceId].v[2];

			facePositions.col(0) = positions.col(faces[faceId].v[0]);
			facePositions.col(1) = positions.col(faces[faceId].v[1]);
			facePositions.col(2) = positions.col(faces[faceId].v[2]);
			Eigen::Matrix<Real, 4, 9> gradI;
			std::array<Eigen::Matrix<Real, 9, 9>, 4> hessI;
			Matrix22 II;
			Matrix22 I = firstFundamentalFormNew(facePositions, true, true, &gradI, &hessI); // correct till here

			Matrix22 Ibar = sheet.getRestShape().getIbars()[faceId];
			Matrix22 Ibarinv = Ibar.inverse();
			Matrix22 IIbar = sheet.getRestShape().getIIbars()[faceId];
			Vector2 dphi = sheet.computeDphi1VectorPerFace(faceId, NULL);
			Eigen::Matrix<Real, 4, 18> gradII; // gradI(i) contains dI(i)/dX
			std::array<Eigen::Matrix<Real, 18, 18>, 4> hessII;
			// if (dphi.norm() > 0)
			// {
			II.noalias() = sffobject->secondFundamentalFormNew(sheet.getBaseMesh(), sheet.getPositionsN3(), sheet.getEdgeDofs(), faceId, &gradII, &hessII);
			// }

			Matrix22 dphidphiTensor;

			dphidphiTensor.noalias() = sheet.computeDphiDphi1VectorPerFaceTensorNew(faceId, 0, NULL, NULL);
			std::array<Matrix22, 3> gradDaDaTensor;
			std::array<Matrix33, 4> hessDaDaTensor;
			Matrix22 dadaTensor;
			std::array<Matrix22, 5> gradDaDphiTensor;
			std::array<Matrix22, 5> gradDphiDaTensor;
			Real lameAlpha = sheet.getLameAlpha();
			Real lameBeta = sheet.getLameBeta();
			dadaTensor = sheet.computeDaDaTensorNew(faceId, 0, &gradDaDaTensor, &hessDaDaTensor);

			Matrix22 dadphiTensor = sheet.computeDaDphi1VectorPerFaceTensorNew(faceId, 0, &gradDaDphiTensor, NULL);
			Matrix22 dphidaTensor = sheet.computeDphi1VectorPerFaceDaTensorNew(faceId, 0, &gradDphiDaTensor, NULL);

			Matrix22 tempM = (Ibarinv * (dphidaTensor + dadphiTensor));
			Matrix22 dM2SVdM;
			Real M2norm = SVNormSquared(sheet.getLameAlpha(), sheet.getLameBeta(), tempM, &dM2SVdM);
			tempM.noalias() = Ibarinv * dadaTensor;
			Matrix22 dM3SVdM;
			Real de4, he4;
			Real M3norm = SVNormSquared(sheet.getLameAlpha(), sheet.getLameBeta(), tempM, &dM3SVdM) / 8.0;

			if (deriv)
			{
				// wrt amps
				Matrix22 dMda = Ibarinv * gradDaDaTensor[localVid];
				de4 = (1 / 8.0) * doubleContraction(dM3SVdM, dMda);
			}
			if (hess)
			{
				// wrt amps
				Matrix22 dMda = Ibarinv * gradDaDaTensor[localVid];
				Matrix22 d2MSVdadM = lameAlpha * dMda.trace() * Matrix22::Identity() + 2 * lameBeta * dMda;
				Matrix22 hessDaDa;
				hessDaDa << hessMat(hessDaDaTensor, localVid, localVid);
				he4 = (1 / 8.0) * (doubleContraction(d2MSVdadM, dMda) + doubleContraction(dM3SVdM, Ibarinv * hessDaDa));
			}
			Real area = 0.5 * sqrt(Ibar.determinant());
			Vector2 w = Ibarinv * dphi;
			Real wnorm = w.norm();
			Real coeff = 4 * lameBeta * (lameAlpha + lameBeta) / (lameAlpha + 2 * lameBeta);
			coeff = coeff / (wnorm * wnorm * Ibar.trace() - w.dot(Ibar * w));
			coeff = coeff / (wnorm * wnorm * Ibar.trace() - w.dot(Ibar * w));
			Real W2s = wnorm * wnorm * II.trace() - w.dot(II * w);
			W2s = W2s * W2s;
			for (int quadId = 0; quadId < sheet.getQuadPoints().size(); quadId++)
			{
				Vector4 tempDeriv;
				Matrix44 tempHess;
				tempDeriv.setZero();
				tempHess.setZero();
				Real energy = 0;
				Vector3 gradA;
				Real a = sheet.computeAmplitudesFromQuadNew(faceId, quadId, NULL, &gradA, NULL, nullptr, nullptr);

				if (dphi.norm() > 1e-4)
				{

					energy += a * a * coeff * W2s;

					Eigen::Matrix<Real, 2, 2> W2sderiv;
					if (deriv || hess)
						W2sderiv = 2 * a * a * coeff * (wnorm * wnorm * Matrix22::Identity() - w * w.transpose()) * (wnorm * wnorm * II.trace() - w.dot(II * w));

					// positions
					if (deriv)
					{
						tempDeriv(0) += 2 * a * gradA(localVid) * coeff * W2s;
						// Matrix22 dW2sdII = W2sderiv; //.reshaped(2, 2);
						Matrix22 dIIdxi;
						for (int i = 0; i < 3; i++)
						{
							int vidComp = idxAt(localVid, i);
							dIIdxi << gradII(0, vidComp), gradII(1, vidComp), gradII(2, vidComp), gradII(3, vidComp);
							tempDeriv(i + 1) += doubleContraction(W2sderiv, dIIdxi);
						}
					}
					if (hess)
					{
						Matrix22 temp = w.norm() * w.norm() * Matrix22::Identity() - w * w.transpose();

						tempHess(0, 0) += 2 * gradA(localVid) * gradA(localVid) * coeff * W2s;
						for (int i = 0; i < 3; i++)
						{
							int vidComp = idxAt(localVid, i);
							Matrix22 dIIdxi;
							dIIdxi << gradII(0, vidComp), gradII(1, vidComp), gradII(2, vidComp), gradII(3, vidComp);
							Real tempval = 2 * a * gradA(localVid) * coeff * 2 * (wnorm * wnorm * II.trace() - w.dot(II * w)) * doubleContraction(temp, dIIdxi);
							tempHess(0, i + 1) += tempval;
							tempHess(i + 1, 0) += tempval;
						}
						// for (unsigned int j = 0; j < 3; j++)
						// {
						// 	Matrix22 dIIdxj;
						// 	dIIdxj << mat(gradII, localVid, j);
						// 	Real h0j = 2 * gradA(localVid) * a * coeff * 2 * (wnorm * wnorm * II.trace() - w.dot(II * w)) * doubleContraction(temp, dIIdxj);
						// 	tempHess(0, j + 1) += h0j;
						// 	tempHess(j + 1, 0) += h0j;
						// }
						for (int j = 0; j < 3; j++)
						{
							Matrix22 d2EdxijdII, dIIdxij;
							dIIdxij << mat(gradII, localVid, j);
							d2EdxijdII = 2 * a * a * coeff * (temp)*doubleContraction(temp, dIIdxij);
							for (int k = 0; k <= j; k++)
							{
								Matrix22 d2IIdxjdxk, dIIdxik;
								dIIdxik << mat(gradII, localVid, k);
								d2IIdxjdxk << hessMat(hessII, idxAt(localVid, j), idxAt(localVid, k));
								Real reqHess = doubleContraction(W2sderiv, d2IIdxjdxk) + doubleContraction(d2EdxijdII, dIIdxik);
								tempHess(j + 1, k + 1) += reqHess;
								if (j != k)
									tempHess(k + 1, j + 1) += reqHess;
							}
						}
					}
				}

				// First term
				// Correct numerically
				Matrix22 M0 = Ibarinv * (I + 0.5 * a * a * dphidphiTensor + 0.5 * dadaTensor - Ibar);
				Matrix22 dMSVdM;
				energy += SVNormSquared(lameAlpha, lameBeta, M0, &dMSVdM);
				Matrix22 dMda;
				if (deriv || hess)
					dMda = Ibarinv * (a * gradA(localVid) * dphidphiTensor + 0.5 * gradDaDaTensor[localVid]);

				if (deriv)
				{
					// wrt amp
					tempDeriv(0) += doubleContraction(dMSVdM, dMda);

					// Derivatives wrt face positions
					for (int i = 0; i < 3; i++)
					{
						Matrix22 dIdxij;
						dIdxij << mat(gradI, localVid, i);
						Matrix22 dM0dxij = Ibarinv * dIdxij;
						tempDeriv(i + 1) += doubleContraction(lameAlpha * M0.trace() * Matrix22::Identity() + 2 * lameBeta * M0, dM0dxij);
					}
					// d || M || dxi = d || M || dM : Ibarinv dIdxi
				}
				if (hess)
				{
					// amps wrt amps
					Matrix22 hessDaDa;
					hessDaDa << hessMat(hessDaDaTensor, localVid, localVid);
					Matrix22 d2Mdada = Ibarinv * (gradA(localVid) * gradA(localVid) * dphidphiTensor + 0.5 * hessDaDa);
					// Matrix22 dMda = Ibarinv * (a * gradA(localVid) * dphidphiTensor + 0.5 * gradDaDaTensor.at(localVid));
					Matrix22 d2MSVdadM = lameAlpha * (dMda).trace() * Matrix22::Identity() + 2 * lameBeta * dMda;
					tempHess(0, 0) += doubleContraction(d2MSVdadM, dMda) + doubleContraction(dMSVdM, d2Mdada);

					// cross terms
					Matrix22 dIdxj;

					// TODO: Add optimizations
					for (unsigned int j = 0; j < 3; j++)
					{
						dIdxj << mat(gradI, localVid, j);
						Matrix22 dMdxj = Ibarinv * dIdxj;

						Real crossHess = doubleContraction(lameAlpha * dMdxj.trace() * Matrix22::Identity() + 2 * lameBeta * dMdxj, dMda);

						tempHess(0, j + 1) += crossHess;
						tempHess(j + 1, 0) += crossHess;
					}
					for (unsigned int j = 0; j < 3; j++)
					{

						Matrix22 dIdxij;
						dIdxij << mat(gradI, localVid, j);
						// d2MSVdxijdM
						Matrix22 d2MSVdxijdM = lameAlpha * doubleContraction(Ibarinv * dIdxij, Matrix22::Identity()) * Matrix22::Identity() + 2 * lameBeta * Ibarinv * dIdxij;
						for (unsigned int k = 0; k <= j; k++)
						{
							Matrix22 dMdxik, d2Mdxijdxik;
							// dMdxik
							Matrix22 dIdxik;
							dIdxik << mat(gradI, localVid, k);
							dMdxik.noalias() = Ibarinv * dIdxik;
							// d2Mdxijdxik
							Matrix22 hessIjk;
							hessIjk << hessMat(hessI, idxAt(localVid, j), idxAt(localVid, k));
							d2Mdxijdxik.noalias() = Ibarinv * hessIjk;
							Real d2Edxijdxik = doubleContraction(d2MSVdxijdM, dMdxik)
							    + doubleContraction(dMSVdM, d2Mdxijdxik);
							tempHess(j + 1, k + 1) += d2Edxijdxik;
							if (k != j)
								tempHess(k + 1, j + 1) += d2Edxijdxik;
						}
					}
				}
				// Third term
				// Numerically correct
				energy += (a * a / 32.0) * M2norm;
				if (deriv)
				{
					tempDeriv(0) += (a * gradA(localVid) / 16.0) * M2norm + (a * a / 32.0) * doubleContraction(dM2SVdM, Ibarinv * (gradDphiDaTensor[localVid] + gradDaDphiTensor[localVid]));
				}
				if (hess)
				{
					Matrix22 dMda, d2MSVdadM;
					dMda = Ibarinv * (gradDphiDaTensor[localVid] + gradDaDphiTensor[localVid]);
					// d2MSVdadM = lameAlpha * dMda.trace() * Matrix22::Identity() + 2 * lameBeta * dMda;
					// wrt amps
					tempHess(0, 0)
					    += (2 * gradA(localVid) * gradA(localVid) / 32.0) * M2norm + 2 * (2 * gradA(localVid) * a / 32.0) * doubleContraction(dM2SVdM, Ibarinv * (gradDphiDaTensor[localVid] + gradDaDphiTensor[localVid]))
					    + (a * a / 32.0) * doubleContraction(lameAlpha * dMda.trace() * Matrix22::Identity() + 2 * lameBeta * dMda, dMda);
				}
				// Fourth term
				// Correct numerically
				energy += M3norm;
				if (deriv)
					tempDeriv(0) += de4;
				if (hess)
					tempHess(0, 0) += he4;
				if (deriv)
					(*deriv) += sheet.getQuadPoints()[quadId].weight * area * sheet.getThickness() / 4.0 * tempDeriv;
				if (hess)
					(*hess) += sheet.getQuadPoints()[quadId].weight * area * sheet.getThickness() / 4.0 * tempHess;
				totalEnergy += sheet.getQuadPoints()[quadId].weight * area * sheet.getThickness() / 4.0 * energy;
			}
		}

		return totalEnergy;
	}

	Real TFWStretch::getEnergyDensityFromQuadPerVtex(int vid, int faceId, int quadId, Sheet& sheet, Vector3* deriv, Matrix33* hess) const
	{
		Real energy = 0;
		const std::vector<trimesh::triangle_t>& faces = sheet.getTriangles();
		const std::vector<trimesh::edge_t>& edges = sheet.getEdges();
		const Matrix3N& positions = sheet.getPositionsN3();
		const Eigen::MatrixXi& edgeOppositeVertices = sheet.getEdgeOppositeVertices();
		const std::vector<trimesh::triangle_t>& vertexOppositeEdges = sheet.getVertexOppositeEdgesPerFace();
		int nfaces = faces.size();
		int nedges = edges.size();
		int nvertices = positions.cols();
		Matrix33 facePositions;
		int vidxi = faces[faceId].v[0];
		int vidxj = faces[faceId].v[1];
		int vidxk = faces[faceId].v[2];
		// facePositions.col(0) = positions.block<3, 1>(0, faces[faceId].v[0]);
		// facePositions.col(1) = positions.block<3, 1>(0, faces[faceId].v[1]);
		// facePositions.col(2) = positions.block<3, 1>(0, faces[faceId].v[2]);
		//
		facePositions.col(0) = positions.col(faces[faceId].v[0]);
		facePositions.col(1) = positions.col(faces[faceId].v[1]);
		facePositions.col(2) = positions.col(faces[faceId].v[2]);
		Eigen::Matrix<Real, 4, 9> gradI;
		std::array<Eigen::Matrix<Real, 9, 9>, 4> hessI;
		Matrix22 I;
		I = firstFundamentalFormNew(facePositions, true, true, &gradI, &hessI); // correct till here
		Matrix22 Ibar = sheet.getRestShape().getIbars()[faceId];
		Matrix22 Ibarinv = Ibar.inverse();
		Matrix22 IIbar = sheet.getRestShape().getIIbars()[faceId];
		Vector3 gradA; // datilda/da
		Vector2 da;
		Matrix23 gradDA; // dda/da
		Real u = sheet.getQuadPoints()[quadId].u;
		Real v = sheet.getQuadPoints()[quadId].v;
		Vector3 ampUVWeights;
		ampUVWeights << 1 - u - v, u, v;
		Real a = sheet.computeAmplitudesFromQuadNew(faceId, quadId, &da, &gradA, &gradDA, nullptr, nullptr);
		Matrix22 gradDphi;
		Vector2 dphi = sheet.computeDphi1VectorPerFace(faceId, &gradDphi);
		Eigen::Matrix<Real, 4, 18> gradII; // gradI(i) contains dI(i)/dX
		std::array<Eigen::Matrix<Real, 18, 18>, 4> hessII;
		Matrix22 II;
		Matrix22 dphidphiTensor;
		dphidphiTensor = sheet.computeDphiDphi1VectorPerFaceTensorNew(faceId, quadId, NULL, NULL);
		Matrix22 dadaTensor;
		dadaTensor = sheet.computeDaDaTensor(faceId, quadId, NULL, NULL);
		Matrix22 dadphiTensor = sheet.computeDaDphi1VectorPerFaceTensor(faceId, quadId, NULL, NULL);
		Matrix22 dphidaTensor = sheet.computeDphi1VectorPerFaceDaTensor(faceId, quadId, NULL, NULL);
		Real lameAlpha = sheet.getLameAlpha();
		Real lameBeta = sheet.getLameBeta();

		if (deriv)
			deriv->setZero();
		if (hess)
			hess->setZero();

		// Second term
		// Derivatives are correct numerically
		// 		if (dphi.norm() > 1e-4)
		// 		{
		// 			std::shared_ptr<SecondFundamentalFormDiscretization> sff;
		// 			sff = std::make_shared<MidedgeAverageFormulation>();
		// 			II.noalias() = sff->secondFundamentalFormNew(sheet.getBaseMesh(), sheet.getPositionsN3(), sheet.getEdgeDofs(), faceId, &gradII, &hessII);

		// #ifdef ARGUS_USE_FLOATS
		// 			typedef TinyAD::Float<4> AReal;
		// #else
		// 			typedef TinyAD::Double<4> AReal;
		// #endif
		// 			auto W2sFunc = [vid, u, v, &Ibarinv, &Ibar, &lameAlpha, &lameBeta](Eigen::Vector<Real, 9> X, Eigen::Vector<Real, 4>* deriv = NULL, Eigen::Matrix<Real, 4, 4>* hess = NULL)
		// 			{
		// 				Eigen::Vector<AReal, 4> activeX = AReal::make_active(X(Eigen::seq(5, 8)));
		// 				AReal a = (1 - u - v) * X(0) + u * X(1) + v * X(2);
		// 				Eigen::Vector<AReal, 2> w = Ibarinv * X(Eigen::seq(3, 4));
		// 				Real coeff = 4 * lameBeta * (lameAlpha + lameBeta) / (lameAlpha + 2 * lameBeta);
		// 				Eigen::Matrix<AReal, 2, 2> II = activeX.reshaped(2, 2);
		// 				AReal num = w.norm() * w.norm() * II.trace() - w.transpose() * II * w;
		// 				AReal den = w.norm() * w.norm() * Ibar.trace() - w.transpose() * Ibar * w;
		// 				AReal e = coeff * a * a * (num / den) * (num / den);
		// 				if (deriv)
		// 					(*deriv) = e.grad;
		// 				if (hess)
		// 					(*hess) = e.Hess;
		// 				return TinyAD::to_passive(e);
		// 			};
		// 			Eigen::Vector<Real, 9> W2sX;
		// 			W2sX(0) = sheet.getAmplitudes()(vidxi);
		// 			W2sX(1) = sheet.getAmplitudes()(vidxj);
		// 			W2sX(2) = sheet.getAmplitudes()(vidxk);
		// 			W2sX(3) = dphi(0);
		// 			W2sX(4) = dphi(1);
		// 			for (int i = 0; i < 2; i++)
		// 			{
		// 				for (int j = 0; j < 2; j++)
		// 				{
		// 					W2sX(i * 2 + j + 5) = II(i, j);
		// 				}
		// 			}

		// 			Eigen::Vector<Real, 4> W2sderiv;
		// 			Eigen::Matrix<Real, 4, 4> W2shess;
		// 			energy += W2sFunc(W2sX, deriv ? &W2sderiv : NULL, hess ? &W2shess : NULL);
		// 			// positions
		// 			if (deriv)
		// 			{
		// 				Matrix22 dW2sdII = W2sderiv.reshaped(2, 2);
		// 				Matrix22 dIIdxi;
		// 				for (int i = 0; i < 3; i++)
		// 				{
		// 					int vidComp = idxAt(vid, i);
		// 					dIIdxi << gradII(0, vidComp), gradII(1, vidComp), gradII(2, vidComp), gradII(3, vidComp);
		// 					(*deriv)(i) += doubleContraction(dW2sdII, dIIdxi);
		// 				}
		// 			}
		// 			if (hess)
		// 			{
		// 				Matrix22 dW2sdII = W2sderiv.reshaped(2, 2);
		// 				Vector2 w = Ibarinv * dphi;
		// 				Real num = w.norm() * w.norm() * II.trace() - w.transpose() * II * w;
		// 				Real den = w.norm() * w.norm() * Ibar.trace() - w.transpose() * Ibar * w;
		// 				Real tvar = (w.norm() * w.norm() * Ibar.trace() - w.transpose() * Ibar * w);
		// 				Real c = 4 * a * a * lameBeta * ((lameAlpha + lameBeta) / (lameAlpha + 2 * lameBeta)) / (tvar);
		// 				c = c / tvar;
		// 				for (int j = 0; j < 3; j++)
		// 				{
		// 					Matrix22 d2EdxijdII, dIIdxij;
		// 					Matrix22 temp = w.norm() * w.norm() * Matrix22::Identity() - w * w.transpose();
		// 					dIIdxij << mat(gradII, vid, j);
		// 					d2EdxijdII = 2 * c * (temp)*doubleContraction(temp, dIIdxij);
		// 					for (int k = 0; k <= j; k++)
		// 					{
		// 						Matrix22 d2IIdxjdxk, dIIdxik;
		// 						dIIdxik << mat(gradII, vid, k);
		// 						d2IIdxjdxk << hessMat(hessII, idxAt(vid, j), idxAt(vid, k));
		// 						Real reqHess = doubleContraction(dW2sdII, d2IIdxjdxk) + doubleContraction(d2EdxijdII, dIIdxik);
		// 						(*hess)(j, k) += reqHess;
		// 						if (j != k)
		// 							(*hess)(k, j) += reqHess;
		// 					}
		// 				}
		// 			}
		// 		}
		if (dphi.norm() > 1e-4)
		{
			std::shared_ptr<SecondFundamentalFormDiscretization> sff;
			sff = std::make_shared<MidedgeAverageFormulation>();
			II.noalias() = sff->secondFundamentalFormNew(sheet.getBaseMesh(), sheet.getPositionsN3(), sheet.getEdgeDofs(), faceId, &gradII, &hessII);
			Vector2 w = Ibarinv * dphi;
			Real wnorm = w.norm();
			Real coeff = 4 * a * a * lameBeta * (lameAlpha + lameBeta) / (lameAlpha + 2 * lameBeta);
			coeff = coeff / (wnorm * wnorm * Ibar.trace() - w.dot(Ibar * w));
			coeff = coeff / (wnorm * wnorm * Ibar.trace() - w.dot(Ibar * w));
			Real W2s = wnorm * wnorm * II.trace() - w.dot(II * w);
			W2s = W2s * W2s;
			energy += coeff * W2s;
			Eigen::Matrix<Real, 2, 2> W2sderiv;
			if (deriv || hess)
				W2sderiv = 2 * coeff * (wnorm * wnorm * Matrix22::Identity() - w * w.transpose()) * (wnorm * wnorm * II.trace() - w.dot(II * w));

			// positions
			if (deriv)
			{
				Matrix22 dW2sdII = W2sderiv.reshaped(2, 2);
				Matrix22 dIIdxi;
				for (int i = 0; i < 3; i++)
				{
					int vidComp = idxAt(vid, i);
					dIIdxi << gradII(0, vidComp), gradII(1, vidComp), gradII(2, vidComp), gradII(3, vidComp);
					(*deriv)(i) += doubleContraction(dW2sdII, dIIdxi);
				}
			}
			if (hess)
			{
				// Vector2 w = Ibarinv * dphi;
				// Real num = w.norm() * w.norm() * II.trace() - w.transpose() * II * w;
				// Real den = w.norm() * w.norm() * Ibar.trace() - w.transpose() * Ibar * w;
				// Real tvar = (w.norm() * w.norm() * Ibar.trace() - w.transpose() * Ibar * w);
				// Real c = 4 * a * a * lameBeta * ((lameAlpha + lameBeta) / (lameAlpha + 2 * lameBeta)) / (tvar);
				// c = c / tvar;
				for (int j = 0; j < 3; j++)
				{
					Matrix22 d2EdxijdII, dIIdxij;
					Matrix22 temp = w.norm() * w.norm() * Matrix22::Identity() - w * w.transpose();
					dIIdxij << mat(gradII, vid, j);
					d2EdxijdII = 2 * coeff * (temp)*doubleContraction(temp, dIIdxij);
					for (int k = 0; k <= j; k++)
					{
						Matrix22 d2IIdxjdxk, dIIdxik;
						dIIdxik << mat(gradII, vid, k);
						d2IIdxjdxk << hessMat(hessII, idxAt(vid, j), idxAt(vid, k));
						Real reqHess = doubleContraction(W2sderiv, d2IIdxjdxk) + doubleContraction(d2EdxijdII, dIIdxik);
						(*hess)(j, k) += reqHess;
						if (j != k)
							(*hess)(k, j) += reqHess;
					}
				}
			}
		}

		// First term
		// Correct numerically
		Matrix22 M0 = Ibarinv * (I + 0.5 * a * a * dphidphiTensor + 0.5 * dadaTensor - Ibar);
		Matrix22 dMSVdM;
		Real MSV;
		MSV = SVNormSquared(lameAlpha, lameBeta, M0, &dMSVdM);
		energy += MSV;
		if (deriv)
		{
			// Derivatives wrt face positions
			for (int i = 0; i < 3; i++)
			{
				Matrix22 dIdxij;
				dIdxij << mat(gradI, vid, i);
				Matrix22 dM0dxij = Ibarinv * dIdxij;
				(*deriv)(i) += doubleContraction(lameAlpha * M0.trace() * Matrix22::Identity() + 2 * lameBeta * M0, dM0dxij);
			}
		}
		if (hess)
		{
			for (unsigned int j = 0; j < 3; j++)
			{
				Matrix22 dIdxij;
				dIdxij << mat(gradI, vid, j);
				for (unsigned int k = 0; k <= j; k++)
				{
					Matrix22 d2MSVdxijdM, dMdxik, d2Mdxijdxik;
					// d2MSVdxijdM
					d2MSVdxijdM = lameAlpha * doubleContraction(Ibarinv * dIdxij, Matrix22::Identity()) * Matrix22::Identity() + 2 * lameBeta * Ibarinv * dIdxij;
					// dMdxik
					Matrix22 dIdxik;
					dIdxik << mat(gradI, vid, k);
					dMdxik = Ibarinv * dIdxik;
					// d2Mdxijdxik
					Matrix22 hessIjk;
					hessIjk << hessMat(hessI, idxAt(vid, j), idxAt(vid, k));
					d2Mdxijdxik = Ibarinv * hessIjk;
					Real d2Edxijdxik = doubleContraction(d2MSVdxijdM, dMdxik)
					    + doubleContraction(dMSVdM, d2Mdxijdxik);
					(*hess)(j, k) += d2Edxijdxik;
					if (k != j)
						(*hess)(k, j) += d2Edxijdxik;
				}
			}
		}
		// Third term
		// Numerically correct
		Matrix22 M2 = (Ibarinv * (dphidaTensor + dadphiTensor));
		energy += (a * a / 32.0) * SVNormSquared(sheet.getLameAlpha(), sheet.getLameBeta(), M2);

		// Fourth term
		// Correct numerically
		Matrix22 M3 = Ibarinv * dadaTensor;
		energy += (1.0 / 8.0) * SVNormSquared(sheet.getLameAlpha(), sheet.getLameBeta(), M3);

		Real area = 0.5 * sqrt(Ibar.determinant());
		if (deriv)
			*deriv = (area * (*deriv)).eval();
		if (hess)
			(*hess) = (area * (*hess)).eval();

		// project

		return area * energy;
	}

	Real TFWStretch::getVertexEnergyNew(Sheet& sheet, int vid, Vector3* deriv, Matrix33* hess) const
	{
		// need to compute deriv and hessian wrt the vertex at index vid
		// just for testing
		Real energy = 0;
		const Eigen::Vector<int, 32>& faces = sheet.sheetInfo.vertexFacesSet.row(vid);
		unsigned int nfaces = 0;
		for (int i = 0; i < faces.rows(); i++)
		{
			if (faces(i) == -1)
				break;
			nfaces++;
		}
		// create buffers to store parallel results
		if (deriv)
			deriv->setZero();
		if (hess)
			hess->setZero();
		for (int i = 0; i < nfaces; i++)
		{
			int faceId = faces(i);
			int localVid = -1;
			for (int ii = 0; ii < 3; ii++)
			{
				int vidii = sheet.getTriangles().at(faceId).v[ii];
				if (vidii == vid)
				{
					localVid = ii;
					break;
				}
			}
			if (localVid == -1)
			{
				std::cout << "Vertex " << localVid << " not found! Face " << faceId << ". Error in tfw_stretching.cpp\n";
				throw "Error!";
			}
			for (int quad = 0; quad < sheet.getQuadPoints().size(); quad++)
			{
				Vector3 tempDeriv;
				Matrix33 tempHess;
				energy += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * TFWStretch::getEnergyDensityFromQuadPerVtex(localVid, faceId, quad, sheet, deriv ? &tempDeriv : NULL, hess ? &tempHess : NULL) / 4.0;
				if (deriv)
					(*deriv) += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * tempDeriv / 4.0;
				if (hess)
					(*hess) += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * tempHess / 4.0;
			}
		}

		// if (hess)
		// 	*hess = argus::util::hessianProjection(*hess, true);
		return energy;
	}
	Real TFWStretch::getVertexEnergy(Sheet& sheet, int vid, Vector3* deriv, Matrix33* hess) const
	{
		Real lameAlpha = sheet.getLameAlpha();
		Real lameBeta = sheet.getLameBeta();
		// need to compute deriv and hessian wrt the vertex at index vid
		// just for testing
		Real totalEnergy = 0;
		const Eigen::Vector<int, 32>& faces = sheet.sheetInfo.vertexFacesSet.row(vid);
		unsigned int nfaces = 0;
		for (int i = 0; i < faces.rows(); i++)
		{
			if (faces(i) == -1)
				break;
			nfaces++;
		}
		// create buffers to store parallel results
		if (deriv)
			deriv->setZero();
		if (hess)
			hess->setZero();
		// cache I and II
		for (int i = 0; i < nfaces; i++)
		{
			int faceId = faces(i);
			int localVid = -1;
			for (int ii = 0; ii < 3; ii++)
			{
				int vidii = sheet.getTriangles().at(faceId).v[ii];
				if (vidii == vid)
				{
					localVid = ii;
					break;
				}
			}
			if (localVid == -1)
			{
				std::cout << "Vertex " << localVid << " not found! Face " << faceId << ". Error in tfw_stretching.cpp\n";
				throw "Error!";
			}
			const std::vector<trimesh::triangle_t>& faces = sheet.getTriangles();
			const std::vector<trimesh::edge_t>& edges = sheet.getEdges();
			const Matrix3N& positions = sheet.getPositionsN3();
			const Eigen::MatrixXi& edgeOppositeVertices = sheet.getEdgeOppositeVertices();
			const std::vector<trimesh::triangle_t>& vertexOppositeEdges = sheet.getVertexOppositeEdgesPerFace();

			Matrix33 facePositions;
			int vidxi = faces[faceId].v[0];
			int vidxj = faces[faceId].v[1];
			int vidxk = faces[faceId].v[2];

			facePositions.col(0) = positions.col(faces[faceId].v[0]);
			facePositions.col(1) = positions.col(faces[faceId].v[1]);
			facePositions.col(2) = positions.col(faces[faceId].v[2]);
			Eigen::Matrix<Real, 4, 9> gradI;
			std::array<Eigen::Matrix<Real, 9, 9>, 4> hessI;
			Matrix22 II;
			Matrix22 I = firstFundamentalFormNew(facePositions, true, true, &gradI, &hessI); // correct till here

			Matrix22 Ibar = sheet.getRestShape().getIbars()[faceId];
			Matrix22 Ibarinv = Ibar.inverse();
			Matrix22 IIbar = sheet.getRestShape().getIIbars()[faceId];
			Vector2 dphi = sheet.computeDphi1VectorPerFace(faceId, NULL);
			Eigen::Matrix<Real, 4, 18> gradII; // gradI(i) contains dI(i)/dX
			std::array<Eigen::Matrix<Real, 18, 18>, 4> hessII;
			if (dphi.norm() > 0)
			{
				II.noalias() = sffobject->secondFundamentalFormNew(sheet.getBaseMesh(), sheet.getPositionsN3(), sheet.getEdgeDofs(), faceId, &gradII, &hessII);
			}

			Matrix22 dphidphiTensor;
			if (sheet.cachedDphiDphiTValid())
			{
				dphidphiTensor.noalias() = sheet.getCachedDphiDpiT().second[faceId];
			}
			else
			{
				dphidphiTensor.noalias() = sheet.computeDphiDphi1VectorPerFaceTensorNew(faceId, 0, NULL, NULL);
			}
			Matrix22 dadaTensor;
			if (sheet.cachedDaDaTValid())
			{
				dadaTensor.noalias() = sheet.getCachedDaDaT().second[faceId];
			}
			else
			{
				dadaTensor.noalias() = sheet.computeDaDaTensor(faceId, 0, NULL, NULL);
			}

			Matrix22 dadphiTensor, dphidaTensor;
			if (sheet.cachedDaDphiTValid())
			{
				dadphiTensor.noalias() = sheet.getCachedDaDpiT().second[faceId];
				dphidaTensor.noalias() = dadphiTensor.transpose();
			}
			else
			{
				dadphiTensor.noalias() = sheet.computeDaDphi1VectorPerFaceTensor(faceId, 0, NULL, NULL);
				dphidaTensor.noalias() = dadphiTensor.transpose();
			}
			Matrix22 tempM = (Ibarinv * (dphidaTensor + dadphiTensor));
			Real M2norm = SVNormSquared(sheet.getLameAlpha(), sheet.getLameBeta(), tempM);
			tempM.noalias() = Ibarinv * dadaTensor;
			Real M3norm = SVNormSquared(sheet.getLameAlpha(), sheet.getLameBeta(), tempM) / 8.0;
			Real area = 0.5 * sqrt(Ibar.determinant());
			for (int quadId = 0; quadId < sheet.getQuadPoints().size(); quadId++)
			{

				Vector3 tempDeriv;
				tempDeriv.setZero();
				Matrix33 tempHess;
				tempHess.setZero();
				Real energy = 0;

				Real a = sheet.computeAmplitudesFromQuadNew(faceId, quadId, NULL, NULL, NULL, nullptr, nullptr);

				if (dphi.norm() > 1e-4)
				{
					Vector2 w = Ibarinv * dphi;
					Real wnorm = w.norm();
					Real coeff = 4 * a * a * lameBeta * (lameAlpha + lameBeta) / (lameAlpha + 2 * lameBeta);
					coeff = coeff / (wnorm * wnorm * Ibar.trace() - w.dot(Ibar * w));
					coeff = coeff / (wnorm * wnorm * Ibar.trace() - w.dot(Ibar * w));
					Real W2s = wnorm * wnorm * II.trace() - w.dot(II * w);
					W2s = W2s * W2s;
					energy += coeff * W2s;
					Eigen::Matrix<Real, 2, 2> W2sderiv;
					if (deriv || hess)
						W2sderiv = 2 * coeff * (wnorm * wnorm * Matrix22::Identity() - w * w.transpose()) * (wnorm * wnorm * II.trace() - w.dot(II * w));

					// positions
					if (deriv)
					{
						// Matrix22 dW2sdII = W2sderiv; //.reshaped(2, 2);
						Matrix22 dIIdxi;
						for (int i = 0; i < 3; i++)
						{
							int vidComp = idxAt(localVid, i);
							dIIdxi << gradII(0, vidComp), gradII(1, vidComp), gradII(2, vidComp), gradII(3, vidComp);
							tempDeriv(i) += doubleContraction(W2sderiv, dIIdxi);
						}
					}
					if (hess)
					{
						Matrix22 temp = w.norm() * w.norm() * Matrix22::Identity() - w * w.transpose();
						for (int j = 0; j < 3; j++)
						{
							Matrix22 d2EdxijdII, dIIdxij;
							dIIdxij << mat(gradII, localVid, j);
							d2EdxijdII = 2 * coeff * (temp)*doubleContraction(temp, dIIdxij);
							for (int k = 0; k <= j; k++)
							{
								Matrix22 d2IIdxjdxk, dIIdxik;
								dIIdxik << mat(gradII, localVid, k);
								d2IIdxjdxk << hessMat(hessII, idxAt(localVid, j), idxAt(localVid, k));
								Real reqHess = doubleContraction(W2sderiv, d2IIdxjdxk) + doubleContraction(d2EdxijdII, dIIdxik);
								tempHess(j, k) += reqHess;
								if (j != k)
									tempHess(k, j) += reqHess;
							}
						}
					}
				}

				// First term
				// Correct numerically
				Matrix22 M0 = Ibarinv * (I + 0.5 * a * a * dphidphiTensor + 0.5 * dadaTensor - Ibar);
				Matrix22 dMSVdM;
				energy += SVNormSquared(lameAlpha, lameBeta, M0, &dMSVdM);
				if (deriv)
				{
					// Derivatives wrt face positions
					for (int i = 0; i < 3; i++)
					{
						Matrix22 dIdxij;
						dIdxij << mat(gradI, localVid, i);
						Matrix22 dM0dxij = Ibarinv * dIdxij;
						tempDeriv(i) += doubleContraction(lameAlpha * M0.trace() * Matrix22::Identity() + 2 * lameBeta * M0, dM0dxij);
					}
					// d || M || dxi = d || M || dM : Ibarinv dIdxi
				}
				if (hess)
				{
					for (unsigned int j = 0; j < 3; j++)
					{
						Matrix22 dIdxij;
						dIdxij << mat(gradI, localVid, j);
						// d2MSVdxijdM
						Matrix22 d2MSVdxijdM = lameAlpha * doubleContraction(Ibarinv * dIdxij, Matrix22::Identity()) * Matrix22::Identity() + 2 * lameBeta * Ibarinv * dIdxij;
						for (unsigned int k = 0; k <= j; k++)
						{
							Matrix22 dMdxik, d2Mdxijdxik;
							// dMdxik
							Matrix22 dIdxik;
							dIdxik << mat(gradI, localVid, k);
							dMdxik.noalias() = Ibarinv * dIdxik;
							// d2Mdxijdxik
							Matrix22 hessIjk;
							hessIjk << hessMat(hessI, idxAt(localVid, j), idxAt(localVid, k));
							d2Mdxijdxik.noalias() = Ibarinv * hessIjk;
							Real d2Edxijdxik = doubleContraction(d2MSVdxijdM, dMdxik)
							    + doubleContraction(dMSVdM, d2Mdxijdxik);
							tempHess(j, k) += d2Edxijdxik;
							if (k != j)
								tempHess(k, j) += d2Edxijdxik;
						}
					}
				}
				// Third term
				// Numerically correct
				energy += (a * a / 32.0) * M2norm;

				// Fourth term
				// Correct numerically
				energy += M3norm;

				if (deriv)
					(*deriv) += sheet.getQuadPoints()[quadId].weight * area * sheet.getThickness() / 4.0 * tempDeriv;
				if (hess)
					(*hess) += sheet.getQuadPoints()[quadId].weight * area * sheet.getThickness() / 4.0 * tempHess;
				totalEnergy += sheet.getQuadPoints()[quadId].weight * area * sheet.getThickness() / 4.0 * energy;
			}
		}
		// if (hess)
		// 	*hess = argus::util::hessianProjection(*hess, true);
		return totalEnergy;
	}

	// slower
	// Real TFWStretch::getVertexEnergy(Sheet& sheet, int vid, Vector3* deriv, Matrix33* hess) const
	// {
	// 	// need to compute deriv and hessian wrt the vertex at index vid
	// 	// just for testing
	// 	Real energy = 0;
	// 	const Eigen::VectorXi& faces = sheet.vbdInfo.vertexFacesSet.row(vid);
	// 	unsigned int nfaces = 0;
	// 	for (int i = 0; i < faces.rows(); i++)
	// 	{
	// 		if (faces(i) == -1)
	// 			break;
	// 		nfaces++;
	// 	}
	// 	// create buffers to store parallel results
	// 	std::vector<Real> energies(nfaces, 0);
	// 	std::vector<Vector3> derivatives(nfaces);
	// 	if (deriv)
	// 		deriv->setZero();
	// 	if (hess)
	// 		hess->setZero();
	// 	if (deriv)
	// 		for (int i = 0; i < nfaces; i++)
	// 		{
	// 			derivatives[i].setZero();
	// 		}
	// 	std::vector<Matrix33> hessians(nfaces);
	// 	if (hess)
	// 		for (int i = 0; i < nfaces; i++)
	// 		{
	// 			hessians[i].setZero();
	// 		}
	// 	// #ifdef USE_OMP
	// 	// #pragma omp parallel for
	// 	// #endif
	// 	for (int i = 0; i < nfaces; i++)
	// 	{
	// 		int faceId = faces(i);
	// 		int localVid = -1;
	// 		for (int ii = 0; ii < 3; ii++)
	// 		{
	// 			int vidii = sheet.getTriangles().at(faceId).v[ii];
	// 			if (vidii == vid)
	// 			{
	// 				localVid = ii;
	// 				break;
	// 			}
	// 		}
	// 		if (localVid == -1)
	// 		{
	// 			std::cout << "Vertex " << localVid << " not found! Face " << faceId << ". Error in tfw_stretching.cpp\n";
	// 			exit(0);
	// 		}
	// 		for (int quad = 0; quad < sheet.getQuadPoints().size(); quad++)
	// 		{
	// 			Vector3 tempDeriv1, tempDeriv2, tempDeriv;
	// 			Matrix33 tempHess, tempHess1, tempHess2;
	// 			/* 	// test hessians
	// 			    Vector3 dx = Vector3::Random().normalized();
	// 			    Real h = 0.1;
	// 			    Vector3 X = sheet.getPositionsN3().col(vid);
	// 			    TFWStretch::getEnergyDensityFromQuadPerVtex(localVid, faceId, quad, sheet, &tempDeriv, &tempHess);
	// 			    std::cout << "Temphess\n"
	// 			              << tempHess << std::endl;
	// 			    for (int ii = 0; ii < 20; ii++)
	// 			    {
	// 			        Vector3 hdx = X + 0.5 * h * dx;
	// 			        sheet.getPositionsN3().col(vid) = hdx;
	// 			        Real e1 = TFWStretch::getEnergyDensityFromQuadPerVtex(localVid, faceId, quad, sheet, &tempDeriv1, &tempHess1);
	// 			        hdx = X - 0.5 * h * dx;
	// 			        sheet.getPositionsN3().col(vid) = hdx;
	// 			        Real e2 = TFWStretch::getEnergyDensityFromQuadPerVtex(localVid, faceId, quad, sheet, &tempDeriv2, &tempHess2);
	// 			        // std::cout << (tempDeriv1 - tempDeriv2).transpose() / h << " - " << (tempHess * dx).transpose() << std::endl;
	// 			        std::cout << "Error Hessians " << (((tempDeriv1 - tempDeriv2) / h) - tempHess * dx).norm() << std::endl;
	// 			        // std::cout << "Error derivative" << std::fabs(((e1 - e2) / h) - dx.dot(tempDeriv)) << std::endl;
	// 			        h /= 10.0;
	// 			    }
	// 			    exit(0); */
	// 			energies[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * TFWStretch::getEnergyDensityFromQuadPerVtex(localVid, faceId, quad, sheet, deriv ? &tempDeriv : NULL, hess ? &tempHess : NULL) / 4.0;
	// 			if (deriv)
	// 				derivatives[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * tempDeriv / 4.0;
	// 			if (hess)
	// 				hessians[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * tempHess / 4.0;
	// 		}
	// 	}
	//
	// 	// consolidate
	// 	for (int i = 0; i < nfaces; i++)
	// 	{
	// 		energy += energies[i];
	// 		if (deriv)
	// 			(*deriv) += derivatives[i];
	// 		if (hess)
	// 			(*hess) += hessians[i];
	// 	}
	// 	return energy;
	// }

	// Function for entire sheet
	Real TFWStretch::getStretchingEnergyDensityFromQuadAmpsPerVtexDphis1VectorPerFace(int faceId, int quadId, Sheet& sheet, VectorN* deriv, MatrixNN* hess) const
	{
		// TODO: Delete
		bool _debug = false;
		const std::vector<trimesh::triangle_t>& faces = sheet.getTriangles();
		const std::vector<trimesh::edge_t>& edges = sheet.getEdges();
		const Matrix3N& positions = sheet.getPositionsN3();
		const Eigen::MatrixXi& edgeOppositeVertices = sheet.getEdgeOppositeVertices();
		const std::vector<trimesh::triangle_t>& vertexOppositeEdges = sheet.getVertexOppositeEdgesPerFace();
		int nfaces = faces.size();
		int nedges = edges.size();
		int nvertices = positions.cols();
		Matrix33 facePositions;
		int vidxi = faces[faceId].v[0];
		int vidxj = faces[faceId].v[1];
		int vidxk = faces[faceId].v[2];
		// facePositions.col(0) = positions.block<3, 1>(0, faces[faceId].v[0]);
		// facePositions.col(1) = positions.block<3, 1>(0, faces[faceId].v[1]);
		// facePositions.col(2) = positions.block<3, 1>(0, faces[faceId].v[2]);
		facePositions.col(0) = positions.col(faces[faceId].v[0]);
		facePositions.col(1) = positions.col(faces[faceId].v[1]);
		facePositions.col(2) = positions.col(faces[faceId].v[2]);
		MatrixNN gradI, gradII; // gradI(i) contains dI(i)/dX
		std::vector<Eigen::Matrix<Real, 9, 9>> hessI;
		Matrix22 I;
		I = firstFundamentalForm(facePositions, hess || deriv ? true : false, hess ? true : false, hess || deriv ? &gradI : NULL, hess ? &hessI : NULL); // correct till here
		std::vector<Eigen::MatrixXd> hessII;
		std::shared_ptr<SecondFundamentalFormDiscretization> sff;
		sff = std::make_shared<MidedgeAverageFormulation>();
		Matrix22 II = sff->secondFundamentalForm(sheet.getBaseMesh(), sheet.getPositionsN3().transpose(), sheet.getEdgeDofs(), faceId, hess || deriv ? &gradII : NULL, hess ? &hessII : NULL);
		Matrix22 Ibar = sheet.getRestShape().getIbars()[faceId];
		Matrix22 Ibarinv = Ibar.inverse();
		Matrix22 IIbar = sheet.getRestShape().getIIbars()[faceId];
		Vector3 gradA; // datilda/da
		Vector2 da;
		Matrix23 gradDA; // dda/da
		Real u = sheet.getQuadPoints()[quadId].u;
		Real v = sheet.getQuadPoints()[quadId].v;
		Vector3 ampUVWeights;
		ampUVWeights << 1 - u - v, u, v;
		Real a = sheet.computeAmplitudesFromQuad(faceId, quadId, &da, hess || deriv ? &gradA : NULL, hess || deriv ? &gradDA : NULL, nullptr, nullptr);
		Matrix22 gradDphi;
		Vector2 dphi = sheet.computeDphi1VectorPerFace(faceId, deriv ? &gradDphi : NULL);
		std::vector<Matrix22> gradDphiDphiTensor;
		std::vector<Matrix22> hessDphiDphiT;
		Matrix22 dphidphiTensor;
		dphidphiTensor = sheet.computeDphiDphi1VectorPerFaceTensor(faceId, quadId, hess || deriv ? &gradDphiDphiTensor : NULL, hess ? &hessDphiDphiT : NULL);
		std::vector<Matrix22> gradDaDaTensor;
		std::vector<Matrix33> hessDaDa;
		Matrix22 dadaTensor;
		dadaTensor = sheet.computeDaDaTensor(faceId, quadId, hess || deriv ? &gradDaDaTensor : NULL, hess ? &hessDaDa : NULL);
		std::vector<Matrix22> gradDaDphiTensor;
		std::vector<Eigen::Matrix<Real, 5, 5>> hessDaDphiTensor, hessDphiDaTensor;
		Matrix22 dadphiTensor = sheet.computeDaDphi1VectorPerFaceTensor(faceId, quadId, deriv || hess ? &gradDaDphiTensor : NULL, hess ? &hessDaDphiTensor : NULL);
		std::vector<Matrix22> gradDphiDaTensor;
		Matrix22 dphidaTensor = sheet.computeDphi1VectorPerFaceDaTensor(faceId, quadId, deriv || hess ? &gradDphiDaTensor : NULL, hess ? &hessDphiDaTensor : NULL);
		Real lameAlpha = sheet.getLameAlpha();
		Real lameBeta = sheet.getLameBeta();
		if (_debug && deriv && faceId == 8)
		{
			std::cout << "getEnergy() quad " << quadId << std::endl;
			// print debinfo
			std::cout << "a " << a << std::endl;
			std::cout << "Dphi " << dphi.transpose() << std::endl;
			std::cout << "I \n"
			          << I << std::endl;
			std::cout << "Ibar \n"
			          << Ibar << std::endl;
			std::cout << "II \n"
			          << II << std::endl;
		}
		Real energy = 0.0;

		if (deriv)
		{
			// Per face
			deriv->resize(3 + 2 + 3 * 3 + 3 * 3); // last 9 for the adjacent face vertices
			deriv->setZero();
		}

		if (hess)
		{
			hess->resize(23, 23);
			hess->setZero();
		}

		// Second term
		// Derivatives are correct numerically
		if (dphi.norm() > 1e-4)
		{
#ifdef ARGUS_USE_FLOATS
			typedef TinyAD::Float<9> AReal;
#else
			typedef TinyAD::Double<9> AReal;
#endif
			auto W2sFunc = [u, v, &Ibarinv, &Ibar, &lameAlpha, &lameBeta](Eigen::Vector<Real, 9> X, Eigen::Vector<Real, 9>* deriv = NULL, Eigen::Matrix<Real, 9, 9>* hess = NULL)
			{
				Eigen::Vector<AReal, 9> activeX = AReal::make_active(X);
				AReal a = (1 - u - v) * activeX(0) + u * activeX(1) + v * activeX(2);
				Eigen::Vector<AReal, 2> w = Ibarinv * activeX(Eigen::seq(3, 4));
				Real coeff = 4 * lameBeta * (lameAlpha + lameBeta) / (lameAlpha + 2 * lameBeta);
				Eigen::Matrix<AReal, 2, 2> II = activeX(Eigen::seq(5, 8)).reshaped(2, 2);
				AReal num = w.norm() * w.norm() * II.trace() - w.transpose() * II * w;
				AReal den = w.norm() * w.norm() * Ibar.trace() - w.transpose() * Ibar * w;
				AReal e = coeff * a * a * (num / den) * (num / den);
				if (deriv)
					(*deriv) = e.grad;
				if (hess)
					(*hess) = e.Hess;
				return TinyAD::to_passive(e);
			};
			Eigen::Vector<Real, 9> W2sX;
			W2sX(0) = sheet.getAmplitudes()(vidxi);
			W2sX(1) = sheet.getAmplitudes()(vidxj);
			W2sX(2) = sheet.getAmplitudes()(vidxk);
			W2sX(3) = dphi(0);
			W2sX(4) = dphi(1);
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					W2sX(i * 2 + j + 5) = II(i, j);
				}
			}
			Eigen::Vector<Real, 9> W2sderiv;
			Eigen::Matrix<Real, 9, 9> W2shess;
			energy += W2sFunc(W2sX, deriv ? &W2sderiv : NULL, hess ? &W2shess : NULL);
			if (deriv)
			{
				// amps and dphis
				for (int i = 0; i < 5; i++)
				{
					deriv->coeffRef(i) += W2sderiv(i);
				}
				// positions
				Matrix22 dW2sdII = W2sderiv(Eigen::seq(5, 8)).reshaped(2, 2);
				for (int i = 0; i < 18; i++)
				{
					Matrix22 dIIdxi;
					dIIdxi << gradII(0, i), gradII(1, i), gradII(2, i), gradII(3, i);
					deriv->coeffRef(5 + i) += doubleContraction(dW2sdII, dIIdxi);
				}
			}
			if (hess)
			{
				// amps, dphi vs amps, dphis
				hess->block(0, 0, 5, 5) += W2shess.block(0, 0, 5, 5);

				Vector2 w = Ibarinv * dphi;
				Real num = w.norm() * w.norm() * II.trace() - w.transpose() * II * w;
				Real den = w.norm() * w.norm() * Ibar.trace() - w.transpose() * Ibar * w;
				// amps vs pos
				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 18; j++)
					{
						Matrix22 dIIdxj;
						dIIdxj << gradII(0, j), gradII(1, j), gradII(2, j), gradII(3, j);
						Matrix22 d2EdaidII = 16 * a * gradA(i) * lameBeta * (lameBeta + lameAlpha)
						    * num * (w.norm() * w.norm() * Matrix22::Identity() - w * w.transpose())
						    / ((lameAlpha + 2 * lameBeta) * (den * den));
						hess->coeffRef(i, 5 + j) += doubleContraction(d2EdaidII, dIIdxj);
						hess->coeffRef(5 + j, i) += doubleContraction(d2EdaidII, dIIdxj);
					}
				}
				// dphis vs pos
				// Ignoring for now, assuming dphis do not change the energy too much

				// pos vs pos
				for (int i = 0; i < 18; i++)
				{
					for (int j = 0; j < 18; j++)
					{
						Matrix22 gradIIi;
						gradIIi << gradII(0, i), gradII(1, i), gradII(2, i), gradII(3, i);
						Matrix22 hessIIij;
						hessIIij << hessII.at(0)(i, j), hessII.at(1)(i, j), hessII.at(2)(i, j), hessII.at(3)(i, j);
						Real coeff11 = 2 * 4 * a * a * lameBeta * (lameAlpha + lameBeta) / (lameAlpha + 2 * lameBeta);
						coeff11 = coeff11 / (den * den);
						Matrix22 coeff1 = coeff11 * (w.norm() * w.norm() * Matrix22::Identity() - w * w.transpose());
						Matrix22 T1, T2, T3, T4;
						T1 = coeff1 * doubleContraction((w.norm() * w.norm() * Matrix22::Identity() - w * w.transpose()), gradIIi);
						T2 = gradIIi;
						T3 = coeff1 * num;
						T4 = hessIIij;
						hess->coeffRef(5 + i, 5 + j) += doubleContraction(T1, T2) + doubleContraction(T3, T4);
					}
				}
			}
		}
		// First term
		// Correct numerically
		Matrix22 M0 = Ibarinv * (I + 0.5 * a * a * dphidphiTensor + 0.5 * dadaTensor - Ibar);
		if (_debug && deriv && faceId == 8)
		{
			std::cout << "M1\n"
			          << M0 << std::endl;
		}
		energy += SVNormSquared(lameAlpha, lameBeta, M0);
		if (deriv)
		{
			// Derivatives wrt face amplitudes
			for (int i = 0; i < 3; i++)
			{
				Matrix22 dM0dai = Ibarinv * (a * gradA(i) * dphidphiTensor + 0.5 * gradDaDaTensor[i]);
				deriv->coeffRef(i) += 2 * (dM0dai * M0).trace() * lameBeta + dM0dai.trace() * M0.trace() * lameAlpha;
			}
			// // Derivatives wrt face dphis
			for (int i = 0; i < 2; i++)
			{
				Matrix22 dM0dDphii = Ibarinv * (0.5 * a * a * gradDphiDphiTensor[i]);
				deriv->coeffRef(3 + i) += 2 * (dM0dDphii * M0).trace() * lameBeta + dM0dDphii.trace() * M0.trace() * lameAlpha;
			}
			// Derivatives wrt face positions
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					Matrix22 dIdxij;
					dIdxij << gradI(0, 3 * i + j), gradI(1, 3 * i + j), gradI(2, 3 * i + j), gradI(3, 3 * i + j);
					Matrix22 dM0dxij = Ibarinv * dIdxij;
					deriv->coeffRef(5 + 3 * i + j) += 2 * (dM0dxij * M0).trace() * lameBeta + dM0dxij.trace() * M0.trace() * lameAlpha;
				}
			}
		}

		if (hess)
		{
			// amps vs amps
			Matrix22 T1, T2, T3, T4;
			T2 = lameAlpha * M0.trace() * Matrix22::Identity() + 2 * lameBeta * M0;
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					T1 << hessDaDa.at(0)(i, j), hessDaDa.at(1)(i, j),
					    hessDaDa.at(2)(i, j), hessDaDa.at(3)(i, j);
					T1 = (Ibarinv * T1 + 2 * gradA(i) * gradA(j) * dphidphiTensor).eval();
					T3 = 0.5 * Ibarinv * (gradDaDaTensor.at(i) + 2 * a * gradA(i) * dphidphiTensor);
					Matrix22 tempT4;
					// T4 = Ibarinv * 0.5 * (2*lameBeta * (T4 + 2 * a * gradA(j) * dphidphiTensor) + lameAlpha *( gradDaDaTensor.at(j)(0,0) + gradDaDaTensor.at(j)(1,1) + 2 *a * gradA(j) * (dphidphiTensor(0,0) + dphidphiTensor(1,1)))*Matrix22::Identity()).eval();
					T4 = lameBeta * Ibarinv * (gradDaDaTensor.at(j) + 2 * a * gradA(j) * dphidphiTensor);
					T4 = (T4 + 0.5 * lameAlpha * (doubleContraction(Ibarinv, gradDaDaTensor.at(j) + 2 * a * gradA(j) * dphidphiTensor)) * Matrix22::Identity()).eval();
					hess->coeffRef(i, j) += d2Edqiqj(T1, T2, T3, T4);
				}
			}
			// dphis vs amps
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					T1 = Ibarinv * a * gradA[j] * gradDphiDphiTensor.at(i);
					T3 = gradDaDaTensor.at(j); //<< gradDaDaTensor.at(0)[j], gradDaDaTensor.at(1)[j],
					                           // gradDaDaTensor.at(2)[j], gradDaDaTensor.at(3)[j];
					T3 = (0.5 * Ibarinv * (T3 + 2 * a * gradA(j) * dphidphiTensor)).eval();
					T4 = (lameBeta * Ibarinv * a * a * gradDphiDphiTensor.at(i) + 0.5 * lameAlpha * doubleContraction(Ibarinv, a * a * gradDphiDphiTensor.at(i)) * Matrix22::Identity()).eval();
					hess->coeffRef(3 + i, j) += d2Edqiqj(T1, T2, T3, T4);
					hess->coeffRef(j, 3 + i) += d2Edqiqj(T1, T2, T3, T4);
				}
			}
			// dphis vs dphis
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					T1 << hessDphiDphiT.at(0)(i, j), hessDphiDphiT.at(1)(i, j),
					    hessDphiDphiT.at(2)(i, j), hessDphiDphiT.at(3)(i, j);
					T1 = (Ibarinv * 0.5 * a * a * T1).eval();
					T3 = 0.5 * Ibarinv * a * a * gradDphiDphiTensor.at(j);
					T4 = a * a * (lameBeta * Ibarinv * gradDphiDphiTensor.at(i) + 0.5 * lameAlpha * doubleContraction(Ibarinv, gradDphiDphiTensor.at(j).transpose()) * Matrix22::Identity()).eval();
					hess->coeffRef(3 + i, 3 + j) += d2Edqiqj(T1, T2, T3, T4);
				}
			}
			// pos vs pos
			for (int i = 0; i < 9; i++)
			{
				for (int j = 0; j < 9; j++)
				{
					T1 << hessI.at(0)(i, j), hessI.at(1)(i, j),
					    hessI.at(2)(i, j), hessI.at(3)(i, j);
					T1 = (Ibarinv * T1).eval();
					T3 << gradI(0, j), gradI(1, j),
					    gradI(2, j), gradI(3, j);
					T3 = (Ibarinv * T2).eval();
					Matrix22 matGradIi;
					matGradIi << gradI(0, i), gradI(1, i),
					    gradI(2, i), gradI(3, i);
					T4 = lameAlpha * doubleContraction(Ibarinv, matGradIi) * Matrix22::Identity() + 2 * lameBeta * Ibarinv * matGradIi;
					hess->coeffRef(5 + i, 5 + j) += d2Edqiqj(T1, T2, T3, T4);
				}
			}
		}
		// // Third term
		// Numerically correct
		Matrix22 M2 = (Ibarinv * (dphidaTensor + dadphiTensor));
		energy += (a * a / 32.0) * SVNormSquared(sheet.getLameAlpha(), sheet.getLameBeta(), M2);
		if (deriv)
		{
			// wrt a*a part
			for (int i = 0; i < 3; i++)
			{
				// This is correct as a*a/32 is not in M2, so the final product will be correct
				Matrix22 dM2dai = ((a * gradA(i) / 32.0) * Ibarinv * (dphidaTensor + dadphiTensor)) + (a * a / 32.0) * Ibarinv * (gradDaDphiTensor[i] + gradDphiDaTensor[i]);
				deriv->coeffRef(i) += 2 * ((dM2dai * M2).trace() * lameBeta) + dM2dai.trace() * M2.trace() * lameAlpha;
			}
			for (int i = 3; i < 5; i++)
			{
				Matrix22 dM2daphi = (a * a / 32.0) * Ibarinv * (gradDaDphiTensor[i] + gradDphiDaTensor[i]);
				deriv->coeffRef(i) += 2 * (dM2daphi * M2).trace() * lameBeta + dM2daphi.trace() * M2.trace() * lameAlpha;
			}
		}
		if (hess)
		{
			Matrix22 T, T1, T2, T3, T4, gradM2i, gradM2j;
			T = M2 * a;
			T2 = lameAlpha * (T.trace()) * Matrix22::Identity() + lameBeta * T * 2;
			// amps vs amps
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					gradM2i = Ibarinv * (gradDphiDaTensor.at(i) + gradDaDphiTensor.at(i));
					gradM2j = Ibarinv * (gradDphiDaTensor.at(j) + gradDaDphiTensor.at(j));
					T1 = gradA(j) * gradM2i + gradA(i) * gradM2j;
					T3 = gradA(j) * M2 + a * gradM2j;
					T4 = lameAlpha * doubleContraction(Ibarinv, gradA(i) * M2 + a * gradM2i) * Matrix22::Identity()
					    + lameBeta * 2 * (gradA(i) * M2 + a * gradM2i);
					hess->coeffRef(i, j) += d2Edqiqj(T1, T2, T3, T4) / 32.0;
				}
			}
			// dphis vs amps
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					gradM2j = Ibarinv * (gradDphiDaTensor.at(j) + gradDaDphiTensor.at(j));
					Matrix22 hessDaDphiTij, hessDphiDaTij;
					hessDaDphiTij << hessDaDphiTensor.at(0)(3 + i, j), hessDaDphiTensor.at(1)(3 + i, j),
					    hessDaDphiTensor.at(2)(3 + i, j), hessDaDphiTensor.at(3)(3 + i, j);
					hessDphiDaTij << hessDphiDaTensor.at(0)(3 + i, j), hessDphiDaTensor.at(1)(3 + i, j),
					    hessDphiDaTensor.at(2)(3 + i, j), hessDphiDaTensor.at(3)(3 + i, j);
					T1 = gradA(j) * Ibarinv * (gradDaDphiTensor.at(3 + i) + gradDphiDaTensor.at(3 + i))
					    + a * Ibarinv * (hessDaDphiTij + hessDphiDaTij);
					T3 = gradA(j) * M2 + a * gradM2j;
					T4 = lameAlpha * a * doubleContraction(Ibarinv, gradDaDphiTensor.at(3 + i).transpose() + gradDphiDaTensor.at(3 + i).transpose()) * Matrix22::Identity()
					    + lameBeta * 2 * a * Ibarinv * (gradDaDphiTensor.at(3 + i) + gradDphiDaTensor.at(3 + i));
					hess->coeffRef(3 + i, j) += d2Edqiqj(T1, T2, T3, T4) / 32.0;
					hess->coeffRef(j, 3 + i) += d2Edqiqj(T1, T2, T3, T4) / 32.0;
				}
			}
			// dphis vs dphis
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					T1 = Matrix22::Zero();
					gradM2j = gradDphiDaTensor.at(3 + j) + gradDaDphiTensor.at(3 + j);
					gradM2i = gradDphiDaTensor.at(3 + i) + gradDaDphiTensor.at(3 + i);
					T3 = a * Ibarinv * gradM2j;
					T4 = lameAlpha * a * doubleContraction(Ibarinv, gradM2i) * Matrix22::Identity() + 2 * lameBeta * a * Ibarinv * gradM2i;
					hess->coeffRef(3 + i, 3 + j) += d2Edqiqj(T1, T2, T3, T4) / 32.0;
				}
			}
		}
		// // Fourth term
		// Correct numerically
		Matrix22 M3 = Ibarinv * dadaTensor;
		energy += (1.0 / 8.0) * SVNormSquared(sheet.getLameAlpha(), sheet.getLameBeta(), M3);
		if (deriv)
		{
			Matrix22 dM3da;
			for (int i = 0; i < 3; i++)
			{
				dM3da = (1.0 / 8.0) * Ibarinv * (gradDaDaTensor[i]);
				deriv->coeffRef(i) += 2 * (dM3da * M3).trace() * lameBeta + dM3da.trace() * M3.trace() * lameAlpha;
			}
		}
		if (hess)
		{
			M3 = (M3 / 8.0).eval();
			Matrix22 T1, T2, T3, T4;
			T2 = lameAlpha * M3.trace() * Matrix22::Identity() + 2 * lameBeta * M3;
			// amps vs amps
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					T1 << hessDaDa.at(0)(i, j), hessDaDa.at(1)(i, j),
					    hessDaDa.at(2)(i, j), hessDaDa.at(3)(i, j);
					T1 = ((1 / 8) * Ibarinv * T1).eval();
					T3 = (1 / 8) * Ibarinv * gradDaDaTensor.at(j);
					T4 = lameAlpha * doubleContraction(Ibarinv, gradDaDaTensor.at(i)) * Matrix22::Identity() + 2 * lameBeta * gradDaDaTensor.at(i);
					hess->coeffRef(i, j) += d2Edqiqj(T1, T2, T3, T4);
				}
			}
		}

		Real area = 0.5 * sqrt(Ibar.determinant());
		if (deriv)
			(*deriv) = area * (*deriv);
		if (hess)
		{
			(*hess) = area * (*hess);
			// project
			(*hess) = argus::util::lowRankApprox(*hess);
		}
		return area * energy;
	}
	Real TFWStretch::getStretchingEnergyAmpsPerVtexDphis1VectorPerFace(Sheet& sheet, VectorN* deriv, std::vector<Eigen::Triplet<Real>>* hess) const
	{
		bool _debug = false;
		// iterate through all the faces in the mesh
		// iterate through each quad point and sum up the energies
		// We can parallelize energy computation for each face, thus, we want to create an array per face for energy and derivatives
		const std::vector<trimesh::triangle_t>& edgeIDOppVertexIDPerface = sheet.getVertexOppositeEdgesPerFace();
		const std::vector<trimesh::triangle_t>& faces = sheet.getTriangles();
		const Eigen::MatrixXi& vtexIDOppositeEdge = sheet.getEdgeOppositeVertices();
		int nfaces = faces.size();
		int nvertices = sheet.getPositionsN3().cols();
		int nedges = sheet.getEdgeCount();
		std::vector<Real> energies(nfaces, 0.0);
		std::vector<VectorN> derivatives;
		std::vector<MatrixNN> hessians;

		if (deriv)
		{
			derivatives.resize(nfaces);
			for (int i = 0; i < nfaces; i++)
			{
				derivatives[i].resize(3 + 2 + 3 * 3 + 3 * 3);
				derivatives[i].setZero();
			}
		}

		if (hess)
		{
			hessians.resize(nfaces);
			for (int i = 0; i < nfaces; i++)
			{
				hessians[i].resize(23, 23);
				hessians[i].setZero();
			}
		}
		// #ifdef USE_OMP
		// #pragma omp parallel for
		// #endif
		for (int i = 0; i < nfaces; i++)
		{
			// if(tensionFaces.find(i)!=tensionFaces.end())
			// continue;
			// for each face iterate through each quad to get the energies
			for (int quad = 0; quad < sheet.getQuadPoints().size(); quad++)
			{
				if (hess && deriv)
				{
					VectorN tempDeriv;
					MatrixNN tempHess;
					energies[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * TFWStretch::getStretchingEnergyDensityFromQuadAmpsPerVtexDphis1VectorPerFace(i, quad, sheet, &tempDeriv, &tempHess) / 4.0;
					derivatives[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * tempDeriv / 4.0;
					hessians[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * tempHess / 4.0;
					// std::cout <<i << " \n" <<derivatives[i] << std::endl;
				}
				else if (deriv)
				{
					VectorN tempDeriv;

					energies[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * TFWStretch::getStretchingEnergyDensityFromQuadAmpsPerVtexDphis1VectorPerFace(i, quad, sheet, &tempDeriv, nullptr) / 4.0;
					derivatives[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * tempDeriv / 4.0;
					// std::cout <<i << " \n" <<derivatives[i] << std::endl;
				}
				else
				{
					energies[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * TFWStretch::getStretchingEnergyDensityFromQuadAmpsPerVtexDphis1VectorPerFace(i, quad, sheet, nullptr, nullptr) / 4.0;
				}
			}
		}
		if (_debug && deriv)
		{
			std::cout << "getEnergy() Gradient for face 8\n";
			std::cout << derivatives[8] << std::endl;
			exit(0);
		}
		Real energy = 0.0;
		for (int i = 0; i < nfaces; i++)
			energy += energies[i];
		// consolidate derivatives
		if (deriv)
		{
			deriv->resize(nvertices + 2 * nfaces + 3 * nvertices);
			deriv->setZero();
			for (int i = 0; i < nfaces; i++)
			{
				deriv->coeffRef(nvertices + 2 * i) += derivatives[i][3];
				deriv->coeffRef(nvertices + 2 * i + 1) += derivatives[i][4];
				for (int j = 0; j < 3; j++)
				{
					int vid = faces[i].v[j];
					int edgeId = edgeIDOppVertexIDPerface[i].v[j];
					int vidTilda = vtexIDOppositeEdge(edgeId, 0) == vid ? vtexIDOppositeEdge(edgeId, 1) : vtexIDOppositeEdge(edgeId, 0);
					deriv->coeffRef(vid) += derivatives[i][j];
					// filling position derivatives
					for (int vc = 0; vc < 3; vc++)
					{
						deriv->coeffRef(nvertices + 2 * nfaces + 3 * vid + vc) += derivatives[i][5 + 3 * j + vc];
						if (vidTilda != -1)
							deriv->coeffRef(nvertices + 2 * nfaces + 3 * vidTilda + vc) += derivatives[i][14 + 3 * j + vc];
					}
				}
			}
		}
		if (hess)
		{
			// hess is a matrix of dimensions 23x23
			// right now we only need the position hessians for vbd
			// clear hess as it might contain triplets from other energies that will linger
			hess->clear();
			int idxj, idxk;
			for (int i = 0; i < nfaces; i++)
			{
				for (int j = 0; j < 3; j++) // hessians involving amps
				{
					idxj = faces[i].v[j];
					for (int k = 0; k < 3; k++)
					{
						idxk = faces[i].v[k];
						hess->push_back(Eigen::Triplet<Real>(idxj, idxk, hessians[i](j, k)));
					}

					// amps vs dphis and viceversa
					for (int k = 0; k < 2; k++)
					{
						idxk = nvertices + 2 * i + k;
						hess->push_back(Eigen::Triplet<Real>(idxj, idxk, hessians[i](j, 3 + k)));
						hess->push_back(Eigen::Triplet<Real>(idxk, idxj, hessians[i](j, 3 + k)));
					}

					// amps vs positions and viceversa
					for (int k = 0; k < 3; k++)
					{
						int vid = faces[i].v[k];
						for (int l = 0; l < 3; l++)
						{
							idxk = nvertices + 2 * nfaces + 3 * vid + l;
							hess->push_back(Eigen::Triplet<Real>(idxj, idxk, hessians[i](j, 5 + 3 * k + l)));
							hess->push_back(Eigen::Triplet<Real>(idxk, idxj, hessians[i](j, 5 + 3 * k + l)));
						}
						// oppsite vertex
						int edgeId = edgeIDOppVertexIDPerface[i].v[j];
						vid = vtexIDOppositeEdge(edgeId, 0) == vid ? vtexIDOppositeEdge(edgeId, 1) : vtexIDOppositeEdge(edgeId, 0);
						if (vid != -1)
						{
							for (int l = 0; l < 3; l++)
							{
								idxk = nvertices + 2 * nfaces + 3 * vid + l;
								hess->push_back(Eigen::Triplet<Real>(idxj, idxk, hessians[i](j, 14 + 3 * k + l)));
								hess->push_back(Eigen::Triplet<Real>(idxk, idxj, hessians[i](j, 14 + 3 * k + l)));
							}
						}
					}
				}

				for (int j = 0; j < 2; j++) // hessians involving dphis
				{
					idxj = nvertices + 2 * i + j;
					// dphis vs dphis
					for (int k = 0; k < 2; k++)
					{
						idxk = nvertices + 2 * i + k;
						hess->push_back(Eigen::Triplet<Real>(idxj, idxk, hessians[i](3 + j, 3 + k)));
					}
					// dphis vs pos
					for (int k = 0; k < 3; k++)
					{
						int vid = faces[i].v[k];
						for (int l = 0; l < 3; l++)
						{
							idxk = nvertices + 2 * nfaces + 3 * vid + l;
							hess->push_back(Eigen::Triplet<Real>(idxj, idxk, hessians[i](3 + j, 5 + 3 * k + l)));
							hess->push_back(Eigen::Triplet<Real>(idxk, idxj, hessians[i](3 + j, 5 + 3 * k + l)));
						}
						// oppsite vertex
						int edgeId = edgeIDOppVertexIDPerface[i].v[j];
						vid = vtexIDOppositeEdge(edgeId, 0) == vid ? vtexIDOppositeEdge(edgeId, 1) : vtexIDOppositeEdge(edgeId, 0);

						if (vid != -1)
						{
							for (int l = 0; l < 3; l++)
							{
								idxk = nvertices + 2 * nfaces + 3 * vid + l;
								hess->push_back(Eigen::Triplet<Real>(idxj, idxk, hessians[i](j, 14 + 3 * k + l)));
								hess->push_back(Eigen::Triplet<Real>(idxk, idxj, hessians[i](j, 14 + 3 * k + l)));
							}
						}
					}
				}
				// it has to be 18 x 18
				// pos vs pos
				for (int j = 0; j < 3; j++) // hessians involving positions
				{
					int vidj = faces[i].v[j];
					for (int jj = 0; jj < 3; jj++)
					{
						idxj = nvertices + 2 * nfaces + 3 * vidj + jj;
						for (int k = 0; k < 3; k++)
						{
							int vidk = faces[i].v[k];
							for (int kk = 0; kk < 3; kk++)
							{
								idxk = nvertices + 2 * nfaces + 3 * vidk + kk;
								hess->push_back(Eigen::Triplet<Real>(idxj, idxk, hessians[i](5 + 3 * j + jj, 5 + 3 * k + kk)));
							}
						}

						// oppsite vertex
						int edgeId = edgeIDOppVertexIDPerface[i].v[j];
						int vidk = vtexIDOppositeEdge(edgeId, 0) == vidj ? vtexIDOppositeEdge(edgeId, 1) : vtexIDOppositeEdge(edgeId, 0);

						if (vidk != -1)
						{
							for (int k = 0; k < 3; k++)
							{
								for (int kk = 0; kk < 3; kk++)
								{
									idxk = nvertices + 2 * nfaces + 3 * vidk + kk;
									hess->push_back(Eigen::Triplet<Real>(idxj, idxk, hessians[i](5 + 3 * j + jj, 14 + 3 * k + kk)));
									hess->push_back(Eigen::Triplet<Real>(idxk, idxj, hessians[i](5 + 3 * j + jj, 14 + 3 * k + kk)));
								}
							}
						}
					}
				}
			}
		}
		return energy;
	}

	void TFWStretch::getEnergy(Sheet& sheet, Real& energy, VectorN* deriv, std::vector<Eigen::Triplet<Real>>* hess) const
	{
		if (hess)
			energy = TFWStretch::getStretchingEnergyAmpsPerVtexDphis1VectorPerFace(sheet, deriv, hess);
		else if (deriv)
			energy = TFWStretch::getStretchingEnergyAmpsPerVtexDphis1VectorPerFace(sheet, deriv, NULL);
		else
			energy = TFWStretch::getStretchingEnergyAmpsPerVtexDphis1VectorPerFace(sheet, NULL, NULL);
	}

	void TFWStretch::getWrinkleShellEnergy(Sheet& sheet, Real& energy, VectorN* deriv, std::vector<Eigen::Triplet<Real>>* hess) const
	{
		if (hess)
			energy = TFWStretch::getStretchingEnergyAmpsPerVtexDphis1VectorPerFace(sheet, deriv, hess);
		else if (deriv)
			energy = TFWStretch::getStretchingEnergyAmpsPerVtexDphis1VectorPerFace(sheet, deriv, NULL);
		else
			energy = TFWStretch::getStretchingEnergyAmpsPerVtexDphis1VectorPerFace(sheet, NULL, NULL);
	}

	void TFWStretch::getWrinkleShellEnergyPerFace(Sheet& sheet, Real& energy, int faceId, VectorN* deriv, MatrixNN* hess) const
	{
		energy = 0;
		VectorN derivatives;
		if (deriv)
		{
			derivatives.resize(3 + 2 + 3 * 3 + 3 * 3);
			derivatives.setZero();
		}

		for (int quad = 0; quad < sheet.getQuadPoints().size(); quad++)
		{
			if (deriv)
			{
				VectorN tempDeriv;

				energy += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * TFWStretch::getStretchingEnergyDensityFromQuadAmpsPerVtexDphis1VectorPerFace(faceId, quad, sheet, &tempDeriv, nullptr) / 4.0;
				derivatives += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * tempDeriv / 4.0;
				// std::cout <<i << " \n" <<derivatives[i] << std::endl;
			}
			else
			{
				energy += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * TFWStretch::getStretchingEnergyDensityFromQuadAmpsPerVtexDphis1VectorPerFace(faceId, quad, sheet, nullptr, nullptr) / 4.0;
			}
		}
		// assuming that this function is used only for dphi per face optimization
		if (deriv)
		{
			deriv->resize(3 + 2 + 3 * 3 + 3 * 3);
			deriv->setZero();
			(*deriv) += derivatives;
		}

		if (hess)
			throw "Not Implemented!";
	}

	void TFWStretch::addForces(Sheet& sheet,
	    Matrix3N& forces) const
	{
		throw "Not Implemented!";
	}
	void TFWStretch::addForcesPerFace(Sheet& sheet, Eigen::Vector<Real, 23>& forces, int faceId) const
	{
		// ensure that forces is initialized with zeros
		bool _debug = false;
		const std::vector<trimesh::triangle_t>& faces = sheet.getTriangles();
		const Matrix3N& positions = sheet.getPositionsN3();
		const VectorN& amps = sheet.getAmpVec();
		const Matrix2N& dphis = sheet.getDphis1VectorPerFace();
		Real alpha = sheet.getLameAlpha();
		Real beta = sheet.getLameBeta();
		int nfaces = sheet.getFaceCount();
		int nvertices = sheet.getVertexCount();

		Matrix33 facePositions;
		int vidxi = faces[faceId].v[0];
		int vidxj = faces[faceId].v[1];
		int vidxk = faces[faceId].v[2];

		facePositions.col(0) = positions.col(faces[faceId].v[0]);
		facePositions.col(1) = positions.col(faces[faceId].v[1]);
		facePositions.col(2) = positions.col(faces[faceId].v[2]);
		Eigen::Matrix<Real, 4, 9> gradI;
		Eigen::Matrix<Real, 4, 18> gradII; // gradI(i) contains dI(i)/dX
		std::vector<Eigen::Matrix<Real, 9, 9>> hessI;
		Matrix22 I;
		I = firstFundamentalFormNew(facePositions, true, false, &gradI, NULL); // correct till here
		std::vector<Eigen::MatrixXd> hessII;
		std::shared_ptr<SecondFundamentalFormDiscretization> sff;
		sff = std::make_shared<MidedgeAverageFormulation>();
		Matrix22 II = sff->secondFundamentalFormNew(sheet.getBaseMesh(), sheet.getPositionsN3(), sheet.getEdgeDofs(), faceId, &gradII, NULL);
		Matrix22 Ibar = sheet.getRestShape().getIbars()[faceId];
		Matrix22 Ibarinv = Ibar.inverse();
		Matrix22 IIbar = sheet.getRestShape().getIIbars()[faceId];

		Real stretchCoeff = sheet.getThickness() / 4.0 * std::sqrt(Ibar.determinant()) * 0.5;
		// if (faceId == 8)
		// {
		// 	std::cout << "Area " << std::sqrt(Ibar.determinant()) << std::endl;
		// 	std::cout << "Thickness " << sheet.getThickness() << std::endl;
		// 	std::cout << stretchCoeff << std::endl;
		// }
		Vector2 dphi = dphis.col(faceId);

		Real a0 = amps(vidxi);
		Real a1 = amps(vidxj);
		Real a2 = amps(vidxk);
		Vector2 da = Vector2({ a1 - a0, a2 - a0 });
		Matrix22 dadaT = da * da.transpose(); // (a1-a0) * (a1-a0), (a2-a0)*(a1-a0), (a2-a0)*(a1-a0), (a2-a0)*(a2-a0)
		Matrix22 derivdadaTa0;
		derivdadaTa0 << 2 * a0 - 2 * a1, 2 * a0 - a1 - a2, 2 * a0 - a1 - a2, 2 * a0 - 2 * a2;
		Matrix22 derivdadaTa1;
		derivdadaTa1 << 2 * a1 - 2 * a0, a2 - a0, a2 - a0, 0;
		Matrix22 derivdadaTa2;
		derivdadaTa2 << 0, a1 - a0, a1 - a0, 2 * a2 - 2 * a0;

		Matrix22 dphidphiT = dphi * dphi.transpose(); // phi0 * phi0, phi0 * phi1, phi1 * phi0, phi1 * phi1
		Matrix22 derivdphidphiTphi0;
		derivdphidphiTphi0 << 2 * dphi(0), dphi(1), dphi(1), 0;
		Matrix22 derivdphidphiTphi1;
		derivdphidphiTphi1 << 0, dphi(0), dphi(0), 2 * dphi(1);

		Matrix22 dadphiT = da * dphi.transpose(); // (a1-a0) * dphi(0), (a1-a0) * dphi(1), (a2-a0) * dphi(0), (a2-a0) * dphi(1)
		Matrix22 derivdadphiTa0;
		derivdadphiTa0 << -dphi(0), -dphi(1), -dphi(0), -dphi(1);
		Matrix22 derivdadphiTa1;
		derivdadphiTa1 << dphi(0), dphi(1), 0, 0;
		Matrix22 derivdadphiTa2;
		derivdadphiTa2 << 0, 0, dphi(0), dphi(1);
		Matrix22 derivdadphiTphi0;
		derivdadphiTphi0 << (a1 - a0), 0, a2 - a0, 0;
		Matrix22 derivdadphiTphi1;
		derivdadphiTphi1 << 0, a1 - a0, 0, a2 - a0;

		// iterate through quads
		for (int quadId = 0; quadId < 3; quadId++)
		{
			Real u = sheet.getQuadPoints()[quadId].u;
			Real v = sheet.getQuadPoints()[quadId].v;
			Real quadWeight = sheet.getQuadPoints()[quadId].weight;
			Real a = (1 - u - v) * a0 + u * a1 + v * a2;
			if (_debug && faceId == 8)
			{
				std::cout << "Add forces quad " << quadId << std::endl;
				// print debinfo
				std::cout << "a " << a << std::endl;
				std::cout << "Dphi " << dphi.transpose() << std::endl;
				std::cout << "I \n"
				          << I << std::endl;
				std::cout << "Ibar \n"
				          << Ibar << std::endl;
				std::cout << "II \n"
				          << II << std::endl;
			}

			// M1
			Matrix22 M1 = Ibarinv * (I - Ibar + 0.5 * a * a * dphidphiT + 0.5 * dadaT);
			if (_debug && faceId == 8)
			{
				std::cout << "M1\n"
				          << M1 << std::endl;
				// std::cout << "Stretch Coeff " << stretchCoeff << " quadWeight " << quadWeight << std::endl;
			}
			Matrix22 dM1SVdM1 = alpha * M1.trace() * Matrix22::Identity() + 2 * beta * M1;
			// amps derivatives
			Matrix22 dM1da0 = Ibarinv * (dphidphiT * (a * (1 - u - v)) + 0.5 * derivdadaTa0);
			Matrix22 dM1da1 = Ibarinv * (dphidphiT * (a * u) + 0.5 * derivdadaTa1);
			Matrix22 dM1da2 = Ibarinv * (dphidphiT * (a * v) + 0.5 * derivdadaTa2);
			forces(0) += stretchCoeff * quadWeight * doubleContraction(dM1SVdM1, dM1da0);
			forces(1) += stretchCoeff * quadWeight * doubleContraction(dM1SVdM1, dM1da1);
			forces(2) += stretchCoeff * quadWeight * doubleContraction(dM1SVdM1, dM1da2);
			// dphi derivatives
			Matrix22 dM1dphi0 = Ibarinv * (0.5 * a * a * derivdphidphiTphi0);
			Matrix22 dM1dphi1 = Ibarinv * (0.5 * a * a * derivdphidphiTphi1);
			forces(3) += stretchCoeff * quadWeight * doubleContraction(dM1SVdM1, dM1dphi0);
			forces(4) += stretchCoeff * quadWeight * doubleContraction(dM1SVdM1, dM1dphi1);
			// wrinkle params
			// pos derivatives
			for (int xii = 0; xii < 9; xii++)
			{
				Matrix22 dM1dxii = Ibarinv * (gradI.col(xii).reshaped(2, 2));
				forces(5 + xii) += stretchCoeff * quadWeight * doubleContraction(dM1SVdM1, dM1dxii);
			}
			if (dphi.norm() > 1e-4)
			{
				Real IItrace = II.trace();
				Real IbarTrace = Ibar.trace();
				Matrix22 Ibarinvinner = Ibarinv.transpose() * Ibarinv;
				Matrix22 IbarinvIIIbarinv = Ibarinv.transpose() * II * Ibarinv;
				// M2 = 4 * a * a * beta * ((alpha + beta)/(alpha + 2*beta))* ((w^tw * II.trace() - w^T II w)/(w^tw * Ibar.trace() - w^T Ibar w))^2
				// w = Ibarinv * dphi
				Real localCoeff = 4 * beta * ((alpha + beta) / (alpha + 2 * beta));
				Vector2 w = Ibarinv * dphi;
				Real num = (w.dot(w) * II.trace() - w.transpose() * II * w);
				Real den = (w.dot(w) * Ibar.trace() - w.transpose() * Ibar * w);
				// amps derivatives
				forces(0) += stretchCoeff * quadWeight * localCoeff * (num / den) * (num / den) * (2 * a * (1 - u - v));
				forces(1) += stretchCoeff * quadWeight * localCoeff * (num / den) * (num / den) * (2 * a * (u));
				forces(2) += stretchCoeff * quadWeight * localCoeff * (num / den) * (num / den) * (2 * a * (v));
				// dphi derivatives
				Matrix22 T0;
				T0 << 2 * dphi(0), dphi(1), dphi(1), 0;
				Matrix22 T1;
				T1 << 0, dphi(0), dphi(0), 2 * dphi(1);
				Real dnumdphi0 = doubleContraction(IItrace * Ibarinvinner - IbarinvIIIbarinv, T0);
				Real ddendphi0 = doubleContraction(IbarTrace * Ibarinvinner - Ibarinv, T0);
				Real dnumdphi1 = doubleContraction(IItrace * Ibarinvinner - IbarinvIIIbarinv, T1);
				Real ddendphi1 = doubleContraction(IbarTrace * Ibarinvinner - Ibarinv, T1);
				forces(3) += stretchCoeff * quadWeight * a * a * localCoeff * 2 * (num / (den * den * den)) * (dnumdphi0 * den - ddendphi0 * num);
				forces(4) += stretchCoeff * quadWeight * a * a * localCoeff * 2 * (num / (den * den * den)) * (dnumdphi1 * den - ddendphi1 * num);
				// pos derivatives
				Real posDerivCoeff = a * a * localCoeff * 2 * num / (den * den);
				Matrix22 dnumdII = w.dot(w) * Matrix22::Identity() - w * w.transpose();
				for (int xii = 0; xii < 9; xii++)
				{
					forces(5 + xii) += stretchCoeff * quadWeight * posDerivCoeff * doubleContraction(dnumdII, gradII.col(xii).reshaped(2, 2));
					forces(14 + xii) += stretchCoeff * quadWeight * posDerivCoeff * doubleContraction(dnumdII, gradII.col(xii + 9).reshaped(2, 2));
				}
			}
			// M3 = (1/std::sqrt(32))* a * Ibarinv * (dadphiT + dadphiTT);
			Matrix22 M3sub = Ibarinv * (dadphiT + dadphiT.transpose());
			Matrix22 dM3subda0 = Ibarinv * (derivdadphiTa0 + derivdadphiTa0.transpose());
			Matrix22 dM3subda1 = Ibarinv * (derivdadphiTa1 + derivdadphiTa1.transpose());
			Matrix22 dM3subda2 = Ibarinv * (derivdadphiTa2 + derivdadphiTa2.transpose());
			Matrix22 dM3subddphi0 = Ibarinv * (derivdadphiTphi0 + derivdadphiTphi0.transpose());
			Matrix22 dM3subddphi1 = Ibarinv * (derivdadphiTphi1 + derivdadphiTphi1.transpose());
			Matrix22 M3 = a * M3sub;
			Matrix22 dM3SVdM3 = alpha * M3.trace() * Matrix22::Identity() + 2 * beta * M3;
			// amps derivatives
			forces(0) += stretchCoeff * quadWeight * (1 / 32.0) * doubleContraction(dM3SVdM3, (1 - u - v) * M3sub + a * dM3subda0);
			forces(1) += stretchCoeff * quadWeight * (1 / 32.0) * doubleContraction(dM3SVdM3, u * M3sub + a * dM3subda1);
			forces(2) += stretchCoeff * quadWeight * (1 / 32.0) * doubleContraction(dM3SVdM3, v * M3sub + a * dM3subda2);
			// dphi derivatives
			forces(3) += stretchCoeff * quadWeight * (1 / 32.0) * doubleContraction(dM3SVdM3, a * dM3subddphi0);
			forces(4) += stretchCoeff * quadWeight * (1 / 32.0) * doubleContraction(dM3SVdM3, a * dM3subddphi1);

			// M4
			Matrix22 M4 = (Ibarinv * dadaT);
			Matrix22 dM4SVdM4 = alpha * M4.trace() * Matrix22::Identity() + 2 * beta * M4;
			// remember 1/8
			// amps derivatives
			forces(0) += stretchCoeff * quadWeight * (1 / 8.0) * doubleContraction(dM4SVdM4, Ibarinv * derivdadaTa0);
			forces(1) += stretchCoeff * quadWeight * (1 / 8.0) * doubleContraction(dM4SVdM4, Ibarinv * derivdadaTa1);
			forces(2) += stretchCoeff * quadWeight * (1 / 8.0) * doubleContraction(dM4SVdM4, Ibarinv * derivdadaTa2);
		}
	}

	void TFWStretch::addForces(Sheet& sheet, VectorN& forces) const
	{
		// ensure that forces is initialized with zeros
		bool _debug = false;
		const std::vector<trimesh::triangle_t>& faces = sheet.getTriangles();
		const Matrix3N& positions = sheet.getPositionsN3();
		const VectorN& amps = sheet.getAmpVec();
		const Matrix2N& dphis = sheet.getDphis1VectorPerFace();
		Real alpha = sheet.getLameAlpha();
		Real beta = sheet.getLameBeta();
		int nfaces = sheet.getFaceCount();
		int nvertices = sheet.getVertexCount();

		std::vector<Eigen::Vector<Real, 23>> faceForces;
		faceForces.resize(nfaces);
#ifdef USE_OMP
#pragma omp parallel for
#endif
		for (int faceId = 0; faceId < nfaces; faceId++)
		{
			// enable this for testing the per face energy functions
			// 		Eigen::Vector<Real, 23> faceDeriv;
			// faceDeriv.setZero();
			// addForcesPerFace(sheet, faceDeriv, faceId);
			// getEnergyPerFace(sheet, &faceDeriv, faceId);
			// original code
			Eigen::Vector<Real, 23> faceDeriv;
			faceDeriv.setZero();

			Matrix33 facePositions;
			int vidxi = faces[faceId].v[0];
			int vidxj = faces[faceId].v[1];
			int vidxk = faces[faceId].v[2];

			facePositions.col(0) = positions.col(faces[faceId].v[0]);
			facePositions.col(1) = positions.col(faces[faceId].v[1]);
			facePositions.col(2) = positions.col(faces[faceId].v[2]);
			Eigen::Matrix<Real, 4, 9> gradI;
			Eigen::Matrix<Real, 4, 18> gradII; // gradI(i) contains dI(i)/dX
			Matrix22 I;
			I = firstFundamentalFormNew(facePositions, true, false, &gradI, NULL); // correct till here
			std::shared_ptr<SecondFundamentalFormDiscretization> sff;
			sff = std::make_shared<MidedgeAverageFormulation>();
			Matrix22 II = sff->secondFundamentalFormNew(sheet.getBaseMesh(), sheet.getPositionsN3(), sheet.getEdgeDofs(), faceId, &gradII, NULL);
			Matrix22 Ibar = sheet.getRestShape().getIbars()[faceId];
			Matrix22 Ibarinv = Ibar.inverse();
			Matrix22 IIbar = sheet.getRestShape().getIIbars()[faceId];

			Real stretchCoeff = sheet.getThickness() / 4.0 * std::sqrt(Ibar.determinant()) * 0.5;
			// if (faceId == 8)
			// {
			// 	std::cout << "Area " << std::sqrt(Ibar.determinant()) << std::endl;
			// 	std::cout << "Thickness " << sheet.getThickness() << std::endl;
			// 	std::cout << stretchCoeff << std::endl;
			// }
			Vector2 dphi = dphis.col(faceId);

			Real a0 = amps(vidxi);
			Real a1 = amps(vidxj);
			Real a2 = amps(vidxk);
			Vector2 da = Vector2({ a1 - a0, a2 - a0 });
			Matrix22 dadaT = da * da.transpose(); // (a1-a0) * (a1-a0), (a2-a0)*(a1-a0), (a2-a0)*(a1-a0), (a2-a0)*(a2-a0)
			Matrix22 derivdadaTa0;
			derivdadaTa0 << 2 * a0 - 2 * a1, 2 * a0 - a1 - a2, 2 * a0 - a1 - a2, 2 * a0 - 2 * a2;
			Matrix22 derivdadaTa1;
			derivdadaTa1 << 2 * a1 - 2 * a0, a2 - a0, a2 - a0, 0;
			Matrix22 derivdadaTa2;
			derivdadaTa2 << 0, a1 - a0, a1 - a0, 2 * a2 - 2 * a0;

			Matrix22 dphidphiT = dphi * dphi.transpose(); // phi0 * phi0, phi0 * phi1, phi1 * phi0, phi1 * phi1
			Matrix22 derivdphidphiTphi0;
			derivdphidphiTphi0 << 2 * dphi(0), dphi(1), dphi(1), 0;
			Matrix22 derivdphidphiTphi1;
			derivdphidphiTphi1 << 0, dphi(0), dphi(0), 2 * dphi(1);

			Matrix22 dadphiT = da * dphi.transpose(); // (a1-a0) * dphi(0), (a1-a0) * dphi(1), (a2-a0) * dphi(0), (a2-a0) * dphi(1)
			Matrix22 derivdadphiTa0;
			derivdadphiTa0 << -dphi(0), -dphi(1), -dphi(0), -dphi(1);
			Matrix22 derivdadphiTa1;
			derivdadphiTa1 << dphi(0), dphi(1), 0, 0;
			Matrix22 derivdadphiTa2;
			derivdadphiTa2 << 0, 0, dphi(0), dphi(1);
			Matrix22 derivdadphiTphi0;
			derivdadphiTphi0 << (a1 - a0), 0, a2 - a0, 0;
			Matrix22 derivdadphiTphi1;
			derivdadphiTphi1 << 0, a1 - a0, 0, a2 - a0;

			// iterate through quads
			for (int quadId = 0; quadId < 3; quadId++)
			{
				Real u = sheet.getQuadPoints()[quadId].u;
				Real v = sheet.getQuadPoints()[quadId].v;
				Real quadWeight = sheet.getQuadPoints()[quadId].weight;
				Real a = (1 - u - v) * a0 + u * a1 + v * a2;
				if (_debug && faceId == 8)
				{
					std::cout << "Add faceDeriv quad " << quadId << std::endl;
					// print debinfo
					std::cout << "a " << a << std::endl;
					std::cout << "Dphi " << dphi.transpose() << std::endl;
					std::cout << "I \n"
					          << I << std::endl;
					std::cout << "Ibar \n"
					          << Ibar << std::endl;
					std::cout << "II \n"
					          << II << std::endl;
				}

				// M1
				Matrix22 M1 = Ibarinv * (I - Ibar + 0.5 * a * a * dphidphiT + 0.5 * dadaT);
				if (_debug && faceId == 8)
				{
					std::cout << "M1\n"
					          << M1 << std::endl;
					// std::cout << "Stretch Coeff " << stretchCoeff << " quadWeight " << quadWeight << std::endl;
				}
				Matrix22 dM1SVdM1 = alpha * M1.trace() * Matrix22::Identity() + 2 * beta * M1;
				// amps derivatives
				Matrix22 dM1da0 = Ibarinv * (dphidphiT * (a * (1 - u - v)) + 0.5 * derivdadaTa0);
				Matrix22 dM1da1 = Ibarinv * (dphidphiT * (a * u) + 0.5 * derivdadaTa1);
				Matrix22 dM1da2 = Ibarinv * (dphidphiT * (a * v) + 0.5 * derivdadaTa2);
				faceDeriv(0) += stretchCoeff * quadWeight * doubleContraction(dM1SVdM1, dM1da0);
				faceDeriv(1) += stretchCoeff * quadWeight * doubleContraction(dM1SVdM1, dM1da1);
				faceDeriv(2) += stretchCoeff * quadWeight * doubleContraction(dM1SVdM1, dM1da2);
				// dphi derivatives
				Matrix22 dM1dphi0 = Ibarinv * (0.5 * a * a * derivdphidphiTphi0);
				Matrix22 dM1dphi1 = Ibarinv * (0.5 * a * a * derivdphidphiTphi1);
				faceDeriv(3) += stretchCoeff * quadWeight * doubleContraction(dM1SVdM1, dM1dphi0);
				faceDeriv(4) += stretchCoeff * quadWeight * doubleContraction(dM1SVdM1, dM1dphi1);
				// wrinkle params
				// pos derivatives
				for (int xii = 0; xii < 9; xii++)
				{
					Matrix22 dM1dxii = Ibarinv * (gradI.col(xii).reshaped(2, 2));
					faceDeriv(5 + xii) += stretchCoeff * quadWeight * doubleContraction(dM1SVdM1, dM1dxii);
				}
				if (dphi.norm() > 1e-4)
				{
					Real IItrace = II.trace();
					Real IbarTrace = Ibar.trace();
					Matrix22 Ibarinvinner = Ibarinv.transpose() * Ibarinv;
					Matrix22 IbarinvIIIbarinv = Ibarinv.transpose() * II * Ibarinv;
					// M2 = 4 * a * a * beta * ((alpha + beta)/(alpha + 2*beta))* ((w^tw * II.trace() - w^T II w)/(w^tw * Ibar.trace() - w^T Ibar w))^2
					// w = Ibarinv * dphi
					Real localCoeff = 4 * beta * ((alpha + beta) / (alpha + 2 * beta));
					Vector2 w = Ibarinv * dphi;
					Real num = (w.dot(w) * II.trace() - w.transpose() * II * w);
					Real den = (w.dot(w) * Ibar.trace() - w.transpose() * Ibar * w);
					// amps derivatives
					faceDeriv(0) += stretchCoeff * quadWeight * localCoeff * (num / den) * (num / den) * (2 * a * (1 - u - v));
					faceDeriv(1) += stretchCoeff * quadWeight * localCoeff * (num / den) * (num / den) * (2 * a * (u));
					faceDeriv(2) += stretchCoeff * quadWeight * localCoeff * (num / den) * (num / den) * (2 * a * (v));
					// dphi derivatives
					Matrix22 T0;
					T0 << 2 * dphi(0), dphi(1), dphi(1), 0;
					Matrix22 T1;
					T1 << 0, dphi(0), dphi(0), 2 * dphi(1);
					Real dnumdphi0 = doubleContraction(IItrace * Ibarinvinner - IbarinvIIIbarinv, T0);
					Real ddendphi0 = doubleContraction(IbarTrace * Ibarinvinner - Ibarinv, T0);
					Real dnumdphi1 = doubleContraction(IItrace * Ibarinvinner - IbarinvIIIbarinv, T1);
					Real ddendphi1 = doubleContraction(IbarTrace * Ibarinvinner - Ibarinv, T1);
					faceDeriv(3) += stretchCoeff * quadWeight * a * a * localCoeff * 2 * (num / (den * den * den)) * (dnumdphi0 * den - ddendphi0 * num);
					faceDeriv(4) += stretchCoeff * quadWeight * a * a * localCoeff * 2 * (num / (den * den * den)) * (dnumdphi1 * den - ddendphi1 * num);
					// pos derivatives
					Real posDerivCoeff = a * a * localCoeff * 2 * num / (den * den);
					Matrix22 dnumdII = w.dot(w) * Matrix22::Identity() - w * w.transpose();
					for (int xii = 0; xii < 9; xii++)
					{
						faceDeriv(5 + xii) += stretchCoeff * quadWeight * posDerivCoeff * doubleContraction(dnumdII, gradII.col(xii).reshaped(2, 2));
						faceDeriv(14 + xii) += stretchCoeff * quadWeight * posDerivCoeff * doubleContraction(dnumdII, gradII.col(xii + 9).reshaped(2, 2));
					}
				}
				// M3 = (1/std::sqrt(32))* a * Ibarinv * (dadphiT + dadphiTT);
				Matrix22 M3sub = Ibarinv * (dadphiT + dadphiT.transpose());
				Matrix22 dM3subda0 = Ibarinv * (derivdadphiTa0 + derivdadphiTa0.transpose());
				Matrix22 dM3subda1 = Ibarinv * (derivdadphiTa1 + derivdadphiTa1.transpose());
				Matrix22 dM3subda2 = Ibarinv * (derivdadphiTa2 + derivdadphiTa2.transpose());
				Matrix22 dM3subddphi0 = Ibarinv * (derivdadphiTphi0 + derivdadphiTphi0.transpose());
				Matrix22 dM3subddphi1 = Ibarinv * (derivdadphiTphi1 + derivdadphiTphi1.transpose());
				Matrix22 M3 = a * M3sub;
				Matrix22 dM3SVdM3 = alpha * M3.trace() * Matrix22::Identity() + 2 * beta * M3;
				// amps derivatives
				// dM3da0 = (1-u-v) *M3sub + a * dM3subda0;
				// dM3da1 = u *M3sub + a * dM3subda1;
				// dM3da2 = v *M3sub + a * dM3subda2;
				faceDeriv(0) += stretchCoeff * quadWeight * (1 / 32.0) * doubleContraction(dM3SVdM3, (1 - u - v) * M3sub + a * dM3subda0);
				faceDeriv(1) += stretchCoeff * quadWeight * (1 / 32.0) * doubleContraction(dM3SVdM3, u * M3sub + a * dM3subda1);
				faceDeriv(2) += stretchCoeff * quadWeight * (1 / 32.0) * doubleContraction(dM3SVdM3, v * M3sub + a * dM3subda2);
				// dphi derivatives
				faceDeriv(3) += stretchCoeff * quadWeight * (1 / 32.0) * doubleContraction(dM3SVdM3, a * dM3subddphi0);
				faceDeriv(4) += stretchCoeff * quadWeight * (1 / 32.0) * doubleContraction(dM3SVdM3, a * dM3subddphi1);

				// M4
				Matrix22 M4 = (Ibarinv * dadaT);
				Matrix22 dM4SVdM4 = alpha * M4.trace() * Matrix22::Identity() + 2 * beta * M4;
				// remember 1/8
				// amps derivatives
				faceDeriv(0) += stretchCoeff * quadWeight * (1 / 8.0) * doubleContraction(dM4SVdM4, Ibarinv * derivdadaTa0);
				faceDeriv(1) += stretchCoeff * quadWeight * (1 / 8.0) * doubleContraction(dM4SVdM4, Ibarinv * derivdadaTa1);
				faceDeriv(2) += stretchCoeff * quadWeight * (1 / 8.0) * doubleContraction(dM4SVdM4, Ibarinv * derivdadaTa2);
			}

			faceForces.at(faceId) = -faceDeriv;
		}

		// debug
		if (_debug)
		{
			std::cout << "addForces() derivative for face 8\n";
			std::cout << faceForces.at(8) << std::endl;
			exit(0);
		}

		const std::vector<trimesh::triangle_t>& edgeIDOppVertexIDPerface = sheet.getVertexOppositeEdgesPerFace();

		const Eigen::MatrixXi& vtexIDOppositeEdge = sheet.getEdgeOppositeVertices();

		// now we have derivatives defined for each face, we need to aggregate derivatives in the forces matrix
		// serial loop
		for (int faceId = 0; faceId < nfaces; faceId++)
		{
			int vidxi = faces[faceId].v[0];
			int vidxj = faces[faceId].v[1];
			int vidxk = faces[faceId].v[2];
			// fill derivatives
			Eigen::Vector<Real, 23> faceForce = faceForces.at(faceId);
			// amp derivatives
			forces(vidxi) += faceForce(0);
			forces(vidxj) += faceForce(1);
			forces(vidxk) += faceForce(2);
			// dphi derivatives
			forces(nvertices + 2 * faceId) += faceForce(3);
			forces(nvertices + 2 * faceId + 1) += faceForce(4);
			// position derivatives
			for (int idxi = 0; idxi < 3; idxi++)
			{
				// indices for current face
				int vidx = faces[faceId].v[idxi];
				forces(nvertices + 2 * nfaces + 3 * vidx) += faceForce(5 + 3 * idxi);
				forces(nvertices + 2 * nfaces + 3 * vidx + 1) += faceForce(5 + 3 * idxi + 1);
				forces(nvertices + 2 * nfaces + 3 * vidx + 2) += faceForce(5 + 3 * idxi + 2);
				// indices for adjacent face vertices
				int edgeId = edgeIDOppVertexIDPerface[faceId].v[idxi];
				int vidTilda = vtexIDOppositeEdge(edgeId, 0) == vidx ? vtexIDOppositeEdge(edgeId, 1) : vtexIDOppositeEdge(edgeId, 0);
				if (vidTilda != -1)
				{
					forces(nvertices + 2 * nfaces + 3 * vidTilda) += faceForce(14 + 3 * idxi);
					forces(nvertices + 2 * nfaces + 3 * vidTilda + 1) += faceForce(14 + 3 * idxi + 1);
					forces(nvertices + 2 * nfaces + 3 * vidTilda + 2) += faceForce(14 + 3 * idxi + 2);
				}
			}
		}
	}
	Real TFWStretch::getEnergyPerFace(Sheet& sheet, Eigen::Vector<Real, 23>* forces, int faceId) const
	{
		// ensure that forces is initialized with zeros
		bool _debug = false;
		const std::vector<trimesh::triangle_t>& faces = sheet.getTriangles();
		const Matrix3N& positions = sheet.getPositionsN3();
		const VectorN& amps = sheet.getAmpVec();
		const Matrix2N& dphis = sheet.getDphis1VectorPerFace();
		Real alpha = sheet.getLameAlpha();
		Real beta = sheet.getLameBeta();
		int nfaces = sheet.getFaceCount();
		int nvertices = sheet.getVertexCount();

		Matrix33 facePositions;
		int vidxi = faces[faceId].v[0];
		int vidxj = faces[faceId].v[1];
		int vidxk = faces[faceId].v[2];

		facePositions.col(0) = positions.col(faces[faceId].v[0]);
		facePositions.col(1) = positions.col(faces[faceId].v[1]);
		facePositions.col(2) = positions.col(faces[faceId].v[2]);
		Eigen::Matrix<Real, 4, 9> gradI;
		Eigen::Matrix<Real, 4, 18> gradII; // gradI(i) contains dI(i)/dX
		std::vector<Eigen::Matrix<Real, 9, 9>> hessI;
		Matrix22 I;
		I = firstFundamentalFormNew(facePositions, true, false, &gradI, NULL); // correct till here
		std::vector<Eigen::MatrixXd> hessII;
		std::shared_ptr<SecondFundamentalFormDiscretization> sff;
		sff = std::make_shared<MidedgeAverageFormulation>();
		Matrix22 II = sff->secondFundamentalFormNew(sheet.getBaseMesh(), sheet.getPositionsN3(), sheet.getEdgeDofs(), faceId, &gradII, NULL);
		Matrix22 Ibar = sheet.getRestShape().getIbars()[faceId];
		Matrix22 Ibarinv = Ibar.inverse();
		Matrix22 IIbar = sheet.getRestShape().getIIbars()[faceId];

		Real stretchCoeff = sheet.getThickness() / 4.0 * std::sqrt(Ibar.determinant()) * 0.5;
		// if (faceId == 8)
		// {
		// 	std::cout << "Area " << std::sqrt(Ibar.determinant()) << std::endl;
		// 	std::cout << "Thickness " << sheet.getThickness() << std::endl;
		// 	std::cout << stretchCoeff << std::endl;
		// }
		Vector2 dphi = dphis.col(faceId);

		Real a0 = amps(vidxi);
		Real a1 = amps(vidxj);
		Real a2 = amps(vidxk);
		Vector2 da = Vector2({ a1 - a0, a2 - a0 });
		Matrix22 dadaT = da * da.transpose(); // (a1-a0) * (a1-a0), (a2-a0)*(a1-a0), (a2-a0)*(a1-a0), (a2-a0)*(a2-a0)
		Matrix22 derivdadaTa0;
		derivdadaTa0 << 2 * a0 - 2 * a1, 2 * a0 - a1 - a2, 2 * a0 - a1 - a2, 2 * a0 - 2 * a2;
		Matrix22 derivdadaTa1;
		derivdadaTa1 << 2 * a1 - 2 * a0, a2 - a0, a2 - a0, 0;
		Matrix22 derivdadaTa2;
		derivdadaTa2 << 0, a1 - a0, a1 - a0, 2 * a2 - 2 * a0;

		Matrix22 dphidphiT = dphi * dphi.transpose(); // phi0 * phi0, phi0 * phi1, phi1 * phi0, phi1 * phi1
		Matrix22 derivdphidphiTphi0;
		derivdphidphiTphi0 << 2 * dphi(0), dphi(1), dphi(1), 0;
		Matrix22 derivdphidphiTphi1;
		derivdphidphiTphi1 << 0, dphi(0), dphi(0), 2 * dphi(1);

		Matrix22 dadphiT = da * dphi.transpose(); // (a1-a0) * dphi(0), (a1-a0) * dphi(1), (a2-a0) * dphi(0), (a2-a0) * dphi(1)
		Matrix22 derivdadphiTa0;
		derivdadphiTa0 << -dphi(0), -dphi(1), -dphi(0), -dphi(1);
		Matrix22 derivdadphiTa1;
		derivdadphiTa1 << dphi(0), dphi(1), 0, 0;
		Matrix22 derivdadphiTa2;
		derivdadphiTa2 << 0, 0, dphi(0), dphi(1);
		Matrix22 derivdadphiTphi0;
		derivdadphiTphi0 << (a1 - a0), 0, a2 - a0, 0;
		Matrix22 derivdadphiTphi1;
		derivdadphiTphi1 << 0, a1 - a0, 0, a2 - a0;
		Real energy = 0;
		// iterate through quads
		for (int quadId = 0; quadId < 3; quadId++)
		{
			Real u = sheet.getQuadPoints()[quadId].u;
			Real v = sheet.getQuadPoints()[quadId].v;
			Real quadWeight = sheet.getQuadPoints()[quadId].weight;
			Real a = (1 - u - v) * a0 + u * a1 + v * a2;
			if (_debug && faceId == 8)
			{
				std::cout << "Add forces quad " << quadId << std::endl;
				// print debinfo
				std::cout << "a " << a << std::endl;
				std::cout << "Dphi " << dphi.transpose() << std::endl;
				std::cout << "I \n"
				          << I << std::endl;
				std::cout << "Ibar \n"
				          << Ibar << std::endl;
				std::cout << "II \n"
				          << II << std::endl;
			}

			// M1
			Matrix22 M1 = Ibarinv * (I - Ibar + 0.5 * a * a * dphidphiT + 0.5 * dadaT);
			energy += stretchCoeff * quadWeight * SVNormSquared(alpha, beta, M1);
			if (_debug && faceId == 8)
			{
				std::cout << "M1\n"
				          << M1 << std::endl;
				// std::cout << "Stretch Coeff " << stretchCoeff << " quadWeight " << quadWeight << std::endl;
			}
			if (forces)
			{
				Matrix22 dM1SVdM1 = alpha * M1.trace() * Matrix22::Identity() + 2 * beta * M1;
				// amps derivatives
				// Matrix22 dM1da0 = Ibarinv * (dphidphiT * (a * (1 - u - v)) + 0.5 * derivdadaTa0);
				// Matrix22 dM1da1 = Ibarinv * (dphidphiT * (a * u) + 0.5 * derivdadaTa1);
				// Matrix22 dM1da2 = Ibarinv * (dphidphiT * (a * v) + 0.5 * derivdadaTa2);
				// forces->coeffRef(0) += stretchCoeff * quadWeight * doubleContraction(dM1SVdM1, dM1da0);
				// forces->coeffRef(1) += stretchCoeff * quadWeight * doubleContraction(dM1SVdM1, dM1da1);
				// forces->coeffRef(2) += stretchCoeff * quadWeight * doubleContraction(dM1SVdM1, dM1da2);
				// dphi derivatives
				Matrix22 dM1dphi0 = Ibarinv * (0.5 * a * a * derivdphidphiTphi0);
				Matrix22 dM1dphi1 = Ibarinv * (0.5 * a * a * derivdphidphiTphi1);
				forces->coeffRef(3) += stretchCoeff * quadWeight * doubleContraction(dM1SVdM1, dM1dphi0);
				forces->coeffRef(4) += stretchCoeff * quadWeight * doubleContraction(dM1SVdM1, dM1dphi1);
				// wrinkle params
				// pos derivatives
				// for (int xii = 0; xii < 9; xii++)
				// {
				// 	Matrix22 dM1dxii = Ibarinv * (gradI.col(xii).reshaped(2, 2));
				// 	forces->coeffRef(5 + xii) += stretchCoeff * quadWeight * doubleContraction(dM1SVdM1, dM1dxii);
				// }
			}
			if (dphi.norm() > 1e-4)
			{
				Real IItrace = II.trace();
				Real IbarTrace = Ibar.trace();
				Matrix22 Ibarinvinner = Ibarinv.transpose() * Ibarinv;
				Matrix22 IbarinvIIIbarinv = Ibarinv.transpose() * II * Ibarinv;
				// M2 = 4 * a * a * beta * ((alpha + beta)/(alpha + 2*beta))* ((w^tw * II.trace() - w^T II w)/(w^tw * Ibar.trace() - w^T Ibar w))^2
				// w = Ibarinv * dphi
				Real localCoeff = 4 * beta * ((alpha + beta) / (alpha + 2 * beta));
				Vector2 w = Ibarinv * dphi;
				Real num = (w.dot(w) * II.trace() - w.transpose() * II * w);
				Real den = (w.dot(w) * Ibar.trace() - w.transpose() * Ibar * w);
				energy += stretchCoeff * quadWeight * localCoeff * a * a * (num * num) / (den * den);
				if (forces)
				{ // amps derivatives
					// forces->coeffRef(0) += stretchCoeff * quadWeight * localCoeff * (num / den) * (num / den) * (2 * a * (1 - u - v));
					// forces->coeffRef(1) += stretchCoeff * quadWeight * localCoeff * (num / den) * (num / den) * (2 * a * (u));
					// forces->coeffRef(2) += stretchCoeff * quadWeight * localCoeff * (num / den) * (num / den) * (2 * a * (v));
					// dphi derivatives
					Matrix22 T0;
					T0 << 2 * dphi(0), dphi(1), dphi(1), 0;
					Matrix22 T1;
					T1 << 0, dphi(0), dphi(0), 2 * dphi(1);
					Real dnumdphi0 = doubleContraction(IItrace * Ibarinvinner - IbarinvIIIbarinv, T0);
					Real ddendphi0 = doubleContraction(IbarTrace * Ibarinvinner - Ibarinv, T0);
					Real dnumdphi1 = doubleContraction(IItrace * Ibarinvinner - IbarinvIIIbarinv, T1);
					Real ddendphi1 = doubleContraction(IbarTrace * Ibarinvinner - Ibarinv, T1);
					forces->coeffRef(3) += stretchCoeff * quadWeight * a * a * localCoeff * 2 * (num / (den * den * den)) * (dnumdphi0 * den - ddendphi0 * num);
					forces->coeffRef(4) += stretchCoeff * quadWeight * a * a * localCoeff * 2 * (num / (den * den * den)) * (dnumdphi1 * den - ddendphi1 * num);
					// pos derivatives
					// Real posDerivCoeff = a * a * localCoeff * 2 * num / (den * den);
					// Matrix22 dnumdII = w.dot(w) * Matrix22::Identity() - w * w.transpose();
					// for (int xii = 0; xii < 9; xii++)
					// {
					// 	forces->coeffRef(5 + xii) += stretchCoeff * quadWeight * posDerivCoeff * doubleContraction(dnumdII, gradII.col(xii).reshaped(2, 2));
					// 	forces->coeffRef(14 + xii) += stretchCoeff * quadWeight * posDerivCoeff * doubleContraction(dnumdII, gradII.col(xii + 9).reshaped(2, 2));
					// }
				}
			}
			// M3 = (1/std::sqrt(32))* a * Ibarinv * (dadphiT + dadphiTT);
			Matrix22 M3sub = Ibarinv * (dadphiT + dadphiT.transpose());
			energy += stretchCoeff * quadWeight * SVNormSquared(alpha, beta, M3sub) * a * a / 32.0;
			if (forces)
			{
				// Matrix22 dM3subda0 = Ibarinv * (derivdadphiTa0 + derivdadphiTa0.transpose());
				// Matrix22 dM3subda1 = Ibarinv * (derivdadphiTa1 + derivdadphiTa1.transpose());
				// Matrix22 dM3subda2 = Ibarinv * (derivdadphiTa2 + derivdadphiTa2.transpose());
				Matrix22 dM3subddphi0 = Ibarinv * (derivdadphiTphi0 + derivdadphiTphi0.transpose());
				Matrix22 dM3subddphi1 = Ibarinv * (derivdadphiTphi1 + derivdadphiTphi1.transpose());
				Matrix22 M3 = a * M3sub;
				Matrix22 dM3SVdM3 = alpha * M3.trace() * Matrix22::Identity() + 2 * beta * M3;
				// amps derivatives

				// forces->coeffRef(0) += stretchCoeff * quadWeight * (1 / 32.0) * doubleContraction(dM3SVdM3, (1 - u - v) * M3sub + a * dM3subda0);
				// forces->coeffRef(1) += stretchCoeff * quadWeight * (1 / 32.0) * doubleContraction(dM3SVdM3, u * M3sub + a * dM3subda1);
				// forces->coeffRef(2) += stretchCoeff * quadWeight * (1 / 32.0) * doubleContraction(dM3SVdM3, v * M3sub + a * dM3subda2);
				// dphi derivatives
				forces->coeffRef(3) += stretchCoeff * quadWeight * (1 / 32.0) * doubleContraction(dM3SVdM3, a * dM3subddphi0);
				forces->coeffRef(4) += stretchCoeff * quadWeight * (1 / 32.0) * doubleContraction(dM3SVdM3, a * dM3subddphi1);
			}
			// M4
			Matrix22 M4 = (Ibarinv * dadaT);
			energy += stretchCoeff * quadWeight * (1 / 8.0) * SVNormSquared(alpha, beta, M4);
			if (forces)
			{
				// Matrix22 dM4SVdM4 = alpha * M4.trace() * Matrix22::Identity() + 2 * beta * M4;
				// // remember 1/8
				// // amps derivatives
				// forces->coeffRef(0) += stretchCoeff * quadWeight * (1 / 8.0) * doubleContraction(dM4SVdM4, Ibarinv * derivdadaTa0);
				// forces->coeffRef(1) += stretchCoeff * quadWeight * (1 / 8.0) * doubleContraction(dM4SVdM4, Ibarinv * derivdadaTa1);
				// forces->coeffRef(2) += stretchCoeff * quadWeight * (1 / 8.0) * doubleContraction(dM4SVdM4, Ibarinv * derivdadaTa2);
			}
		}
		return energy;
	}
	// delete
	const void TFWStretch::printDebInfo(Sheet& sheet, int faceId, int quadId) const
	{
		std::cout << "printing debinfo for face " << faceId << " quad " << quadId << std::endl;
		const std::vector<trimesh::triangle_t>& faces = sheet.getTriangles();
		const std::vector<trimesh::edge_t>& edges = sheet.getEdges();
		const Matrix3N& positions = sheet.getPositionsN3();
		const Eigen::MatrixXi& edgeOppositeVertices = sheet.getEdgeOppositeVertices();
		const std::vector<trimesh::triangle_t>& vertexOppositeEdges = sheet.getVertexOppositeEdgesPerFace();
		int nfaces = faces.size();
		int nedges = edges.size();
		int nvertices = positions.cols();
		Matrix33 facePositions;
		int vidxi = faces[faceId].v[0];
		int vidxj = faces[faceId].v[1];
		int vidxk = faces[faceId].v[2];
		// facePositions.col(0) = positions.block<3, 1>(0, faces[faceId].v[0]);
		// facePositions.col(1) = positions.block<3, 1>(0, faces[faceId].v[1]);
		// facePositions.col(2) = positions.block<3, 1>(0, faces[faceId].v[2]);
		facePositions.col(0) = positions.col(faces[faceId].v[0]);
		facePositions.col(1) = positions.col(faces[faceId].v[1]);
		facePositions.col(2) = positions.col(faces[faceId].v[2]);
		MatrixNN gradI, gradII; // gradI(i) contains dI(i)/dX
		std::vector<Eigen::Matrix<Real, 9, 9>> hessI;
		Matrix22 I;
		I = firstFundamentalForm(facePositions, false, false, NULL, NULL); // correct till here
		std::vector<Eigen::MatrixXd> hessII;
		std::shared_ptr<SecondFundamentalFormDiscretization> sff;
		sff = std::make_shared<MidedgeAverageFormulation>();
		Matrix22 II = sff->secondFundamentalForm(sheet.getBaseMesh(), sheet.getPositionsN3().transpose(), sheet.getEdgeDofs(), faceId, NULL, NULL);
		Matrix22 Ibar = sheet.getRestShape().getIbars()[faceId];
		Matrix22 Ibarinv = Ibar.inverse();
		Matrix22 IIbar = sheet.getRestShape().getIIbars()[faceId];
		Vector3 gradA; // datilda/da
		Vector2 da;
		Matrix23 gradDA; // dda/da
		Real u = sheet.getQuadPoints()[quadId].u;
		Real v = sheet.getQuadPoints()[quadId].v;
		Vector3 ampUVWeights;
		ampUVWeights << 1 - u - v, u, v;
		Real a = sheet.computeAmplitudesFromQuad(faceId, quadId, &da, NULL, NULL, nullptr, nullptr);
		Matrix22 gradDphi;
		Vector2 dphi = sheet.computeDphi1VectorPerFace(faceId, NULL);
		std::vector<Matrix22> gradDphiDphiTensor;
		std::vector<Matrix22> hessDphiDphiT;
		Matrix22 dphidphiTensor;
		dphidphiTensor = sheet.computeDphiDphi1VectorPerFaceTensor(faceId, quadId, NULL, NULL);
		std::vector<Matrix22> gradDaDaTensor;
		std::vector<Matrix33> hessDaDa;
		Matrix22 dadaTensor;
		dadaTensor = sheet.computeDaDaTensor(faceId, quadId, NULL, NULL);
		std::vector<Matrix22> gradDaDphiTensor;
		std::vector<Eigen::Matrix<Real, 5, 5>> hessDaDphiTensor, hessDphiDaTensor;
		Matrix22 dadphiTensor = sheet.computeDaDphi1VectorPerFaceTensor(faceId, quadId, NULL, NULL);
		std::vector<Matrix22> gradDphiDaTensor;
		Matrix22 dphidaTensor = sheet.computeDphi1VectorPerFaceDaTensor(faceId, quadId, NULL, NULL);
		Real lameAlpha = sheet.getLameAlpha();
		Real lameBeta = sheet.getLameBeta();

		Real energy = 0.0;

		// Second term
		Real w2se = 0;
		// Derivatives are correct numerically
		if (dphi.norm() > 1e-4)
		{
#ifdef ARGUS_USE_FLOATS
			typedef TinyAD::Float<9> AReal;
#else
			typedef TinyAD::Double<9> AReal;
#endif
			auto W2sFunc = [u, v, &Ibarinv, &Ibar, &lameAlpha, &lameBeta](Eigen::Vector<Real, 9> X, Eigen::Vector<Real, 9>* deriv = NULL, Eigen::Matrix<Real, 9, 9>* hess = NULL)
			{
				Eigen::Vector<AReal, 9> activeX = AReal::make_active(X);
				AReal a = (1 - u - v) * activeX(0) + u * activeX(1) + v * activeX(2);
				Eigen::Vector<AReal, 2> w = Ibarinv * activeX(Eigen::seq(3, 4));
				Real coeff = 4 * lameBeta * (lameAlpha + lameBeta) / (lameAlpha + 2 * lameBeta);
				Eigen::Matrix<AReal, 2, 2> II = activeX(Eigen::seq(5, 8)).reshaped(2, 2);
				AReal num = w.norm() * w.norm() * II.trace() - w.transpose() * II * w;
				AReal den = w.norm() * w.norm() * Ibar.trace() - w.transpose() * Ibar * w;
				AReal e = coeff * a * a * (num / den) * (num / den);
				if (deriv)
					(*deriv) = e.grad;
				if (hess)
					(*hess) = e.Hess;
				return TinyAD::to_passive(e);
			};
			Eigen::Vector<Real, 9> W2sX;
			W2sX(0) = sheet.getAmplitudes()(vidxi);
			W2sX(1) = sheet.getAmplitudes()(vidxj);
			W2sX(2) = sheet.getAmplitudes()(vidxk);
			W2sX(3) = dphi(0);
			W2sX(4) = dphi(1);
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					W2sX(i * 2 + j + 5) = II(i, j);
				}
			}
			Eigen::Vector<Real, 9> W2sderiv;
			Eigen::Matrix<Real, 9, 9> W2shess;
			w2se = W2sFunc(W2sX, NULL, NULL);
			energy += w2se;
		}

		// First term
		// Correct numerically
		Matrix22 M0 = Ibarinv * (I + 0.5 * a * a * dphidphiTensor + 0.5 * dadaTensor - Ibar);
		std::cout << "M0 \n"
		          << M0 << std::endl;
		std::cout << "W1 \n"
		          << w2se << std::endl;

		// // Third term
		// Numerically correct
		Matrix22 M2 = (Ibarinv * (dphidaTensor + dadphiTensor));
		energy += (a * a / 32.0) * SVNormSquared(sheet.getLameAlpha(), sheet.getLameBeta(), M2);
		std::cout << "M2 \n"
		          << M2 << std::endl;
		// // Fourth term
		// Correct numerically
		Matrix22 M3 = Ibarinv * dadaTensor;
		std::cout << "M3 \n"
		          << M3 << std::endl;
		std::cout << "source amps \n";
		for (int i = 0; i < 3; i++)
		{
			int vidx = faces[faceId].v[i];
			Real ai = sheet.getAmpVec()(vidx);
			std::cout << "a" << i << " " << ai << std::endl;
		}
		std::cout << "quad " << quadId << " weights u: " << u << " v: " << v << " amp " << a << std::endl;

		Real area = 0.5 * sqrt(Ibar.determinant());
		std::cout << "I \n"
		          << I << std::endl;
		std::cout << "Ibar \n"
		          << Ibar << std::endl;
		std::cout << "II \n"
		          << II << std::endl;
		std::cout << "dphi " << dphi.transpose() << std::endl;
		std::cout << "da \n"
		          << da.transpose() << std::endl;
		std::cout << "dphidphiT \n"
		          << dphidphiTensor << std::endl;
		std::cout << "dadphiT \n"
		          << dadphiTensor << std::endl;
		std::cout << "dadaT \n"
		          << dadaTensor << std::endl;
	}
};
