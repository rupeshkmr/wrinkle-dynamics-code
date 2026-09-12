#include <argus/common.hpp>
#include <argus/tfw_bending.hpp>
#include <argus/util.hpp>
#include <argus/geometry_functions.hpp>
#include "../../vendor/SecondFundamentalForm/SecondFundamentalFormDiscretization.h"
#include "halfedge/trimesh_types.h"

#define idxAt(i, j) 3 * i + j
#define mat(A, i, j) A(0, idxAt(i, j)), A(1, idxAt(i, j)), A(2, idxAt(i, j)), A(3, idxAt(i, j))
#define hessMat(A, i, j) A[0](i, j), A[1](i, j), A[2](i, j), A[3](i, j)

namespace argus
{
#ifdef ARGUS_CHECKPOINT
	void TFWBending::saveInfo()
	{
		nlohmann::json data;
		data["Energy"] = "TFW Bending";
		data["Is external"] = isExternal();
		simConf["Energies"].push_back(data);
	}
#endif

	Real TFWBending::getEnergyDensityFromQuadPerVtex1(int vid, int faceId, int quadId, Sheet& sheet, VectorN* deriv, MatrixNN* hess) const
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
			Eigen::Matrix<AReal3, 2, 2> M1, M2;
			M1 = Ibar.inverse() * (II - IIbar);
			AReal3 e1 = (0.5 * lA * (M1.trace() * M1.trace()) + lB * (M1 * M1).trace());
			energy += e1;
			M2 = Ibar.inverse() * (dphi * dphi.transpose());
			AReal3 e2 = (a * a * 0.5) * (0.5 * lA * (M2.trace() * M2.trace()) + lB * (M2 * M2).trace());
			energy += e2;

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

	Real TFWBending::getVertexEnergy1(Sheet& sheet, int vid, std::vector<int>& fids, VectorN* deriv, MatrixNN* hess) const
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

		Real thickness = sheet.getThickness();
		Real coeff = thickness * thickness * thickness;
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
				coeff = coeff * sheet.getQuadPoints()[quad].weight / 12.0;
				energies[i] += coeff * TFWBending::getEnergyDensityFromQuadPerVtex1(localVid, faceId, quad, sheet, deriv ? &tempDeriv : NULL, hess ? &tempHess : NULL);
				if (deriv)
				{
					tempDeriv = (coeff * tempDeriv).eval();
					derivatives[i](0) += tempDeriv(0);
					derivatives[i](1 + 2 * i) += tempDeriv(1);
					derivatives[i](1 + 2 * i + 1) += tempDeriv(2);
				}
				if (hess)
				{
					tempHess = (coeff * tempHess).eval();
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
	Real TFWBending::getEnergyDensityFromQuadPerVtex1(int vid, int faceId, int quadId, Sheet& sheet, Real* deriv, Real* hess) const
	{
		Real energy = 0;
		const std::vector<trimesh::triangle_t>& faces = sheet.getTriangles();
		const std::vector<trimesh::edge_t>& edges = sheet.getEdges();
		const Matrix3N& positions = sheet.getPositionsN3();
		const Eigen::MatrixXi& edgeOppositeVertices = sheet.getEdgeOppositeVertices();
		const std::vector<trimesh::triangle_t>& vertexOppositeEdges = sheet.getVertexOppositeEdgesPerFace();
		int nfaces = faces.size();
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
		// Eigen::Matrix<Real, 4, 9> gradI;
		// Eigen::Matrix<Real, 4, 18> gradII; // gradI(i) contains dI(i)/dX
		// std::array<Eigen::Matrix<Real, 9, 9>, 4> hessI;
		// std::array<Eigen::Matrix<Real, 18, 18>, 4> hessII;

		std::shared_ptr<SecondFundamentalFormDiscretization> sff;
		sff = std::make_shared<MidedgeAverageFormulation>();
		Matrix22 II = sff->secondFundamentalFormNew(sheet.getBaseMesh(), sheet.getPositionsN3(), sheet.getEdgeDofs(), faceId, NULL, NULL);
		Matrix22 Ibar = sheet.getRestShape().getIbars()[faceId];
		Matrix22 Ibarinv = Ibar.inverse();
		Matrix22 IIbar = sheet.getRestShape().getIIbars()[faceId];
		Vector3 gradA; // datilda/da
		Vector2 da;
		// Matrix23 gradDA; // dda/da
		Real u = sheet.getQuadPoints()[quadId].u;
		Real v = sheet.getQuadPoints()[quadId].v;

		Real a = sheet.computeAmplitudesFromQuad(faceId, quadId, &da, &gradA, NULL, nullptr, nullptr);

		Matrix22 dphidphiTensor;
		dphidphiTensor = sheet.computeDphiDphi1VectorPerFaceTensorNew(faceId, quadId, NULL, NULL);

		Real lameAlpha = sheet.getLameAlpha();
		Real lameBeta = sheet.getLameBeta();
		if (deriv)
			(*deriv) = 0;
		if (hess)
			(*hess) = 0;
		// First term
		// Correct numerically
		Matrix22 M0 = Ibarinv * (II - IIbar);
		Matrix22 dMSVdM;
		Real MSV;
		MSV = SVNormSquared(lameAlpha, lameBeta, M0, &dMSVdM);
		energy += MSV;

		// Second term
		Matrix22 M1 = Ibarinv * dphidphiTensor;
		Real M1norm = SVNormSquared(lameAlpha, lameBeta, M1, NULL);
		energy += a * a * 0.5 * M1norm;
		if (deriv)
		{
			(*deriv) += a * gradA(vid) * M1norm;
		}
		if (hess)
		{
			(*hess) += gradA(vid) * gradA(vid) * M1norm;
		}
		Real area = 0.5 * sqrt(Ibar.determinant());
		if (deriv)
			(*deriv) = area * (*deriv);
		if (hess)
			(*hess) = area * (*hess);
		if ((*hess) < 0)
			(*hess) *= -1;
		return area * energy;
	}

	Real TFWBending::getVertexEnergy1New(Sheet& sheet, int vid, Real* deriv, Real* hess) const
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
			(*deriv) = 0;
		if (hess)
			(*hess) = 0;

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
				std::cout << "Vertex not found! Face " << faceId << ". Error in tfw_bending.cpp\n";
				throw "Error!";
			}
			Real thickness = sheet.getThickness();
			Real coeff = thickness * thickness * thickness / 12.0;
			for (int quad = 0; quad < sheet.getQuadPoints().size(); quad++)
			{
				Real tempDeriv, tempHess;
				energy += sheet.getQuadPoints()[quad].weight * coeff * TFWBending::getEnergyDensityFromQuadPerVtex1(localVid, faceId, quad, sheet, deriv ? &tempDeriv : NULL, hess ? &tempHess : NULL);
				if (deriv)
					(*deriv) += sheet.getQuadPoints()[quad].weight * coeff * tempDeriv;
				if (hess)
					(*hess) += sheet.getQuadPoints()[quad].weight * coeff * tempHess;
			}
		}
		return energy;
	}
	Real TFWBending::getVertexEnergy1(Sheet& sheet, int vid, Real* deriv, Real* hess) const
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

			Vector3 gradA;
			Real lameAlpha = sheet.getLameAlpha();
			Real lameBeta = sheet.getLameBeta();
			Matrix22 II;
			Matrix22 IIbar = sheet.getRestShape().getIIbars()[faceId];
			if (sheet.cachedIIsValid())
			{
				II.noalias() = sheet.getCachedIIs().second[faceId];
			}
			else
			{
				// std::shared_ptr<SecondFundamentalFormDiscretization> sff;
				// sff = std::make_shared<MidedgeAverageFormulation>();
				II.noalias() = sffobject->secondFundamentalFormNew(sheet.getBaseMesh(), sheet.getPositionsN3(), sheet.getEdgeDofs(), faceId, NULL, NULL);
			}
			// First term
			// Correct numerically
			Matrix22 M0 = Ibarinv * (II - IIbar);
			Real e1 = SVNormSquared(lameAlpha, lameBeta, M0, NULL);
			Real thickness = sheet.getThickness();
			Real coeff = thickness * thickness * thickness / 12.0;
			Real area = 0.5 * sqrt(Ibar.determinant());
			Matrix22 M1 = Ibarinv * dphidphiTensor;
			Real M1norm = SVNormSquared(lameAlpha, lameBeta, M1, NULL);
			for (int quadId = 0; quadId < sheet.getQuadPoints().size(); quadId++)
			{
				Real energy = 0;
				Real tempDeriv = 0;
				Real tempHess = 0;
				Real a = sheet.computeAmplitudesFromQuad(faceId, quadId, NULL, &gradA, NULL, nullptr, nullptr);
				// First term
				energy += e1;
				// Second term

				energy += a * a * 0.5 * M1norm;
				if (deriv)
				{
					tempDeriv += a * gradA(localVid) * M1norm;
				}
				if (hess)
				{
					tempHess += gradA(localVid) * gradA(localVid) * M1norm;
				}

				totalEnergy += sheet.getQuadPoints()[quadId].weight * energy * coeff * area;
				if (deriv)
					(*deriv) += sheet.getQuadPoints()[quadId].weight * tempDeriv * coeff * area;
				if (hess)
					(*hess) += sheet.getQuadPoints()[quadId].weight * tempHess * coeff * area;
			}
		}
		// if (*hess < 0)
		// (*hess) = -1 * (*hess);
		// consolidate
		return totalEnergy;
	}

	Real TFWBending::getEnergyDensityFromQuadPerVtex4(int vid, int faceId, int quadId, Sheet& sheet, Vector4* deriv, Matrix44* hess) const
	{
		Real energy = 0;
		const std::vector<trimesh::triangle_t>& faces = sheet.getTriangles();
		const std::vector<trimesh::edge_t>& edges = sheet.getEdges();
		const Matrix3N positions = sheet.getPositionsN3();
		const Eigen::MatrixXi edgeOppositeVertices = sheet.getEdgeOppositeVertices();
		const std::vector<trimesh::triangle_t> vertexOppositeEdges = sheet.getVertexOppositeEdgesPerFace();
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
		std::vector<Matrix33> hessDaDa;
		Matrix22 dadaTensor;
		dadaTensor = sheet.computeDaDaTensor(faceId, quadId, &gradDaDaTensor, &hessDaDa);
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
		Matrix22 M0 = Ibarinv * (II - IIbar);

		Matrix22 dMSVdM;
		Real MSV;
		MSV = SVNormSquared(lameAlpha, lameBeta, M0, &dMSVdM);
		energy += MSV;
		// Derivatives wrt face positions
		if (deriv)
			for (int i = 0; i < 3; i++)
			{
				Matrix22 dIIdxij;
				dIIdxij << mat(gradII, vid, i);
				Matrix22 dM0dxij = Ibarinv * dIIdxij;
				(*deriv)(i + 1) += doubleContraction(lameAlpha * M0.trace() * Matrix22::Identity() + 2 * lameBeta * M0, dM0dxij);
			}
		if (hess)
			for (unsigned int j = 0; j < 3; j++)
			{
				Matrix22 dIIdxij;
				dIIdxij << mat(gradII, vid, j);
				for (unsigned int k = 0; k <= j; k++)
				{
					Matrix22 d2MSVdxijdM, dMdxik, d2Mdxijdxik;
					// d2MSVdxijdM
					d2MSVdxijdM = lameAlpha * doubleContraction(Ibarinv * dIIdxij, Matrix22::Identity()) * Matrix22::Identity() + 2 * lameBeta * Ibarinv * dIIdxij;
					// dMdxik
					Matrix22 dIIdxik;
					dIIdxik << mat(gradII, vid, k);
					dMdxik = Ibarinv * dIIdxik;
					// d2Mdxijdxik
					Matrix22 hessIIjk;
					hessIIjk << hessMat(hessII, idxAt(vid, j), idxAt(vid, k));
					d2Mdxijdxik = Ibarinv * hessIIjk;
					Real d2Edxijdxik = doubleContraction(d2MSVdxijdM, dMdxik)
					    + doubleContraction(dMSVdM, d2Mdxijdxik);
					(*hess)(j + 1, k + 1) += d2Edxijdxik;
					if (j != k)
						(*hess)(k + 1, j + 1) += d2Edxijdxik;
				}
			}

		// Second term
		Matrix22 M1 = Ibarinv * dphidphiTensor;
		Real M1norm = SVNormSquared(lameAlpha, lameBeta, M1, NULL);
		energy += a * a * 0.5 * M1norm;
		if (deriv)
		{
			(*deriv)(0) += a * gradA(vid) * M1norm;
		}
		if (hess)
		{
			(*hess)(0, 0) += gradA(vid) * gradA(vid) * M1norm;
		}
		// project
		// if (hess)
		// 	*hess = argus::util::hessianProjection((*hess), true);
		Real area = 0.5 * sqrt(Ibar.determinant());
		if (deriv)
			*deriv = (area * (*deriv)).eval();
		if (hess)
			(*hess) = (area * (*hess)).eval();

		return area * energy;
	}

	Real TFWBending::getVertexEnergy4New(Sheet& sheet, int vid, Vector4* deriv, Matrix44* hess) const
	{
		// need to compute deriv and hessian wrt the vertex at index vid
		// just for testing
		Real energy = 0;
		Eigen::VectorXi faces = sheet.sheetInfo.vertexFacesSet.row(vid);
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
				std::cout << "Vertex not found! Face " << faceId << ". Error in tfw_stretching.cpp\n";
				exit(0);
			}
			Real thickness = sheet.getThickness();
			Real coeff = thickness * thickness * thickness / 12.0;

			for (int quad = 0; quad < sheet.getQuadPoints().size(); quad++)
			{
				Vector4 tempDeriv;
				Matrix44 tempHess;
				energies[i] += sheet.getQuadPoints()[quad].weight * coeff * TFWBending::getEnergyDensityFromQuadPerVtex4(localVid, faceId, quad, sheet, deriv ? &tempDeriv : NULL, hess ? &tempHess : NULL);
				if (deriv)
					derivatives[i] += sheet.getQuadPoints()[quad].weight * coeff * tempDeriv;
				if (hess)
					hessians[i] += sheet.getQuadPoints()[quad].weight * coeff * tempHess;
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
	Real TFWBending::getVertexEnergy4(Sheet& sheet, int vid, Vector4* deriv, Matrix44* hess) const
	{
		// need to compute deriv and hessian wrt the vertex at index vid
		// just for testing
		Real lameAlpha = sheet.getLameAlpha();
		Real lameBeta = sheet.getLameBeta();

		// need to compute deriv and hessian wrt the vertex at index vid
		// just for testing
		Real totalEnergy = 0;
		Eigen::Vector<int, 32> faces = sheet.sheetInfo.vertexFacesSet.row(vid);
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
				std::cerr << "Vertex not found! Error in tfw_bending.cpp\n";
				throw "Error!";
			}
			Matrix22 Ibar = sheet.getRestShape().getIbars()[faceId];
			Matrix22 Ibarinv = Ibar.inverse();
			Matrix22 IIbar = sheet.getRestShape().getIIbars()[faceId];
			Eigen::Matrix<Real, 4, 18> gradII; // gradI(i) contains dI(i)/dX
			std::array<Eigen::Matrix<Real, 18, 18>, 4> hessII;
			Matrix22 II = sffobject->secondFundamentalFormNew(sheet.getBaseMesh(), sheet.getPositionsN3(), sheet.getEdgeDofs(), faceId, &gradII, &hessII);
			Matrix22 dphidphiTensor;
			if (sheet.cachedDphiDphiTValid())
			{
				dphidphiTensor.noalias() = sheet.getCachedDphiDpiT().second[faceId];
			}
			else
			{
				dphidphiTensor.noalias() = sheet.computeDphiDphi1VectorPerFaceTensor(faceId, 0, NULL, NULL);
			}
			Real area = 0.5 * sqrt(Ibar.determinant());
			// First term
			// Correct numerically
			Matrix22 M0 = Ibarinv * (II - IIbar);

			Matrix22 dMSVdM;
			Real e0 = SVNormSquared(lameAlpha, lameBeta, M0, &dMSVdM);

			Vector3 d0;
			d0.setZero();
			Matrix33 h0;
			h0.setZero();
			// Derivatives wrt face positions
			if (deriv)
				for (int i = 0; i < 3; i++)
				{
					Matrix22 dIIdxij;
					dIIdxij << mat(gradII, localVid, i);
					Matrix22 dM0dxij = Ibarinv * dIIdxij;
					d0(i) = doubleContraction(lameAlpha * M0.trace() * Matrix22::Identity() + 2 * lameBeta * M0, dM0dxij);
				}
			if (hess)
				for (unsigned int j = 0; j < 3; j++)
				{
					Matrix22 dIIdxij;
					dIIdxij << mat(gradII, localVid, j);
					// d2MSVdxijdM
					Matrix22 d2MSVdxijdM = lameAlpha * doubleContraction(Ibarinv * dIIdxij, Matrix22::Identity()) * Matrix22::Identity() + 2 * lameBeta * Ibarinv * dIIdxij;
					for (unsigned int k = 0; k <= j; k++)
					{
						Matrix22 dMdxik, d2Mdxijdxik;
						// dMdxik
						Matrix22 dIIdxik;
						dIIdxik << mat(gradII, localVid, k);
						dMdxik = Ibarinv * dIIdxik;
						// d2Mdxijdxik
						Matrix22 hessIIjk;
						hessIIjk << hessMat(hessII, idxAt(localVid, j), idxAt(localVid, k));
						d2Mdxijdxik = Ibarinv * hessIIjk;
						Real d2Edxijdxik = doubleContraction(d2MSVdxijdM, dMdxik)
						    + doubleContraction(dMSVdM, d2Mdxijdxik);
						h0(j, k) += d2Edxijdxik;
						if (j != k)
							h0(k, j) += d2Edxijdxik;
					}
				}
			Real thickness = sheet.getThickness();
			Real coeff = area * thickness * thickness * thickness / 12.0; // std::pow(sheet.getThickness(), 3) / 12.0;
			Matrix22 tempM = Ibarinv * dphidphiTensor;
			Real M1svnorm = (SVNormSquared(lameAlpha, lameBeta, tempM));

			for (int quadId = 0; quadId < sheet.getQuadPoints().size(); quadId++)
			{
				Vector4 tempDeriv;
				tempDeriv.setZero();
				Matrix44 tempHess;
				tempHess.setZero();
				Real energy = 0;
				Vector3 gradA;
				Real a = sheet.computeAmplitudesFromQuadNew(faceId, quadId, NULL, &gradA, NULL, nullptr, nullptr);

				energy = e0;
				if (deriv)
				{
					tempDeriv.segment<3>(1) = d0;
				}
				if (hess)
					tempHess.block(1, 1, 3, 3) = h0;
				// Second term
				energy += a * a * 0.5 * M1svnorm;
				if (deriv)
				{
					tempDeriv(0) += gradA(localVid) * a * M1svnorm;
				}
				if (hess)
				{
					tempHess(0, 0) += gradA(localVid) * gradA(localVid) * M1svnorm;
				}

				if (deriv)
					(*deriv) += coeff * tempDeriv * sheet.getQuadPoints()[quadId].weight;
				if (hess)
					(*hess) += coeff * tempHess * sheet.getQuadPoints()[quadId].weight;
				totalEnergy += coeff * energy * sheet.getQuadPoints()[quadId].weight;
			}
		}
		// if (hess)
		// 	*hess = argus::util::hessianProjection(*hess, true);

		return totalEnergy;
	}
	Real TFWBending::getEnergyDensityFromQuadPerVtex(int vid, int faceId, int quadId, Sheet& sheet, Vector3* deriv, Matrix33* hess) const
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
		Eigen::Matrix<Real, 4, 18> gradII; // gradI(i) contains dI(i)/dX
		std::array<Eigen::Matrix<Real, 18, 18>, 4> hessII;
		std::shared_ptr<SecondFundamentalFormDiscretization> sff;
		sff = std::make_shared<MidedgeAverageFormulation>();
		Matrix22 II = sff->secondFundamentalFormNew(sheet.getBaseMesh(), sheet.getPositionsN3(), sheet.getEdgeDofs(), faceId, &gradII, &hessII);
		Matrix22 Ibar = sheet.getRestShape().getIbars()[faceId];
		Matrix22 Ibarinv = Ibar.inverse();
		Matrix22 IIbar = sheet.getRestShape().getIIbars()[faceId];
		Real u = sheet.getQuadPoints()[quadId].u;
		Real v = sheet.getQuadPoints()[quadId].v;

		Real a = sheet.computeAmplitudesFromQuad(faceId, quadId, NULL, NULL, NULL, nullptr, nullptr);
		Matrix22 gradDphi;
		Vector2 dphi = sheet.computeDphi1VectorPerFace(faceId, &gradDphi);
		Matrix22 dphidphiTensor = sheet.computeDphiDphi1VectorPerFaceTensorNew(faceId, quadId, NULL, NULL);
		Matrix22 dadaTensor = sheet.computeDaDaTensor(faceId, quadId, NULL, NULL);
		Matrix22 dadphiTensor = sheet.computeDaDphi1VectorPerFaceTensor(faceId, quadId, NULL, NULL);
		Matrix22 dphidaTensor = sheet.computeDphi1VectorPerFaceDaTensor(faceId, quadId, NULL, NULL);
		Real lameAlpha = sheet.getLameAlpha();
		Real lameBeta = sheet.getLameBeta();
		if (deriv)
			deriv->setZero();
		if (hess)
			hess->setZero();

		// First term
		// Correct numerically
		Matrix22 M0 = Ibarinv * (II - IIbar);
		Matrix22 dMSVdM;
		Real MSV;
		MSV = SVNormSquared(lameAlpha, lameBeta, M0, &dMSVdM);
		energy += MSV;
		// Derivatives wrt face positions
		if (deriv)
			for (int i = 0; i < 3; i++)
			{
				Matrix22 dIIdxij;
				dIIdxij << mat(gradII, vid, i);
				Matrix22 dM0dxij = Ibarinv * dIIdxij;
				(*deriv)(i) += doubleContraction(lameAlpha * M0.trace() * Matrix22::Identity() + 2 * lameBeta * M0, dM0dxij);
			}
		if (hess)
			for (unsigned int j = 0; j < 3; j++)
			{
				Matrix22 dIIdxij;
				dIIdxij << mat(gradII, vid, j);
				for (unsigned int k = 0; k <= j; k++)
				{
					Matrix22 d2MSVdxijdM, dMdxik, d2Mdxijdxik;
					// d2MSVdxijdM
					d2MSVdxijdM = lameAlpha * doubleContraction(Ibarinv * dIIdxij, Matrix22::Identity()) * Matrix22::Identity() + 2 * lameBeta * Ibarinv * dIIdxij;
					// dMdxik
					Matrix22 dIIdxik;
					dIIdxik << mat(gradII, vid, k);
					dMdxik = Ibarinv * dIIdxik;
					// d2Mdxijdxik
					Matrix22 hessIIjk;
					hessIIjk << hessMat(hessII, idxAt(vid, j), idxAt(vid, k));
					d2Mdxijdxik = Ibarinv * hessIIjk;
					Real d2Edxijdxik = doubleContraction(d2MSVdxijdM, dMdxik)
					    + doubleContraction(dMSVdM, d2Mdxijdxik);
					(*hess)(j, k) += d2Edxijdxik;
					if (j != k)
						(*hess)(k, j) += d2Edxijdxik;
				}
			}

		// Second term
		Matrix22 M1 = Ibarinv * dphidphiTensor;
		energy += a * a * 0.5 * (SVNormSquared(lameAlpha, lameBeta, M1));

		Real area = 0.5 * sqrt(Ibar.determinant());
		if (deriv)
			*deriv = (area * (*deriv)).eval();
		if (hess)
			(*hess) = (area * (*hess)).eval();
		// project
		// if (hess)
		// 	(*hess) = argus::util::hessianProjection((*hess), true);
		return area * energy;
	}

	Real TFWBending::getVertexEnergyNew(Sheet& sheet, int vid, Vector3* deriv, Matrix33* hess) const
	{
		// need to compute deriv and hessian wrt the vertex at index vid
		// just for testing
		Real energy = 0;
		Eigen::Vector<int, 32> faces = sheet.sheetInfo.vertexFacesSet.row(vid);
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
		Real thickness = sheet.getThickness();
		Real coeff = thickness * thickness * thickness / 12.0; // std::pow(sheet.getThickness(), 3) / 12.0;
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
				std::cerr << "Vertex not found! Error in tfw_bending.cpp\n";
				throw "Error!";
			}
			for (int quad = 0; quad < sheet.getQuadPoints().size(); quad++)
			{
				Vector3 tempDeriv1, tempDeriv2, tempDeriv;
				Matrix33 tempHess, tempHess1, tempHess2;
				energy += sheet.getQuadPoints()[quad].weight * coeff * TFWBending::getEnergyDensityFromQuadPerVtex(localVid, faceId, quad, sheet, deriv ? &tempDeriv : NULL, hess ? &tempHess : NULL);
				if (deriv)
					(*deriv) += sheet.getQuadPoints()[quad].weight * coeff * tempDeriv;
				if (hess)
					(*hess) += sheet.getQuadPoints()[quad].weight * coeff * tempHess;
			}
		}
		// if (hess)
		// 	*hess = argus::util::hessianProjection(*hess, true);
		return energy;
	}
	Real TFWBending::getVertexEnergy(Sheet& sheet, int vid, Vector3* deriv, Matrix33* hess) const
	{
		// need to compute deriv and hessian wrt the vertex at index vid
		// just for testing
		Real lameAlpha = sheet.getLameAlpha();
		Real lameBeta = sheet.getLameBeta();
		// need to compute deriv and hessian wrt the vertex at index vid
		// just for testing
		Real totalEnergy = 0;
		Eigen::Vector<int, 32> faces = sheet.sheetInfo.vertexFacesSet.row(vid);
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
		Real thickness = sheet.getThickness();
		Real coeff = thickness * thickness * thickness / 12.0; // std::pow(sheet.getThickness(), 3) / 12.0;
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
				std::cerr << "Vertex not found! Error in tfw_bending.cpp\n";
				throw "Error!";
			}
			Matrix22 Ibar = sheet.getRestShape().getIbars()[faceId];
			Matrix22 Ibarinv = Ibar.inverse();
			Matrix22 IIbar = sheet.getRestShape().getIIbars()[faceId];
			Eigen::Matrix<Real, 4, 18> gradII; // gradI(i) contains dI(i)/dX
			std::array<Eigen::Matrix<Real, 18, 18>, 4> hessII;
			Matrix22 II = sffobject->secondFundamentalFormNew(sheet.getBaseMesh(), sheet.getPositionsN3(), sheet.getEdgeDofs(), faceId, &gradII, &hessII);
			Matrix22 dphidphiTensor;
			if (sheet.cachedDphiDphiTValid())
			{
				dphidphiTensor.noalias() = sheet.getCachedDphiDpiT().second[faceId];
			}
			else
			{
				dphidphiTensor.noalias() = sheet.computeDphiDphi1VectorPerFaceTensor(faceId, 0, NULL, NULL);
			}
			Real area = 0.5 * sqrt(Ibar.determinant());
			// First term
			// Correct numerically
			Matrix22 M0 = Ibarinv * (II - IIbar);
			Matrix22 dMSVdM;
			Real e0 = SVNormSquared(lameAlpha, lameBeta, M0, &dMSVdM);
			Vector3 d0;
			d0.setZero();
			Matrix33 h0;
			h0.setZero();
			// Derivatives wrt face positions
			if (deriv)
				for (int i = 0; i < 3; i++)
				{
					Matrix22 dIIdxij;
					dIIdxij << mat(gradII, localVid, i);
					Matrix22 dM0dxij = Ibarinv * dIIdxij;
					d0(i) = doubleContraction(lameAlpha * M0.trace() * Matrix22::Identity() + 2 * lameBeta * M0, dM0dxij);
				}
			if (hess)
				for (unsigned int j = 0; j < 3; j++)
				{
					Matrix22 dIIdxij;
					dIIdxij << mat(gradII, localVid, j);
					// d2MSVdxijdM
					Matrix22 d2MSVdxijdM = lameAlpha * doubleContraction(Ibarinv * dIIdxij, Matrix22::Identity()) * Matrix22::Identity() + 2 * lameBeta * Ibarinv * dIIdxij;
					for (unsigned int k = 0; k <= j; k++)
					{
						Matrix22 dMdxik, d2Mdxijdxik;
						// dMdxik
						Matrix22 dIIdxik;
						dIIdxik << mat(gradII, localVid, k);
						dMdxik = Ibarinv * dIIdxik;
						// d2Mdxijdxik
						Matrix22 hessIIjk;
						hessIIjk << hessMat(hessII, idxAt(localVid, j), idxAt(localVid, k));
						d2Mdxijdxik = Ibarinv * hessIIjk;
						Real d2Edxijdxik = doubleContraction(d2MSVdxijdM, dMdxik)
						    + doubleContraction(dMSVdM, d2Mdxijdxik);
						h0(j, k) += d2Edxijdxik;
						if (j != k)
							h0(k, j) += d2Edxijdxik;
					}
				}
			Real thickness = sheet.getThickness();
			Real coeff = area * thickness * thickness * thickness / 12.0; // std::pow(sheet.getThickness(), 3) / 12.0;
			Matrix22 tempM = Ibarinv * dphidphiTensor;
			Real M1svnorm = (SVNormSquared(lameAlpha, lameBeta, tempM));
			for (int quadId = 0; quadId < sheet.getQuadPoints().size(); quadId++)
			{
				Vector3 tempDeriv;
				tempDeriv.setZero();
				Matrix33 tempHess;
				tempHess.setZero();
				Real energy = 0;
				Real a = sheet.computeAmplitudesFromQuadNew(faceId, quadId, NULL, NULL, NULL, nullptr, nullptr);

				energy += e0;
				if (deriv)
				{
					tempDeriv = d0;
				}
				if (hess)
					tempHess = h0;
				// Second term
				energy += a * a * 0.5 * M1svnorm;

				if (deriv)
					(*deriv) += coeff * tempDeriv * sheet.getQuadPoints()[quadId].weight;
				if (hess)
					(*hess) += coeff * tempHess * sheet.getQuadPoints()[quadId].weight;
				totalEnergy += coeff * energy * sheet.getQuadPoints()[quadId].weight;
			}
		}
		// if (hess)
		// 	*hess = argus::util::hessianProjection(*hess, true);
		return totalEnergy;
	}

	Real getBendingEnergyDensityFromQuadAmpsPerVtexDphis1VectorPerFace(int faceId, int quadId, Sheet& sheet, VectorN* deriv, MatrixNN* hess)
	{
		const std::vector<trimesh::triangle_t>& faces = sheet.getTriangles();
		const std::vector<trimesh::edge_t>& edges = sheet.getEdges();
		const Matrix3N& positions = sheet.getPositionsN3();
		const Eigen::MatrixXi& edgeOppositeVertices = sheet.getEdgeOppositeVertices();
		const std::vector<trimesh::triangle_t>& vertexOppositeEdges = sheet.getVertexOppositeEdgesPerFace();
		int nfaces = faces.size();
		int nedges = edges.size();
		int nvertices = positions.cols();

		std::shared_ptr<SecondFundamentalFormDiscretization> sff;
		sff = std::make_shared<MidedgeAverageFormulation>();
		MatrixNN gradII;
		std::vector<Eigen::MatrixXd> hessII;
		Matrix22 II = sff->secondFundamentalForm(sheet.getBaseMesh(), sheet.getPositionsN3().transpose(), sheet.getEdgeDofs(), faceId, deriv || hess ? &gradII : NULL, hess ? &hessII : NULL);
		Matrix22 Ibar = sheet.getRestShape().getIbars()[faceId];
		Matrix22 Ibarinv = Ibar.inverse();
		Matrix22 IIbar = sheet.getRestShape().getIIbars()[faceId];

		Vector3 gradA; // datilda/da
		Vector2 da;
		Matrix23 gradDA; // dda/da
		Real a = sheet.computeAmplitudesFromQuad(faceId, quadId, &da, &gradA, &gradDA, nullptr, nullptr);
		std::vector<Matrix22> gradDphiDphiTensor, hessDphiDphiTensor;

		Matrix22 dphidphiTensor = sheet.computeDphiDphi1VectorPerFaceTensor(faceId, quadId, deriv || hess ? &gradDphiDphiTensor : NULL, hess ? &hessDphiDphiTensor : NULL);
		if (deriv)
		{
			deriv->resize(3 + 2 + 3 * 3 + 3 * 3);
			deriv->setZero();
		}
		if (hess)
		{
			hess->resize(3 + 2 + 3 * 3 + 3 * 3, 23);
			hess->setZero();
		}
		Real energy = 0.0;
		Real lameAlpha = sheet.getLameAlpha();
		Real lameBeta = sheet.getLameBeta();
		// First term
		// Numerically correct
		Matrix22 M0 = Ibar.inverse() * (II - IIbar);
		energy += SVNormSquared(sheet.getLameAlpha(), sheet.getLameBeta(), M0);
		if (deriv)
		{
			for (int i = 0; i < 9; i++)
			{
				// For current triangles vertices
				Matrix22 dM0dxi;
				dM0dxi << gradII(0, i), gradII(1, i), gradII(2, i), gradII(3, i);
				dM0dxi = (Ibar.inverse() * dM0dxi).eval();
				deriv->coeffRef(5 + i) += 2 * (dM0dxi * M0).trace() * lameBeta + dM0dxi.trace() * M0.trace() * lameAlpha;

				// For opposite vertices
				Matrix22 dM0dxitilda;
				dM0dxitilda << gradII(0, 9 + i), gradII(1, 9 + i), gradII(2, 9 + i), gradII(3, 9 + i);
				dM0dxitilda = (Ibar.inverse() * dM0dxitilda).eval();
				deriv->coeffRef(14 + i) += 2 * (dM0dxitilda * M0).trace() * lameBeta + dM0dxitilda.trace() * M0.trace() * lameAlpha;
			}
		}
		if (hess)
		{
			Matrix22 T1, T2, T3, T4;
			// pos vs pos
			for (int i = 0; i < 18; i++)
			{
				for (int j = 0; j < 18; j++)
				{
					Matrix22 hessIIij, gradIIi, gradIIj;
					hessIIij << hessII.at(0)(i, j), hessII.at(1)(i, j), hessII.at(2)(i, j), hessII.at(3)(i, j);
					gradIIi << gradII(0, i), gradII(1, i), gradII(2, i), gradII(3, i);
					gradIIj << gradII(0, j), gradII(1, j), gradII(2, j), gradII(3, j);
					T1 = Ibarinv * hessIIij;
					T2 = lameAlpha * M0.trace() * Matrix22::Identity() + 2 * lameBeta * M0;
					T3 = Ibarinv * gradIIj;
					T4 = lameAlpha * (doubleContraction(Ibarinv, gradIIi)) * Matrix22::Identity() + 2 * lameBeta * Ibarinv * gradIIi;
					hess->coeffRef(5 + i, 5 + j) += doubleContraction(T1, T2) + doubleContraction(T3, T4);
				}
			}
		}

		// Second term
		// Correct numerically
		Matrix22 M1 = Ibar.inverse() * dphidphiTensor;
		Real SV = SVNormSquared(sheet.getLameAlpha(), sheet.getLameBeta(), M1);
		energy += (a * a / 2.0) * SV;
		if (deriv)
		{
			// Derivatives wrt face amplitudes
			for (int i = 0; i < 3; i++)
			{
				deriv->coeffRef(i) += a * gradA(i) * SV;
			}
			// // Derivatives wrt face dphis
			for (int i = 0; i < 2; i++)
			{
				Matrix22 dM1dDphii = Ibar.inverse() * (gradDphiDphiTensor[i]);
				deriv->coeffRef(3 + i) += 0.5 * a * a * (2 * (dM1dDphii * M1).trace() * lameBeta + dM1dDphii.trace() * M1.trace() * lameAlpha);
				// Matrix22 dM1dphii = Ibar.inverse() * gradDphiDphiTensor[i];
				// deriv->coeffRef(3+i) += (a*a/2)* (lameAlpha * M1.trace() * dM1dphii.trace() + 2 * lameBeta * (M1.transpose() * dM1dphii).trace());
			}
		}
		if (hess)
		{
			Matrix22 T1, T2, T3, T4, N;
			// amps vs amps
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					T3 = gradA(j) * Ibarinv * dphidphiTensor;
					T4 = lameAlpha * gradA(i) * doubleContraction(Ibarinv, dphidphiTensor) * Matrix22::Identity()
					    + 2 * lameBeta * gradA(i) * Ibarinv * dphidphiTensor;
					hess->coeffRef(i, j) += doubleContraction(T3, T4);
				}
			}
			N = Ibarinv * dphidphiTensor;
			// dphis vs amps
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					T1 = gradA(j) * Ibarinv * gradDphiDphiTensor.at(i);
					T2 = lameAlpha * a * N.trace() * Matrix22::Identity() + 2 * lameBeta * a * N;
					T3 = gradA(j) * Ibarinv * dphidphiTensor;
					T4 = lameAlpha * a * doubleContraction(Ibarinv, gradDphiDphiTensor.at(i).transpose()) * Matrix22::Identity()
					    + 2 * lameBeta * a * Ibarinv * gradDphiDphiTensor.at(i);
					hess->coeffRef(3 + i, j) += doubleContraction(T1, T2) + doubleContraction(T3, T4);
					hess->coeffRef(j, 3 + i) += doubleContraction(T1, T2) + doubleContraction(T3, T4);
				}
			}
			// dphis vs dphis
			T2 = lameAlpha * a * N.trace() * Matrix22::Identity() + 2 * lameBeta * a * N;
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					Matrix22 hessDphiDphiij;
					hessDphiDphiij << hessDphiDphiTensor.at(0)(i, j), hessDphiDphiTensor.at(1)(i, j),
					    hessDphiDphiTensor.at(2)(i, j), hessDphiDphiTensor.at(3)(i, j);
					T1 = a * Ibarinv * hessDphiDphiij;
					T3 = a * Ibarinv * gradDphiDphiTensor.at(j);
					T4 = lameAlpha * (doubleContraction(Ibarinv, gradDphiDphiTensor.at(i).transpose())) * Matrix22::Identity()
					    + lameBeta * 2 * a * gradDphiDphiTensor.at(i);
					hess->coeffRef(3 + i, 3 + j) += doubleContraction(T1, T2) + doubleContraction(T3, T4);
				}
			}
		}

		Real area = 0.5 * sqrt(Ibar.determinant());
		if (deriv)
			(*deriv) = area * (*deriv);
		if (hess)
		{
			(*hess) = area * (*hess);
			// project hessian
			// (*hess) = argus::util::lowRankApprox(*hess);
		}

		return area * energy;
	}

	Real getBendingEnergyAmpsPerVtexDphis1VectorPerFace(Sheet& sheet, VectorN* deriv, std::vector<Eigen::Triplet<Real>>* hess)
	{
		// iterate through all the faces in the mesh
		// iterate through each quad point and sum up the energies
		// We can parallelize energy computation for each face, thus, we want to create an array per face for energy and derivatives
		std::vector<trimesh::triangle_t> edgeIDOppVertexIDPerface = sheet.getVertexOppositeEdgesPerFace();
		std::vector<trimesh::triangle_t> faces = sheet.getTriangles();
		Eigen::MatrixXi vtexIDOppositeEdge = sheet.getEdgeOppositeVertices();
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

#ifdef USE_OMP
#pragma omp parallel for
#endif
		for (int i = 0; i < nfaces; i++)
		{
			// if(tensionFaces.find(i) != tensionFaces.end())
			//     continue;
			// for each face iterate through each quad to get the energies
			for (int quad = 0; quad < sheet.getQuadPoints().size(); quad++)
			{
				if (hess)
				{
					VectorN tempDeriv;
					MatrixNN tempHess;
					energies[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * sheet.getThickness() * sheet.getThickness() * getBendingEnergyDensityFromQuadAmpsPerVtexDphis1VectorPerFace(i, quad, sheet, &tempDeriv, &tempHess) / 12.0;
					derivatives[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * sheet.getThickness() * sheet.getThickness() * tempDeriv / 12.0;
					hessians[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * sheet.getThickness() * sheet.getThickness() * tempHess / 12.0;
				}
				else if (deriv)
				{
					VectorN tempDeriv;
					energies[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * sheet.getThickness() * sheet.getThickness() * getBendingEnergyDensityFromQuadAmpsPerVtexDphis1VectorPerFace(i, quad, sheet, &tempDeriv, nullptr) / 12.0;
					derivatives[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * sheet.getThickness() * sheet.getThickness() * tempDeriv / 12.0;
				}
				else
				{
					energies[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * sheet.getThickness() * sheet.getThickness() * getBendingEnergyDensityFromQuadAmpsPerVtexDphis1VectorPerFace(i, quad, sheet, nullptr, nullptr) / 12.0;
				}
			}
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
				deriv->coeffRef(nvertices + 2 * i) = derivatives[i][3];
				deriv->coeffRef(nvertices + 2 * i + 1) = derivatives[i][4];
				for (int j = 0; j < 3; j++)
				{
					int vid = faces[i].v[j];
					int edgeId = edgeIDOppVertexIDPerface[i].v[j];
					int vidTilda = -1;

					vidTilda = vtexIDOppositeEdge(edgeId, 0) == vid ? vtexIDOppositeEdge(edgeId, 1) : vtexIDOppositeEdge(edgeId, 0);
					deriv->coeffRef(vid) += derivatives[i][j];
					// deriv->coeffRef(nvertices + ) += derivatives[i][3+j];
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

	void TFWBending::getEnergy(Sheet& sheet, Real& energy, VectorN* deriv, std::vector<Eigen::Triplet<Real>>* hess) const
	{
		if (hess)
			energy = getBendingEnergyAmpsPerVtexDphis1VectorPerFace(sheet, deriv, hess);
		else if (deriv)
			energy = getBendingEnergyAmpsPerVtexDphis1VectorPerFace(sheet, deriv, NULL);
		else
			energy = getBendingEnergyAmpsPerVtexDphis1VectorPerFace(sheet, NULL, NULL);
	}

	void TFWBending::getWrinkleShellEnergy(Sheet& sheet, Real& energy, std::vector<bool> params, VectorN* deriv, std::vector<Eigen::Triplet<Real>>* hess) const
	{
		if (hess)
			energy = getBendingEnergyAmpsPerVtexDphis1VectorPerFace(sheet, deriv, hess);
		else if (deriv)
			energy = getBendingEnergyAmpsPerVtexDphis1VectorPerFace(sheet, deriv, NULL);
		else
			energy = getBendingEnergyAmpsPerVtexDphis1VectorPerFace(sheet, NULL, NULL);
	}

	void TFWBending::getWrinkleShellEnergyPerFace(Sheet& sheet, Real& energy, int faceId, VectorN* deriv, MatrixNN* hess) const
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

				energy += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * sheet.getThickness() * sheet.getThickness() * getBendingEnergyDensityFromQuadAmpsPerVtexDphis1VectorPerFace(faceId, quad, sheet, &tempDeriv, nullptr) / 12.0;
				derivatives += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * sheet.getThickness() * sheet.getThickness() * tempDeriv / 12.0;
				// std::cout <<i << " \n" <<derivatives[i] << std::endl;
			}
			else
			{
				energy += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * sheet.getThickness() * sheet.getThickness() * getBendingEnergyDensityFromQuadAmpsPerVtexDphis1VectorPerFace(faceId, quad, sheet, nullptr, nullptr) / 12.0;
			}
		}
		// assuming that this function is used only for dphi per face optimization
		if (deriv)
		{
			deriv->resize(3 + 2 + 3 * 3 + 3 * 3);
			// (*deriv)(0) = derivatives(3);
			// (*deriv)(1) = derivatives(4);
			(*deriv).setZero();
			(*deriv) += derivatives;
		}
		if (hess)
			throw "Not Implemented!";
	}

	void TFWBending::addForces(Sheet& sheet,
	    Matrix3N& forces) const
	{
		throw "Not Implemented!";
	}
	void TFWBending::addForcesPerFace(Sheet& sheet, Eigen::Vector<Real, 23>& forces, int faceId) const
	{
		// ensure that forces is initialized with zeros
		// TODO: Delete
		bool _debug = false;
		// forces is set to zero
		const std::vector<trimesh::triangle_t>& faces = sheet.getTriangles();
		const Matrix3N& positions = sheet.getPositionsN3();
		const VectorN& amps = sheet.getAmpVec();
		const Matrix2N& dphis = sheet.getDphis1VectorPerFace();
		Real alpha = sheet.getLameAlpha();
		Real beta = sheet.getLameBeta();
		Real hcube = sheet.getThickness() * sheet.getThickness() * sheet.getThickness();
		Matrix33 facePositions;
		int vidxi = faces[faceId].v[0];
		int vidxj = faces[faceId].v[1];
		int vidxk = faces[faceId].v[2];

		facePositions.col(0) = positions.col(faces[faceId].v[0]);
		facePositions.col(1) = positions.col(faces[faceId].v[1]);
		facePositions.col(2) = positions.col(faces[faceId].v[2]);
		Eigen::Matrix<Real, 4, 18> gradII; // gradI(i) contains dI(i)/dX
		std::shared_ptr<SecondFundamentalFormDiscretization> sff;
		sff = std::make_shared<MidedgeAverageFormulation>();
		Matrix22 II = sff->secondFundamentalFormNew(sheet.getBaseMesh(), sheet.getPositionsN3(), sheet.getEdgeDofs(), faceId, &gradII, NULL);
		Matrix22 Ibar = sheet.getRestShape().getIbars()[faceId];
		Matrix22 Ibarinv = Ibar.inverse();
		Matrix22 IIbar = sheet.getRestShape().getIIbars()[faceId];

		Real bendCoeff = hcube / 12.0 * std::sqrt(Ibar.determinant()) * 0.5;
		// if (faceId == 8)
		// {
		// 	std::cout << "Area " << std::sqrt(Ibar.determinant()) << std::endl;
		// 	std::cout << "Thickness " << sheet.getThickness() << std::endl;
		// 	std::cout << bendCoeff << std::endl;
		// }
		Vector2 dphi = dphis.col(faceId);

		Real a0 = amps(vidxi);
		Real a1 = amps(vidxj);
		Real a2 = amps(vidxk);

		Matrix22 dphidphiT = dphi * dphi.transpose(); // phi0 * phi0, phi0 * phi1, phi1 * phi0, phi1 * phi1
		Matrix22 derivdphidphiTphi0;
		derivdphidphiTphi0 << 2 * dphi(0), dphi(1), dphi(1), 0;
		Matrix22 derivdphidphiTphi1;
		derivdphidphiTphi1 << 0, dphi(0), dphi(0), 2 * dphi(1);

		Matrix22 M1 = Ibarinv * (II - IIbar);
		Matrix22 dM1SVdM1 = alpha * M1.trace() * Matrix22::Identity() + 2 * beta * M1;

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

				std::cout << "Ibar \n"
				          << Ibar << std::endl;
				std::cout << "II \n"
				          << II << std::endl;
			}

			// M1
			if (_debug && faceId == 8)
			{
				std::cout << "M1\n"
				          << M1 << std::endl;
				// std::cout << "Stretch Coeff " << bendCoeff << " quadWeight " << quadWeight << std::endl;
			}

			// wrinkle params
			// pos derivatives
			for (int xii = 0; xii < 18; xii++)
			{
				Matrix22 dM1dxii = Ibarinv * (gradII.col(xii).reshaped(2, 2));
				forces(5 + xii) += bendCoeff * quadWeight * doubleContraction(dM1SVdM1, dM1dxii);
			}

			// M2
			Matrix22 M2 = a * Ibarinv * dphidphiT;
			Matrix22 dM2SVdM2 = alpha * M2.trace() * Matrix22::Identity() + 2 * beta * M2;

			// amps derivatives
			Matrix22 dM2da0 = Ibarinv * (dphidphiT * ((1 - u - v)));
			Matrix22 dM2da1 = Ibarinv * (dphidphiT * (u));
			Matrix22 dM2da2 = Ibarinv * (dphidphiT * (v));
			forces(0) += 0.5 * bendCoeff * quadWeight * doubleContraction(dM2SVdM2, dM2da0);
			forces(1) += 0.5 * bendCoeff * quadWeight * doubleContraction(dM2SVdM2, dM2da1);
			forces(2) += 0.5 * bendCoeff * quadWeight * doubleContraction(dM2SVdM2, dM2da2);
			// dphi derivatives
			Matrix22 dM2dphi0 = Ibarinv * (a * derivdphidphiTphi0);
			Matrix22 dM2dphi1 = Ibarinv * (a * derivdphidphiTphi1);
			forces(3) += 0.5 * bendCoeff * quadWeight * doubleContraction(dM2SVdM2, dM2dphi0);
			forces(4) += 0.5 * bendCoeff * quadWeight * doubleContraction(dM2SVdM2, dM2dphi1);
		}
	}

	void TFWBending::addForces(Sheet& sheet, VectorN& forces) const
	{
		// TODO: Delete
		bool _debug = false;
		// forces is set to zero
		int nfaces = sheet.getFaceCount();
		int nvertices = sheet.getVertexCount();
		const std::vector<trimesh::triangle_t>& faces = sheet.getTriangles();
		const Matrix3N& positions = sheet.getPositionsN3();
		const VectorN& amps = sheet.getAmpVec();
		const Matrix2N& dphis = sheet.getDphis1VectorPerFace();
		Real alpha = sheet.getLameAlpha();
		Real beta = sheet.getLameBeta();
		Real hcube = sheet.getThickness() * sheet.getThickness() * sheet.getThickness();
		std::vector<Eigen::Vector<Real, 23>> faceForces;
		faceForces.resize(nfaces);

#ifdef USE_OMP
#pragma omp parallel for
#endif
		for (int faceId = 0; faceId < nfaces; faceId++)
		{
			// Enable for debugging the per face functions
			// Eigen::Vector<Real, 23> faceDeriv;
			// faceDeriv.setZero();
			// addForcesPerFace(sheet, faceDeriv, faceId);
			// getEnergyPerFace(sheet, &faceDeriv, faceId);
			Eigen::Vector<Real, 23> faceDeriv;
			faceDeriv.setZero();
			Matrix33 facePositions;
			int vidxi = faces[faceId].v[0];
			int vidxj = faces[faceId].v[1];
			int vidxk = faces[faceId].v[2];

			facePositions.col(0) = positions.col(faces[faceId].v[0]);
			facePositions.col(1) = positions.col(faces[faceId].v[1]);
			facePositions.col(2) = positions.col(faces[faceId].v[2]);
			Eigen::Matrix<Real, 4, 18> gradII; // gradI(i) contains dI(i)/dX
			std::shared_ptr<SecondFundamentalFormDiscretization> sff;
			sff = std::make_shared<MidedgeAverageFormulation>();
			Matrix22 II = sff->secondFundamentalFormNew(sheet.getBaseMesh(), sheet.getPositionsN3(), sheet.getEdgeDofs(), faceId, &gradII, NULL);
			Matrix22 Ibar = sheet.getRestShape().getIbars()[faceId];
			Matrix22 Ibarinv = Ibar.inverse();
			Matrix22 IIbar = sheet.getRestShape().getIIbars()[faceId];

			Real bendCoeff = hcube / 12.0 * std::sqrt(Ibar.determinant()) * 0.5;
			// if (faceId == 8)
			// {
			// 	std::cout << "Area " << std::sqrt(Ibar.determinant()) << std::endl;
			// 	std::cout << "Thickness " << sheet.getThickness() << std::endl;
			// 	std::cout << bendCoeff << std::endl;
			// }
			Vector2 dphi = dphis.col(faceId);

			Real a0 = amps(vidxi);
			Real a1 = amps(vidxj);
			Real a2 = amps(vidxk);

			Matrix22 dphidphiT = dphi * dphi.transpose(); // phi0 * phi0, phi0 * phi1, phi1 * phi0, phi1 * phi1
			Matrix22 derivdphidphiTphi0;
			derivdphidphiTphi0 << 2 * dphi(0), dphi(1), dphi(1), 0;
			Matrix22 derivdphidphiTphi1;
			derivdphidphiTphi1 << 0, dphi(0), dphi(0), 2 * dphi(1);

			Matrix22 M1 = Ibarinv * (II - IIbar);
			Matrix22 dM1SVdM1 = alpha * M1.trace() * Matrix22::Identity() + 2 * beta * M1;

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

					std::cout << "Ibar \n"
					          << Ibar << std::endl;
					std::cout << "II \n"
					          << II << std::endl;
				}

				// M1
				if (_debug && faceId == 8)
				{
					std::cout << "M1\n"
					          << M1 << std::endl;
					// std::cout << "Stretch Coeff " << bendCoeff << " quadWeight " << quadWeight << std::endl;
				}

				// wrinkle params
				// pos derivatives
				for (int xii = 0; xii < 18; xii++)
				{
					Matrix22 dM1dxii = Ibarinv * (gradII.col(xii).reshaped(2, 2));
					faceDeriv(5 + xii) += bendCoeff * quadWeight * doubleContraction(dM1SVdM1, dM1dxii);
				}

				// M2
				Matrix22 M2 = a * Ibarinv * dphidphiT;
				Matrix22 dM2SVdM2 = alpha * M2.trace() * Matrix22::Identity() + 2 * beta * M2;

				// amps derivatives
				Matrix22 dM2da0 = Ibarinv * (dphidphiT * ((1 - u - v)));
				Matrix22 dM2da1 = Ibarinv * (dphidphiT * (u));
				Matrix22 dM2da2 = Ibarinv * (dphidphiT * (v));
				faceDeriv(0) += 0.5 * bendCoeff * quadWeight * doubleContraction(dM2SVdM2, dM2da0);
				faceDeriv(1) += 0.5 * bendCoeff * quadWeight * doubleContraction(dM2SVdM2, dM2da1);
				faceDeriv(2) += 0.5 * bendCoeff * quadWeight * doubleContraction(dM2SVdM2, dM2da2);
				// dphi derivatives
				Matrix22 dM2dphi0 = Ibarinv * (a * derivdphidphiTphi0);
				Matrix22 dM2dphi1 = Ibarinv * (a * derivdphidphiTphi1);
				faceDeriv(3) += 0.5 * bendCoeff * quadWeight * doubleContraction(dM2SVdM2, dM2dphi0);
				faceDeriv(4) += 0.5 * bendCoeff * quadWeight * doubleContraction(dM2SVdM2, dM2dphi1);
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
	Real TFWBending::getEnergyPerFace(Sheet& sheet, Eigen::Vector<Real, 23>* forces, int faceId) const
	{
		bool _debug = false;
		// forces is set to zero
		const std::vector<trimesh::triangle_t>& faces = sheet.getTriangles();
		const Matrix3N& positions = sheet.getPositionsN3();
		const VectorN& amps = sheet.getAmpVec();
		const Matrix2N& dphis = sheet.getDphis1VectorPerFace();
		Real alpha = sheet.getLameAlpha();
		Real beta = sheet.getLameBeta();
		Real hcube = sheet.getThickness() * sheet.getThickness() * sheet.getThickness();
		Matrix33 facePositions;
		int vidxi = faces[faceId].v[0];
		int vidxj = faces[faceId].v[1];
		int vidxk = faces[faceId].v[2];

		facePositions.col(0) = positions.col(faces[faceId].v[0]);
		facePositions.col(1) = positions.col(faces[faceId].v[1]);
		facePositions.col(2) = positions.col(faces[faceId].v[2]);
		Eigen::Matrix<Real, 4, 18> gradII; // gradI(i) contains dI(i)/dX
		std::shared_ptr<SecondFundamentalFormDiscretization> sff;
		sff = std::make_shared<MidedgeAverageFormulation>();
		Matrix22 II = sff->secondFundamentalFormNew(sheet.getBaseMesh(), sheet.getPositionsN3(), sheet.getEdgeDofs(), faceId, &gradII, NULL);
		Matrix22 Ibar = sheet.getRestShape().getIbars()[faceId];
		Matrix22 Ibarinv = Ibar.inverse();
		Matrix22 IIbar = sheet.getRestShape().getIIbars()[faceId];

		Real bendCoeff = hcube / 12.0 * std::sqrt(Ibar.determinant()) * 0.5;
		// if (faceId == 8)
		// {
		// 	std::cout << "Area " << std::sqrt(Ibar.determinant()) << std::endl;
		// 	std::cout << "Thickness " << sheet.getThickness() << std::endl;
		// 	std::cout << bendCoeff << std::endl;
		// }
		Vector2 dphi = dphis.col(faceId);

		Real a0 = amps(vidxi);
		Real a1 = amps(vidxj);
		Real a2 = amps(vidxk);

		Matrix22 dphidphiT = dphi * dphi.transpose(); // phi0 * phi0, phi0 * phi1, phi1 * phi0, phi1 * phi1
		Matrix22 derivdphidphiTphi0;
		derivdphidphiTphi0 << 2 * dphi(0), dphi(1), dphi(1), 0;
		Matrix22 derivdphidphiTphi1;
		derivdphidphiTphi1 << 0, dphi(0), dphi(0), 2 * dphi(1);

		Matrix22 M1 = Ibarinv * (II - IIbar);

		Matrix22 dM1SVdM1 = alpha * M1.trace() * Matrix22::Identity() + 2 * beta * M1;
		Real energy = 0;
		energy += bendCoeff * SVNormSquared(alpha, beta, M1);
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

				std::cout << "Ibar \n"
				          << Ibar << std::endl;
				std::cout << "II \n"
				          << II << std::endl;
			}

			// M1
			if (_debug && faceId == 8)
			{
				std::cout << "M1\n"
				          << M1 << std::endl;
				// std::cout << "Stretch Coeff " << bendCoeff << " quadWeight " << quadWeight << std::endl;
			}

			// wrinkle params
			// pos derivatives
			if (forces)
			{
				// for (int xii = 0; xii < 18; xii++)
				// {
				// 	// TODO Optimize here when required
				// 	Matrix22 dM1dxii = Ibarinv * (gradII.col(xii).reshaped(2, 2));
				// 	forces->coeffRef(5 + xii) += bendCoeff * quadWeight * doubleContraction(dM1SVdM1, dM1dxii);
				// }
			}

			// M2
			Matrix22 M2 = a * Ibarinv * dphidphiT;
			energy += 0.5 * SVNormSquared(alpha, beta, M2);
			if (forces)
			{
				Matrix22 dM2SVdM2 = alpha * M2.trace() * Matrix22::Identity() + 2 * beta * M2;

				// amps derivatives
				// Matrix22 dM2da0 = Ibarinv * (dphidphiT * ((1 - u - v)));
				// Matrix22 dM2da1 = Ibarinv * (dphidphiT * (u));
				// Matrix22 dM2da2 = Ibarinv * (dphidphiT * (v));
				// forces->coeffRef(0) += 0.5 * bendCoeff * quadWeight * doubleContraction(dM2SVdM2, dM2da0);
				// forces->coeffRef(1) += 0.5 * bendCoeff * quadWeight * doubleContraction(dM2SVdM2, dM2da1);
				// forces->coeffRef(2) += 0.5 * bendCoeff * quadWeight * doubleContraction(dM2SVdM2, dM2da2);
				// dphi derivatives
				Matrix22 dM2dphi0 = Ibarinv * (a * derivdphidphiTphi0);
				Matrix22 dM2dphi1 = Ibarinv * (a * derivdphidphiTphi1);
				forces->coeffRef(3) += 0.5 * bendCoeff * quadWeight * doubleContraction(dM2SVdM2, dM2dphi0);
				forces->coeffRef(4) += 0.5 * bendCoeff * quadWeight * doubleContraction(dM2SVdM2, dM2dphi1);
			}
		}
		return energy;
	}

};
