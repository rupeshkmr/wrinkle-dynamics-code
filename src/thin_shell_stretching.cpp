#include <argus/thin_shell_stretching.hpp>
#include <argus/util.hpp>
#include <argus/geometry_functions.hpp>
#include "../../vendor/SecondFundamentalForm/SecondFundamentalFormDiscretization.h"
#include "../../vendor/SecondFundamentalForm/MidedgeAverageFormulation.h"
/* #include "../vendor/MeshLib/IntrinsicGeometry.h"
#include "../../vendor/MeshLib/MeshConnectivity.h"
#include "../../vendor/MeshLib/MeshGeometry.h" */
#define idxAt(i, j) 3 * i + j
#define mat(A, i, j) A(0, idxAt(i, j)), A(1, idxAt(i, j)), A(2, idxAt(i, j)), A(3, idxAt(i, j))
#define hessMat(A, i, j) A.at(0)(i, j), A.at(1)(i, j), A.at(2)(i, j), A.at(3)(i, j)
namespace argus
{
#ifdef ARGUS_CHECKPOINT
	void ThinShellStretching::saveInfo()
	{
		nlohmann::json data;
		data["Energy"] = "StVK Thin Shell Stretching";
		data["Is external"] = isExternal();
		simConf["Energies"].push_back(data);
	}
#endif
	// function for per vertex derivatives and hessians amps
	Real ThinShellStretching::getEnergyDensityFromQuadPerVtex1(int vid, int faceId, int quadId, Sheet& sheet, Real* deriv, Real* hess) const
	{
		throw "Not a wrinkle model";
	}

	Real ThinShellStretching::getVertexEnergy1(Sheet& sheet, int vid, Real* deriv, Real* hess) const
	{
		throw "Not a wrinkle model";
	}

	// function for per vertex derivatives and hessians amps + pos
	Real ThinShellStretching::getEnergyDensityFromQuadPerVtex4(int vid, int faceId, int quadId, Sheet& sheet, Vector4* deriv, Matrix44* hess) const
	{
		throw "Not a wrinkle model";
	}

	Real ThinShellStretching::getVertexEnergy4(Sheet& sheet, int vid, Vector4* deriv, Matrix44* hess) const
	{
		throw "Not a wrinkle model";
	}

	Real ThinShellStretching::getEnergyDensityFromQuadPerVtex(int vid, int faceId, int quadId, Sheet& sheet, Vector3* deriv, Matrix33* hess) const
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
		facePositions.col(0) = positions.block<3, 1>(0, faces[faceId].v[0]);
		facePositions.col(1) = positions.block<3, 1>(0, faces[faceId].v[1]);
		facePositions.col(2) = positions.block<3, 1>(0, faces[faceId].v[2]);
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
		Real lameAlpha = sheet.getLameAlpha();
		Real lameBeta = sheet.getLameBeta();

		if (deriv)
			deriv->setZero();
		if (hess)
			hess->setZero();

		// First term
		// Correct numerically
		Matrix22 M0 = Ibarinv * (I - Ibar);
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

		Real area = 0.5 * sqrt(Ibar.determinant());
		if (deriv)
			*deriv = (area * (*deriv)).eval();
		if (hess)
			(*hess) = (area * (*hess)).eval();
		// project
		if (hess)
			*hess = argus::util::hessianProjection(*hess, true);
		return area * energy;
	}
	Real ThinShellStretching::getVertexEnergy(Sheet& sheet, int vid, Vector3* deriv, Matrix33* hess) const
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

			Real area = 0.5 * sqrt(Ibar.determinant());

			Vector3 tempDeriv;
			tempDeriv.setZero();
			Matrix33 tempHess;
			tempHess.setZero();
			Real energy = 0;

			// First term
			// Correct numerically
			Matrix22 M0 = Ibarinv * (I - Ibar);
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
			for (int quadId = 0; quadId < sheet.getQuadPoints().size(); quadId++)
			{

				if (deriv)
					(*deriv) += sheet.getQuadPoints()[quadId].weight * area * sheet.getThickness() / 4.0 * tempDeriv;
				if (hess)
					(*hess) += sheet.getQuadPoints()[quadId].weight * area * sheet.getThickness() / 4.0 * tempHess;
				totalEnergy += sheet.getQuadPoints()[quadId].weight * area * sheet.getThickness() / 4.0 * energy;
			}
		}

		return totalEnergy;
	}
	/* Real ThinShellStretching::getVertexEnergy(Sheet& sheet, int vid, Vector3* deriv, Matrix33* hess) const
	{
	    // need to compute deriv and hessian wrt the vertex at index vid
	    // just for testing
	    Real energy = 0;
	    const Eigen::VectorXi& faces = sheet.vbdInfo.vertexFacesSet.row(vid);
	    unsigned int nfaces = 0;
	    for (int i = 0; i < faces.rows(); i++)
	    {
	        if (faces(i) == -1)
	            break;
	        nfaces++;
	    }
	    // create buffers to store parallel results
	    std::vector<Real> energies(nfaces, 0);
	    std::vector<Vector3> derivatives(nfaces);
	    if (deriv)
	        deriv->setZero();
	    if (hess)
	        hess->setZero();
	    if (deriv)
	        for (int i = 0; i < nfaces; i++)
	        {
	            derivatives[i].setZero();
	        }
	    std::vector<Matrix33> hessians(nfaces);
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
	            std::cout << "Vertex " << localVid << " not found! Face " << faceId << " Vertex ID: " << vid << ". Error in tfw_stretching.cpp\n";
	            exit(0);
	        }
	        for (int quad = 0; quad < sheet.getQuadPoints().size(); quad++)
	        {
	            Vector3 tempDeriv1, tempDeriv2, tempDeriv;
	            Matrix33 tempHess, tempHess1, tempHess2;
	            energies[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * ThinShellStretching::getEnergyDensityFromQuadPerVtex(localVid, faceId, quad, sheet, deriv ? &tempDeriv : NULL, hess ? &tempHess : NULL) / 4.0;
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
 */
	// Function for entire sheet
	Real ThinShellStretching::getStretchingEnergyDensityFromQuadAmpsPerVtexDphis1VectorPerFace(int faceId, int quadId, Sheet& sheet, VectorN* deriv, MatrixNN* hess) const
	{
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
		facePositions.col(0) = positions.block<3, 1>(0, faces[faceId].v[0]);
		facePositions.col(1) = positions.block<3, 1>(0, faces[faceId].v[1]);
		facePositions.col(2) = positions.block<3, 1>(0, faces[faceId].v[2]);
		MatrixNN gradI, gradII; // gradI(i) contains dI(i)/dX
		std::vector<Eigen::Matrix<Real, 9, 9>> hessI;
		Matrix22 I;
		I = firstFundamentalForm(facePositions, hess || deriv ? true : false, hess ? true : false, hess || deriv ? &gradI : NULL, hess ? &hessI : NULL); // correct till here
		Matrix22 Ibar = sheet.getRestShape().getIbars()[faceId];

		Real area = 0.5 * sqrt(Ibar.determinant());
		if (deriv)
			(*deriv) = area * (*deriv);
		if (hess)
		{
			(*hess) = area * (*hess);
			// project
			(*hess) = argus::util::lowRankApprox(*hess);
		}
		return area * 0;
	}
	Real ThinShellStretching::getStretchingEnergyAmpsPerVtexDphis1VectorPerFace(Sheet& sheet, VectorN* deriv, std::vector<Eigen::Triplet<Real>>* hess) const
	{
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
#ifdef USE_OMP
#pragma omp parallel for
#endif
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
					energies[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * ThinShellStretching::getStretchingEnergyDensityFromQuadAmpsPerVtexDphis1VectorPerFace(i, quad, sheet, &tempDeriv, &tempHess) / 4.0;
					derivatives[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * tempDeriv / 4.0;
					hessians[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * tempHess / 4.0;
					// std::cout <<i << " \n" <<derivatives[i] << std::endl;
				}
				else if (deriv)
				{
					VectorN tempDeriv;

					energies[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * ThinShellStretching::getStretchingEnergyDensityFromQuadAmpsPerVtexDphis1VectorPerFace(i, quad, sheet, &tempDeriv, nullptr) / 4.0;
					derivatives[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * tempDeriv / 4.0;
					// std::cout <<i << " \n" <<derivatives[i] << std::endl;
				}
				else
				{
					energies[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * ThinShellStretching::getStretchingEnergyDensityFromQuadAmpsPerVtexDphis1VectorPerFace(i, quad, sheet, nullptr, nullptr) / 4.0;
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

	void ThinShellStretching::getEnergy(Sheet& sheet, Real& energy, VectorN* deriv, std::vector<Eigen::Triplet<Real>>* hess) const
	{
		if (hess)
			energy = ThinShellStretching::getStretchingEnergyAmpsPerVtexDphis1VectorPerFace(sheet, deriv, hess);
		else if (deriv)
			energy = ThinShellStretching::getStretchingEnergyAmpsPerVtexDphis1VectorPerFace(sheet, deriv, NULL);
		else
			energy = ThinShellStretching::getStretchingEnergyAmpsPerVtexDphis1VectorPerFace(sheet, NULL, NULL);
	}

	void ThinShellStretching::getWrinkleShellEnergy(Sheet& sheet, Real& energy, VectorN* deriv, std::vector<Eigen::Triplet<Real>>* hess) const
	{
		if (hess)
			energy = ThinShellStretching::getStretchingEnergyAmpsPerVtexDphis1VectorPerFace(sheet, deriv, hess);
		else if (deriv)
			energy = ThinShellStretching::getStretchingEnergyAmpsPerVtexDphis1VectorPerFace(sheet, deriv, NULL);
		else
			energy = ThinShellStretching::getStretchingEnergyAmpsPerVtexDphis1VectorPerFace(sheet, NULL, NULL);
	}

	void ThinShellStretching::getWrinkleShellEnergyPerFace(Sheet& sheet, Real& energy, int faceId, VectorN* deriv, MatrixNN* hess) const
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

				energy += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * ThinShellStretching::getStretchingEnergyDensityFromQuadAmpsPerVtexDphis1VectorPerFace(faceId, quad, sheet, &tempDeriv, nullptr) / 4.0;
				derivatives += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * tempDeriv / 4.0;
				// std::cout <<i << " \n" <<derivatives[i] << std::endl;
			}
			else
			{
				energy += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * ThinShellStretching::getStretchingEnergyDensityFromQuadAmpsPerVtexDphis1VectorPerFace(faceId, quad, sheet, nullptr, nullptr) / 4.0;
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

	Real ThinShellStretching::getEnergy(Sheet& sheet) const
	{
		Real totalEnergy = 0;
		Real lameAlpha = sheet.getLameAlpha();
		Real lameBeta = sheet.getLameBeta();
		Real thickness = sheet.getThickness();
		unsigned int nfaces = sheet.getFaceCount();
		const std::vector<trimesh::triangle_t>& faces = sheet.getTriangles();
		const Matrix3N& positions = sheet.getPositionsN3();
		const std::vector<Matrix22>& Ibars = sheet.getRestShape().getIbars();
		Real coeff = thickness / 4.0;

		// iterate through all faces
		for (int faceId = 0; faceId < nfaces; faceId++)
		{
			// each face has three vertices, so size 3x3 matrix for storing derivatives
			// deriv.col(i) = dU/dxi
			Matrix33 deriv;
			deriv.setZero();
			Matrix33 facePositions;
			facePositions.col(0) = positions.col(faces[faceId].v[0]);
			facePositions.col(1) = positions.col(faces[faceId].v[1]);
			facePositions.col(2) = positions.col(faces[faceId].v[2]);
			Eigen::Matrix<Real, 4, 9> Ideriv;
			Matrix22 I = firstFundamentalFormNew(facePositions, &Ideriv);
			Matrix22 Ibarinv = Ibars[faceId].inverse();
			Matrix22 M = Ibarinv * (I - Ibars[faceId]);
			totalEnergy += 0.5 * std::sqrt(Ibars[faceId].determinant()) * coeff * SVNormSquared(lameAlpha, lameBeta, M);
		}
		return totalEnergy;
	}

	void ThinShellStretching::addForces(Sheet& sheet,
	    Matrix3N& forces) const
	{
		Real totalEnergy = 0;
		Real lameAlpha = sheet.getLameAlpha();
		Real lameBeta = sheet.getLameBeta();
		Real thickness = sheet.getThickness();
		unsigned int nfaces = sheet.getFaceCount();
		const std::vector<trimesh::triangle_t>& faces = sheet.getTriangles();
		const Matrix3N& positions = sheet.getPositionsN3();
		const std::vector<Matrix22>& Ibars = sheet.getRestShape().getIbars();
		std::vector<Eigen::Vector<Real, 9>> facesDerivatives;
		facesDerivatives.resize(nfaces);
// iterate through all faces
#ifdef USE_OMP
#pragma omp parallel for
#endif
		for (int faceId = 0; faceId < nfaces; faceId++)
		{
			// each face has three vertices, so size 3x3 matrix for storing derivatives
			// deriv.col(i) = dU/dxi
			// Matrix33 deriv;
			// deriv.setZero();
			Eigen::Vector<Real, 9> faceDeriv;
			faceDeriv.setZero();
			Matrix33 facePositions;
			facePositions.col(0) = positions.col(faces[faceId].v[0]);
			facePositions.col(1) = positions.col(faces[faceId].v[1]);
			facePositions.col(2) = positions.col(faces[faceId].v[2]);
			Eigen::Matrix<Real, 4, 9> Ideriv;
			Matrix22 I = firstFundamentalFormNew(facePositions, &Ideriv);
			// Matrix22 I = firstFundamentalFormNew(facePositions, true, false, &Ideriv, NULL); // correct till here
			Real coeff = 0.5 * std::sqrt(Ibars[faceId].determinant()) * thickness / 4.0;

			Matrix22 Ibarinv = Ibars[faceId].inverse();
			Matrix22 M = Ibarinv * (I - Ibars[faceId]);
			Matrix22 dMSVdM = lameAlpha * M.trace() * Matrix22::Identity() + lameBeta * 2 * M;
			for (int xii = 0; xii < 9; xii++)
			{
				Matrix22 dM1dxii = Ibarinv * (Ideriv.col(xii).reshaped(2, 2));
				faceDeriv(xii) += coeff * doubleContraction(dMSVdM, dM1dxii);
			}
			facesDerivatives[faceId] = faceDeriv;
			// Matrix22 dM1SVdM1 = alpha * M1.trace() * Matrix22::Identity() + 2 * beta * M1;

			// iterate through three vertices and fill in the forces
			// for (int localvid = 0; localvid < 3; localvid++)
			// {
			// 	// for localvid
			// 	std::array<Matrix22, 3> dIdx;
			// 	// dIdx0
			// 	dIdx[0] = Ideriv.block<4, 1>(0, 3 * localvid + 0).reshaped(2, 2);
			// 	// dIdx1
			// 	dIdx[1] = Ideriv.block<4, 1>(0, 3 * localvid + 1).reshaped(2, 2);
			// 	// dIdx2
			// 	dIdx[2] = Ideriv.block<4, 1>(0, 3 * localvid + 2).reshaped(2, 2);

			// 	// get energy derivative
			// 	// d||M||_SV^2 / dxij = d||M||/dM : Ibarinv^T dIdxij
			// 	forces.col(faces[faceId].v[localvid])(0) -= coeff * doubleContraction(dMSVdM, Ibarinv.transpose() * dIdx[0]);
			// 	forces.col(faces[faceId].v[localvid])(1) -= coeff * doubleContraction(dMSVdM, Ibarinv.transpose() * dIdx[1]);
			// 	forces.col(faces[faceId].v[localvid])(2) -= coeff * doubleContraction(dMSVdM, Ibarinv.transpose() * dIdx[2]);
			// }
		}
		// aggregate
		for (int i = 0; i < nfaces; i++)
		{
			forces.col(faces[i].v[0]) -= facesDerivatives[i].segment(0, 3);
			forces.col(faces[i].v[1]) -= facesDerivatives[i].segment(3, 3);
			forces.col(faces[i].v[2]) -= facesDerivatives[i].segment(6, 3);
		}
	}

};
