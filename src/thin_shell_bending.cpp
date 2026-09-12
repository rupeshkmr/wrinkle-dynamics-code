#include <argus/thin_shell_bending.hpp>
#include <argus/util.hpp>
#define idxAt(i, j) 3 * i + j
#define mat(A, i, j) A(0, idxAt(i, j)), A(1, idxAt(i, j)), A(2, idxAt(i, j)), A(3, idxAt(i, j))
#define hessMat(A, i, j) A.at(0)(i, j), A.at(1)(i, j), A.at(2)(i, j), A.at(3)(i, j)

namespace argus
{
#ifdef ARGUS_CHECKPOINT
	void ThinShellBending::saveInfo()
	{
		nlohmann::json data;
		data["Energy"] = "TFW Bending";
		data["Is external"] = isExternal();
		simConf["Energies"].push_back(data);
	}
#endif
	Real ThinShellBending::getEnergyDensityFromQuadPerVtex1(int vid, int faceId, int quadId, Sheet& sheet, Real* deriv, Real* hess) const
	{
		throw "Not wrinkle energy!";
	}

	Real ThinShellBending::getVertexEnergy1(Sheet& sheet, int vid, Real* deriv, Real* hess) const
	{
		throw "Not wrinkle energy!";
	}

	Real ThinShellBending::getEnergyDensityFromQuadPerVtex4(int vid, int faceId, int quadId, Sheet& sheet, Vector4* deriv, Matrix44* hess) const
	{
		throw "Not wrinkle energy!";
	}

	Real ThinShellBending::getVertexEnergy4(Sheet& sheet, int vid, Vector4* deriv, Matrix44* hess) const

	{
		throw "Not wrinkle energy!";
	}

	Real ThinShellBending::getEnergyDensityFromQuadPerVtex(int vid, int faceId, int quadId, Sheet& sheet, Vector3* deriv, Matrix33* hess) const
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

	/* 	Real ThinShellBending::getVertexEnergy(Sheet& sheet, int vid, Vector3* deriv, Matrix33* hess) const
	    {
	        // need to compute deriv and hessian wrt the vertex at index vid
	        // just for testing
	        Real energy = 0;
	        Eigen::VectorXi faces = sheet.vbdInfo.vertexFacesSet.row(vid);
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
	        Real coeff = std::pow(sheet.getThickness(), 3) / 12.0;
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
	                std::cout << "Vertex not found! Error in thin_shell_bending.cpp\n";
	                exit(0);
	            }
	            for (int quad = 0; quad < sheet.getQuadPoints().size(); quad++)
	            {
	                Vector3 tempDeriv1, tempDeriv2, tempDeriv;
	                Matrix33 tempHess, tempHess1, tempHess2;
	                energies[i] += sheet.getQuadPoints()[quad].weight * coeff * ThinShellBending::getEnergyDensityFromQuadPerVtex(localVid, faceId, quad, sheet, deriv ? &tempDeriv : NULL, hess ? &tempHess : NULL);
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
	 */
	Real ThinShellBending::getVertexEnergy(Sheet& sheet, int vid, Vector3* deriv, Matrix33* hess) const
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

			for (int quadId = 0; quadId < sheet.getQuadPoints().size(); quadId++)
			{
				if (deriv)
					(*deriv) += coeff * d0 * sheet.getQuadPoints()[quadId].weight;
				if (hess)
					(*hess) += coeff * h0 * sheet.getQuadPoints()[quadId].weight;
				totalEnergy += coeff * e0 * sheet.getQuadPoints()[quadId].weight;
			}
		}

		return totalEnergy;
	}

	Real ThinShellBending::getBendingEnergyDensityFromQuadAmpsPerVtexDphis1VectorPerFace(int faceId, int quadId, Sheet& sheet, VectorN* deriv, MatrixNN* hess) const
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

	Real ThinShellBending::getBendingEnergyAmpsPerVtexDphis1VectorPerFace(Sheet& sheet, VectorN* deriv, std::vector<Eigen::Triplet<Real>>* hess) const
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
					energies[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * sheet.getThickness() * sheet.getThickness() * ThinShellBending::getBendingEnergyDensityFromQuadAmpsPerVtexDphis1VectorPerFace(i, quad, sheet, &tempDeriv, &tempHess) / 12.0;
					derivatives[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * sheet.getThickness() * sheet.getThickness() * tempDeriv / 12.0;
					hessians[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * sheet.getThickness() * sheet.getThickness() * tempHess / 12.0;
				}
				else if (deriv)
				{
					VectorN tempDeriv;
					energies[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * sheet.getThickness() * sheet.getThickness() * ThinShellBending::getBendingEnergyDensityFromQuadAmpsPerVtexDphis1VectorPerFace(i, quad, sheet, &tempDeriv, nullptr) / 12.0;
					derivatives[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * sheet.getThickness() * sheet.getThickness() * tempDeriv / 12.0;
				}
				else
				{
					energies[i] += sheet.getQuadPoints()[quad].weight * sheet.getThickness() * sheet.getThickness() * sheet.getThickness() * ThinShellBending::getBendingEnergyDensityFromQuadAmpsPerVtexDphis1VectorPerFace(i, quad, sheet, nullptr, nullptr) / 12.0;
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
	void ThinShellBending::getEnergy(Sheet& sheet, Real& energy, VectorN* deriv, std::vector<Eigen::Triplet<Real>>* hess) const
	{
		if (hess)
			energy = ThinShellBending::getBendingEnergyAmpsPerVtexDphis1VectorPerFace(sheet, deriv, hess);
		else if (deriv)
			energy = ThinShellBending::getBendingEnergyAmpsPerVtexDphis1VectorPerFace(sheet, deriv, NULL);
		else
			energy = ThinShellBending::getBendingEnergyAmpsPerVtexDphis1VectorPerFace(sheet, NULL, NULL);
	}

	void ThinShellBending::getWrinkleShellEnergy(Sheet& sheet, Real& energy, std::vector<bool> params, VectorN* deriv, std::vector<Eigen::Triplet<Real>>* hess) const
	{
		throw "Function getWrinkleShellEnergy not defined for ThinShellBending!";
	}

	void ThinShellBending::getWrinkleShellEnergyPerFace(Sheet& sheet, Real& energy, int faceId, VectorN* deriv, MatrixNN* hess) const
	{

		throw "Function getWrinkleShellEnergyPerFace not defined for ThinShellBending!";
	}

	Real ThinShellBending::getEnergy(Sheet& sheet) const
	{
		Real totalEnergy = 0;
		Real lameAlpha = sheet.getLameAlpha();
		Real lameBeta = sheet.getLameBeta();
		unsigned int nfaces = sheet.getFaceCount();
		const std::vector<trimesh::triangle_t>& faces = sheet.getTriangles();
		const Matrix3N& positions = sheet.getPositionsN3();
		const std::vector<Matrix22>& Ibars = sheet.getRestShape().getIbars();
		const std::vector<Matrix22>& IIbars = sheet.getRestShape().getIIbars();
		std::shared_ptr<SecondFundamentalFormDiscretization> sff;
		sff = std::make_shared<MidedgeAverageFormulation>();
		Real thickness = sheet.getThickness();

		Real coeff = thickness * thickness * thickness / 12.0;

		// iterate through all faces
		for (int faceId = 0; faceId < nfaces; faceId++)
		{
			// each face has three vertices, so size 3x3 matrix for storing derivatives
			// deriv.col(i) = dU/dxi
			Matrix22 Ibarinv = Ibars[faceId].inverse();
			Matrix22 II = sff->secondFundamentalForm(sheet.getBaseMesh(), sheet.getPositionsN3().transpose(), sheet.getEdgeDofs(), faceId, NULL, NULL);
			// Matrix22 II = sff->secondFundamentalFormNew(sheet.getBaseMesh(), sheet.getPositionsN3(), sheet.getEdgeDofs(), faceId, NULL, NULL);
			Matrix22 M = Ibarinv * (II - IIbars[faceId]);
			totalEnergy += coeff * SVNormSquared(lameAlpha, lameBeta, M) * std::sqrt(Ibars[faceId].determinant()) * 0.5;
		}
		return totalEnergy;
	}

	void ThinShellBending::addForces(Sheet& sheet,
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
		const std::vector<Matrix22>& IIbars = sheet.getRestShape().getIIbars();
		std::shared_ptr<SecondFundamentalFormDiscretization> sff;
		sff = std::make_shared<MidedgeAverageFormulation>();
		const std::vector<trimesh::triangle_t>& edgeIDOppVertexIDPerface = sheet.getVertexOppositeEdgesPerFace();
		const Eigen::MatrixXi& vtexIDOppositeEdge = sheet.getEdgeOppositeVertices();
		Real coeff = thickness * thickness * thickness / 12.0;
		std::vector<Eigen::Vector<Real, 18>> facesDerivatives;
		facesDerivatives.resize(nfaces);
// iterate through all faces
#ifdef USE_OMP
#pragma omp parallel for
#endif
		for (int faceId = 0; faceId < nfaces; faceId++)
		{
			// each face has three vertices, so size 3x3 matrix for storing derivatives
			// deriv.col(i) = dU/dxi
			Eigen::Vector<Real, 18> faceDeriv;
			faceDeriv.setZero();
			Real localCoeff = 0.5 * std::sqrt(Ibars[faceId].determinant()) * coeff;
			// Eigen::Matrix<Real, 4, 18> IIderiv;
			// Matrix22 II = sff->secondFundamentalFormNew(sheet.getBaseMesh(), sheet.getPositionsN3(), sheet.getEdgeDofs(), faceId, &IIderiv, NULL);
			Eigen::Matrix<Real, 4, 18> IIderiv;
			Matrix22 II = sff->secondFundamentalFormNew(sheet.getBaseMesh(), sheet.getPositionsN3(), sheet.getEdgeDofs(), faceId, &IIderiv, NULL);

			Matrix22 Ibarinv = Ibars[faceId].inverse();
			Matrix22 M = Ibarinv * (II - IIbars[faceId]);
			Matrix22 dMSVdM = lameAlpha * M.trace() * Matrix22::Identity() + lameBeta * 2 * M;
			// TODO: check for correctedness
			// iterate through three vertices and fill in the forces
			for (int localvid = 0; localvid < 3; localvid++)
			{
				// for localvid
				std::array<Matrix22, 6> dIIdx;
				// dIdx0
				dIIdx[0] = IIderiv.block<4, 1>(0, 3 * localvid + 0).reshaped(2, 2);
				// dIdx1
				dIIdx[1] = IIderiv.block<4, 1>(0, 3 * localvid + 1).reshaped(2, 2);
				// dIdx2
				dIIdx[2] = IIderiv.block<4, 1>(0, 3 * localvid + 2).reshaped(2, 2);

				// dIdx0
				dIIdx[3] = IIderiv.block<4, 1>(0, 9 + 3 * localvid + 0).reshaped(2, 2);
				// dIdx1
				dIIdx[4] = IIderiv.block<4, 1>(0, 9 + 3 * localvid + 1).reshaped(2, 2);
				// dIdx2
				dIIdx[5] = IIderiv.block<4, 1>(0, 9 + 3 * localvid + 2).reshaped(2, 2);

				// get energy derivative
				// d||M||_SV^2 / dxij = d||M||/dM : Ibarinv^T dIdxij
				faceDeriv(3 * localvid) -= localCoeff * doubleContraction(dMSVdM, Ibarinv.transpose() * dIIdx[0]);
				faceDeriv(3 * localvid + 1) -= localCoeff * doubleContraction(dMSVdM, Ibarinv.transpose() * dIIdx[1]);
				faceDeriv(3 * localvid + 2) -= localCoeff * doubleContraction(dMSVdM, Ibarinv.transpose() * dIIdx[2]);
				faceDeriv(3 * localvid + 9) -= localCoeff * doubleContraction(dMSVdM, Ibarinv.transpose() * dIIdx[3]);
				faceDeriv(3 * localvid + 10) -= localCoeff * doubleContraction(dMSVdM, Ibarinv.transpose() * dIIdx[4]);
				faceDeriv(3 * localvid + 11) -= localCoeff * doubleContraction(dMSVdM, Ibarinv.transpose() * dIIdx[5]);
			}
			facesDerivatives[faceId] = faceDeriv;
		}
		for (int faceId = 0; faceId < nfaces; faceId++)
		{
			forces.col(faces[faceId].v[0]) += facesDerivatives[faceId].segment(0, 3);
			forces.col(faces[faceId].v[1]) += facesDerivatives[faceId].segment(3, 3);
			forces.col(faces[faceId].v[2]) += facesDerivatives[faceId].segment(6, 3);
			// fill forces for opposite vertex to opposite edge
			for (int localvid = 0; localvid < 3; localvid++)
			{
				int edgeId = edgeIDOppVertexIDPerface[faceId].v[localvid];
				int vidTilda = -1;
				vidTilda = vtexIDOppositeEdge(edgeId, 0) == faces[faceId].v[localvid] ? vtexIDOppositeEdge(edgeId, 1) : vtexIDOppositeEdge(edgeId, 0);
				if (vidTilda != -1)
				{
					forces.col(vidTilda) += facesDerivatives[faceId].segment(9 + 3 * localvid, 3);
				}
			}
		}
	}

};
