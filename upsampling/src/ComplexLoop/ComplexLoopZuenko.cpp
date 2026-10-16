#include "../../include/CommonTools.h"
#include "../../include/ComplexLoop/ComplexLoopZuenko.h"
#include <iostream>
#include <cassert>
#include <memory>
#include <chrono>
#include <cstdlib>
#include <string>


std::complex<double> ComplexLoopZuenko::interpZ(const std::vector<std::complex<double>>& zList, const std::vector<Eigen::Vector3d>& gradThetaList, std::vector<double>& coords, const std::vector<Eigen::Vector3d>& pList)
{
	Eigen::Vector3d P = Eigen::Vector3d::Zero();
	
	int n = pList.size();
	for (int i = 0; i < n; i++)
	{
		P += coords[i] * pList[i];
	}
	Eigen::Vector3d gradThetaAve = Eigen::Vector3d::Zero();
	for (int i = 0; i < n; i++)
	{
		gradThetaAve += coords[i] * gradThetaList[i];
	}

	std::complex<double> z = 0;

	std::complex<double> I = std::complex<double>(0, 1);

	std::vector<double> deltaThetaList(n);
	for (int i = 0; i < n; i++)
	{
		deltaThetaList[i] = (P - pList[i]).dot(gradThetaAve);
	}

	for (int i = 0; i < n; i++)
	{
		z += coords[i] * zList[i] * std::complex<double>(std::cos(deltaThetaList[i]), std::sin(deltaThetaList[i]));
	}

	return z;
}

void ComplexLoopZuenko::interpZ(
	const std::vector<std::complex<double>>& zList1, const std::vector<Eigen::Vector3d>& gradThetaList1,
	const std::vector<std::complex<double>>& zList2, const std::vector<Eigen::Vector3d>& gradThetaList2,
	const std::vector<double>& coords, const std::vector<Eigen::Vector3d>& pList,
	std::complex<double>& z1, std::complex<double>& z2)
{
	Eigen::Vector3d P = Eigen::Vector3d::Zero();
	
	int n = pList.size();
	for (int i = 0; i < n; i++)
	{
		P += coords[i] * pList[i];
	}

	std::vector<double> deltaThetaList1(n);
	std::vector<double> deltaThetaList2(n);
	
	Eigen::Vector3d gradThetaAve1 = Eigen::Vector3d::Zero();
	Eigen::Vector3d gradThetaAve2 = Eigen::Vector3d::Zero();
	for (int i = 0; i < n; i++)
	{
		gradThetaAve1 += coords[i] * gradThetaList1[i];
		gradThetaAve2 += coords[i] * gradThetaList2[i];
	}

	for (int i = 0; i < n; i++)
	{
		deltaThetaList1[i] = (P - pList[i]).dot(gradThetaAve1);
		deltaThetaList2[i] = (P - pList[i]).dot(gradThetaAve2);
	}

	z1 = 0;
	z2 = 0;
	for (int i = 0; i < n; i++)
	{
		z1 += coords[i] * zList1[i] * std::complex<double>(std::cos(deltaThetaList1[i]), std::sin(deltaThetaList1[i]));
		z2 += coords[i] * zList2[i] * std::complex<double>(std::cos(deltaThetaList2[i]), std::sin(deltaThetaList2[i]));
	}
}

std::complex<double> ComplexLoopZuenko::computeZandGradZ(const Eigen::VectorXd& omega, const std::vector<std::complex<double>>& zvals, int fid, const Eigen::Vector3d& bary, Eigen::Vector3cd* gradz, Mesh &mesh)
{
	Eigen::Matrix3d gradBary;
	Eigen::Vector3d P0 = _mesh.GetVertPos(mesh.GetFaceVerts(fid)[0]);
	Eigen::Vector3d P1 = _mesh.GetVertPos(mesh.GetFaceVerts(fid)[1]);
	Eigen::Vector3d P2 = _mesh.GetVertPos(mesh.GetFaceVerts(fid)[2]);
	
	computeBaryGradient(P0, P1, P2, bary, gradBary);

	std::complex<double> z = 0;

	if (gradz)
	{
		gradz->setZero();
	}
	
	
	std::complex<double> tau = std::complex<double>(0, 1);

	for (int i = 0; i < 3; i++)
	{
		double alphai = bary(i);
		int vid = mesh.GetFaceVerts(fid)[i];
		int eid0 = mesh.GetFaceEdges(fid)[i];
		int eid1 = mesh.GetFaceEdges(fid)[(i + 2) % 3];

		double w1 = omega(eid0);
		double w2 = omega(eid1);

		if (vid == mesh.GetEdgeVerts(eid0)[1])
			w1 *= -1;
		if (vid == mesh.GetEdgeVerts(eid1)[1])
			w2 *= -1;

		double deltaTheta = w1 * bary((i + 1) % 3) + w2 * bary((i + 2) % 3);
		std::complex<double> expip = std::complex<double>(std::cos(deltaTheta), std::sin(deltaTheta));

		z += alphai * zvals[vid] * expip;

		if(gradz)
			(*gradz) += zvals[vid] * expip * (gradBary.row(i) + tau * alphai * w1 * gradBary.row((i + 1) % 3) + tau * alphai * w2 * gradBary.row((i + 2) % 3));
	}
	return z;
}

Eigen::Vector3d ComplexLoopZuenko::computeGradThetaFromOmegaPerface(const Eigen::VectorXd& omega, int fid, int vInF, Mesh &mesh)
{
    int vid = mesh.GetFaceVerts(fid)[vInF];
    int eid0 = mesh.GetFaceEdges(fid)[vInF];
    int eid1 = mesh.GetFaceEdges(fid)[(vInF + 2) % 3];
    Eigen::Vector3d r0 = _mesh.GetVertPos(mesh.GetEdgeVerts(eid0)[1]) - _mesh.GetVertPos(mesh.GetEdgeVerts(eid0)[0]);
    Eigen::Vector3d r1 = _mesh.GetVertPos(mesh.GetEdgeVerts(eid1)[1]) - _mesh.GetVertPos(mesh.GetEdgeVerts(eid1)[0]);

    Eigen::Matrix2d Iinv, I;
    I << r0.dot(r0), r0.dot(r1), r1.dot(r0), r1.dot(r1);
    Iinv = I.inverse();

    Eigen::Vector2d rhs;
    double w1 = omega(eid0);
    double w2 = omega(eid1);
    rhs << w1, w2;

    Eigen::Vector2d u = Iinv * rhs;
    return u[0] * r0 + u[1] * r1;
}

Eigen::Vector3d ComplexLoopZuenko::computeBaryGradThetaFromOmegaPerface(const Eigen::VectorXd& omega, int fid, const Eigen::Vector3d& bary, Mesh &mesh)
{
    Eigen::Vector3d gradTheta = Eigen::Vector3d::Zero();
    for(int i = 0; i < 3; i++)
    {
        gradTheta += bary[i] * computeGradThetaFromOmegaPerface(omega, fid, i, mesh);
    }
    return gradTheta;
}

void ComplexLoopZuenko::updateLoopedZvals(
	const Eigen::VectorXd& omega, const std::vector<std::complex<double>>& zvals, std::vector<std::complex<double>>& upZvals, 
	Mesh &mesh)
{
	int V = mesh.GetVertCount();
	int E = mesh.GetEdgeCount();
	
	upZvals.resize(V + E);
	
	// Even verts
	#pragma omp parallel for schedule(static)
	for (int vi = 0; vi < V; ++vi)
	{
		if (mesh.IsVertBoundary(vi))
		{
			if(_isFixBnd) {
				upZvals[vi] = zvals[vi];
			}
			else
			{
				std::vector<int> boundary(2);
				boundary[0] = mesh.GetVertEdges(vi).front();
				boundary[1] = mesh.GetVertEdges(vi).back();

				std::vector<std::complex<double>> zp(2);
				std::vector<Eigen::Vector3d> gradthetap(2);
				std::vector<double> coords = { 1. / 2, 1. / 2 };
				std::vector<Eigen::Vector3d> pList(2);

				for (int j = 0; j < boundary.size(); ++j)
				{
					int edge = boundary[j];
					int face = mesh.GetEdgeFaces(edge)[0];
					int viInface = mesh.GetVertIndexInFace(face, vi);

					int viInEdge = mesh.GetVertIndexInEdge(edge, vi);
					int vj = mesh.GetEdgeVerts(edge)[(viInEdge + 1) % 2];

					int vjInface = mesh.GetVertIndexInFace(face, vj);

					Eigen::Vector3d bary = Eigen::Vector3d::Zero();
					bary(viInface) = 3. / 4;
					bary(vjInface) = 1. / 4;

                    pList[j] = 3. / 4 * _mesh.GetVertPos(vi) + 1. / 4 * _mesh.GetVertPos(vj);
                    gradthetap[j] = computeGradThetaFromOmegaPerface(omega, face, viInface, mesh);
                    zp[j] = computeZandGradZ(omega, zvals, face, bary, NULL, mesh);
				}
				upZvals[vi] = interpZ(zp, gradthetap, coords, pList);
			}
		}
		else
		{
			const std::vector<int>& vFaces = mesh.GetVertFaces(vi);
			int nNeiFaces = vFaces.size();

			Scalar alpha = 0.375;
			if (nNeiFaces == 3) alpha /= 2;
			else                    alpha /= nNeiFaces;

			double beta = nNeiFaces / 2. * alpha;

			std::vector<std::complex<double>> zp(nNeiFaces);
			std::vector<Eigen::Vector3d> gradthetap(nNeiFaces);
			std::vector<double> coords;
			coords.resize(nNeiFaces, 1. / nNeiFaces);
			std::vector<Eigen::Vector3d> pList(nNeiFaces);

			for (int k = 0; k < nNeiFaces; ++k)
			{
				int face = vFaces[k];
				int viInface = mesh.GetVertIndexInFace(face, vi);
				Eigen::Vector3d bary;
				bary.setConstant(beta);
				bary(viInface) = 1 - 2 * beta;

				pList[k] = Eigen::Vector3d::Zero();
				for (int i = 0; i < 3; i++)
				{
					pList[k] += bary(i) * _mesh.GetVertPos(mesh.GetFaceVerts(face)[i]);
				}

                zp[k] = computeZandGradZ(omega, zvals, face, bary, NULL, mesh);
                gradthetap[k] = computeBaryGradThetaFromOmegaPerface(omega, face, bary, mesh);
			}
			upZvals[vi] = interpZ(zp, gradthetap, coords, pList);
		}
	}

	// Odd verts
	#pragma omp parallel for schedule(static)
	for (int edge = 0; edge < E; ++edge)
	{
		int row = edge + V;
		if (mesh.IsEdgeBoundary(edge))
		{
			int face = mesh.GetEdgeFaces(edge)[0];
			int eindexFace = mesh.GetEdgeIndexInFace(face, edge);
			Eigen::Vector3d bary;
			bary.setConstant(0.5);
			bary((eindexFace + 2) % 3) = 0;
            
			upZvals[row] = computeZandGradZ(omega, zvals, face, bary, NULL, mesh);
		}
		else
		{
			std::vector<std::complex<double>> zp(2);
			std::vector<Eigen::Vector3d> gradthetap(2);
			std::vector<double> coords = { 1. / 2, 1. / 2 };
			std::vector<Eigen::Vector3d> pList(2);

			for (int j = 0; j < 2; ++j)
			{
				int face = mesh.GetEdgeFaces(edge)[j];
				int offset = mesh.GetEdgeIndexInFace(face, edge);

				Eigen::Vector3d bary;
				bary.setConstant(3. / 8.);
				bary((offset + 2) % 3) = 0.25;

				pList[j] = Eigen::Vector3d::Zero();
				for (int i = 0; i < 3; i++)
				{
					pList[j] += bary(i) * _mesh.GetVertPos(mesh.GetFaceVerts(face)[i]);
				}

                zp[j] = computeZandGradZ(omega, zvals, face, bary, NULL, mesh);
                gradthetap[j] = computeBaryGradThetaFromOmegaPerface(omega, face, bary, mesh);
			}

			upZvals[row] = interpZ(zp, gradthetap, coords, pList);
		}
	}
}

void ComplexLoopZuenko::updateLoopedZvals(
	const Eigen::VectorXd& omega1, const std::vector<std::complex<double>>& zvals1, std::vector<std::complex<double>>& upZvals1, 
	const Eigen::VectorXd& omega2, const std::vector<std::complex<double>>& zvals2, std::vector<std::complex<double>>& upZvals2, 
	Mesh &mesh)
{
	int V = mesh.GetVertCount();
	int E = mesh.GetEdgeCount();
	
	upZvals1.resize(V + E);
	upZvals2.resize(V + E);
	
	// Even verts
	#pragma omp parallel for schedule(static)
	for (int vi = 0; vi < V; ++vi)
	{
		if (mesh.IsVertBoundary(vi))
		{
			if(_isFixBnd) {
				upZvals1[vi] = zvals1[vi];
				upZvals2[vi] = zvals2[vi];
			}
			else
			{
				std::vector<int> boundary(2);
				boundary[0] = mesh.GetVertEdges(vi).front();
				boundary[1] = mesh.GetVertEdges(vi).back();

				// For zvals1/omega1
				std::vector<std::complex<double>> zp1(2);
				std::vector<Eigen::Vector3d> gradthetap1(2);
				std::vector<double> coords = { 1. / 2, 1. / 2 };
				std::vector<Eigen::Vector3d> pList(2);

				// For zvals2/omega2
				std::vector<std::complex<double>> zp2(2);
				std::vector<Eigen::Vector3d> gradthetap2(2);

				for (int j = 0; j < boundary.size(); ++j)
				{
					int edge = boundary[j];
					int face = mesh.GetEdgeFaces(edge)[0];
					int viInface = mesh.GetVertIndexInFace(face, vi);

					int viInEdge = mesh.GetVertIndexInEdge(edge, vi);
					int vj = mesh.GetEdgeVerts(edge)[(viInEdge + 1) % 2];

					int vjInface = mesh.GetVertIndexInFace(face, vj);

					Eigen::Vector3d bary = Eigen::Vector3d::Zero();
					bary(viInface) = 3. / 4;
					bary(vjInface) = 1. / 4;

                    pList[j] = 3. / 4 * _mesh.GetVertPos(vi) + 1. / 4 * _mesh.GetVertPos(vj);
                    
                    // For omega1/zvals1
                    gradthetap1[j] = computeGradThetaFromOmegaPerface(omega1, face, viInface, mesh);
                    zp1[j] = computeZandGradZ(omega1, zvals1, face, bary, NULL, mesh);

                    // For omega2/zvals2
                    gradthetap2[j] = computeGradThetaFromOmegaPerface(omega2, face, viInface, mesh);
                    zp2[j] = computeZandGradZ(omega2, zvals2, face, bary, NULL, mesh);
				}
				interpZ(zp1, gradthetap1, zp2, gradthetap2, coords, pList, upZvals1[vi], upZvals2[vi]);
			}
		}
		else
		{
			const std::vector<int>& vFaces = mesh.GetVertFaces(vi);
			int nNeiFaces = vFaces.size();

			// Fig5 left [Wang et al. 2006]
			Scalar alpha = 0.375;
			if (nNeiFaces == 3) alpha /= 2;
			else                    alpha /= nNeiFaces;

			double beta = nNeiFaces / 2. * alpha;

			// For zvals1/omega1
			std::vector<std::complex<double>> zp1(nNeiFaces);
			std::vector<Eigen::Vector3d> gradthetap1(nNeiFaces);
			
			// For zvals2/omega2
			std::vector<std::complex<double>> zp2(nNeiFaces);
			std::vector<Eigen::Vector3d> gradthetap2(nNeiFaces);
			
			std::vector<double> coords;
			coords.resize(nNeiFaces, 1. / nNeiFaces);
			std::vector<Eigen::Vector3d> pList(nNeiFaces);

			for (int k = 0; k < nNeiFaces; ++k)
			{
				int face = vFaces[k];
				int viInface = mesh.GetVertIndexInFace(face, vi);
				Eigen::Vector3d bary;
				bary.setConstant(beta);
				bary(viInface) = 1 - 2 * beta;

				pList[k] = Eigen::Vector3d::Zero();
				for (int i = 0; i < 3; i++)
				{
					pList[k] += bary(i) * _mesh.GetVertPos(mesh.GetFaceVerts(face)[i]);
				}

                // For omega1/zvals1
                zp1[k] = computeZandGradZ(omega1, zvals1, face, bary, NULL, mesh);
                gradthetap1[k] = computeBaryGradThetaFromOmegaPerface(omega1, face, bary, mesh);
            
                // For omega2/zvals2
                zp2[k] = computeZandGradZ(omega2, zvals2, face, bary, NULL, mesh);
                gradthetap2[k] = computeBaryGradThetaFromOmegaPerface(omega2, face, bary, mesh);
			}
			interpZ(zp1, gradthetap1, zp2, gradthetap2, coords, pList, upZvals1[vi], upZvals2[vi]);
		}
	}

	// Odd verts
	#pragma omp parallel for schedule(static)
	for (int edge = 0; edge < E; ++edge)
	{
		int row = edge + V;
		if (mesh.IsEdgeBoundary(edge))
		{
			int face = mesh.GetEdgeFaces(edge)[0];
			int eindexFace = mesh.GetEdgeIndexInFace(face, edge);
			Eigen::Vector3d bary;
			bary.setConstant(0.5);
			bary((eindexFace + 2) % 3) = 0;
            
			upZvals1[row] = computeZandGradZ(omega1, zvals1, face, bary, NULL, mesh);
			upZvals2[row] = computeZandGradZ(omega2, zvals2, face, bary, NULL, mesh);
		}
		else
		{
			// For zvals1/omega1
			std::vector<std::complex<double>> zp1(2);
			std::vector<Eigen::Vector3d> gradthetap1(2);
			
			// For zvals2/omega2
			std::vector<std::complex<double>> zp2(2);
			std::vector<Eigen::Vector3d> gradthetap2(2);
			
			std::vector<double> coords = { 1. / 2, 1. / 2 };
			std::vector<Eigen::Vector3d> pList(2);

			for (int j = 0; j < 2; ++j)
			{
				int face = mesh.GetEdgeFaces(edge)[j];
				int offset = mesh.GetEdgeIndexInFace(face, edge);

				Eigen::Vector3d bary;
				bary.setConstant(3. / 8.);
				bary((offset + 2) % 3) = 0.25;

				pList[j] = Eigen::Vector3d::Zero();
				for (int i = 0; i < 3; i++)
				{
					pList[j] += bary(i) * _mesh.GetVertPos(mesh.GetFaceVerts(face)[i]);
				}

                // For omega1/zvals1
                zp1[j] = computeZandGradZ(omega1, zvals1, face, bary, NULL, mesh);
                gradthetap1[j] = computeBaryGradThetaFromOmegaPerface(omega1, face, bary, mesh);
             
                // For omega2/zvals2
                zp2[j] = computeZandGradZ(omega2, zvals2, face, bary, NULL, mesh);
                gradthetap2[j] = computeBaryGradThetaFromOmegaPerface(omega2, face, bary, mesh);
			}

			interpZ(zp1, gradthetap1, zp2, gradthetap2, coords, pList, upZvals1[row], upZvals2[row]);
		}
	}
}

// void ComplexLoopZuenko::Subdivide(
// 	const Eigen::VectorXd& omega, 
// 	const std::vector<std::complex<double>>& zvals, 
// 	Eigen::VectorXd& omegaNew, 
// 	std::vector<std::complex<double>>& upZvals, 
// 	int level,
// 	std::vector<SparseMatrixX> &S0List, 
// 	std::vector<SparseMatrixX> &S1List,
// 	std::vector<Mesh> &meshList	
// )
// {
// 	static const bool timingEnabled = []() {
// 		const char* v = std::getenv("CWF_TIMING");
// 		return (v != nullptr) && (std::string(v) != "0");
// 	}();

// 	using Clock = std::chrono::high_resolution_clock;
// 	auto totalStart = Clock::now();
// 	int startVerts = _mesh.GetVertCount();
// 	int startEdges = _mesh.GetEdgeCount();
// 	int startFaces = _mesh.GetFaceCount();
	
// 	int nverts = _mesh.GetVertCount();
// 	omegaNew = omega;
// 	upZvals = zvals;


// 	MatrixX X;
// 	_mesh.GetPos(X);

// 	Eigen::VectorXd amp(nverts);
	
// 	for (int i = 0; i < nverts; i++)
// 	{
// 		amp(i) = std::abs(zvals[i]);
// 	}

// 	for (int l = 0; l < level; ++l) 
// 	{
// 		auto iterStart = Clock::now();
// 		auto t0 = iterStart;

// 		auto t1 = Clock::now();
// 		auto t2 = Clock::now();

// 		X = S0List[l] * X;
// 		auto t3 = Clock::now();
// 		amp = S0List[l] * amp;
// 		auto t4 = Clock::now();

// 		std::vector<std::complex<double>> upZvalsNew;
// 		//std::vector<Eigen::Vector3cd> upGradZvals;

// 		updateLoopedZvals(omegaNew, upZvals, upZvalsNew);
// 		auto t5 = Clock::now();

// 		//updateLoopedZvalsNew(meshNew, omegaNew, upZvals, upZvalsNew);
// 		omegaNew = S1List[l] * omegaNew;
// 		auto t6 = Clock::now();

// 		std::vector<Vector3> points;
// 		ConvertToVector3(X, points);
// 		auto t7 = Clock::now();

// 		// std::vector< std::vector<int> > edgeToVert;
// 		// GetSubdividedEdges(edgeToVert);
// 		auto t8 = Clock::now();

// 		// std::vector< std::vector<int> > faceToVert;
// 		// GetSubdividedFaces(faceToVert);
// 		auto t9 = Clock::now();

// 		_mesh.Populate(points, meshList[l+1]);
// 		auto t10 = Clock::now();

// 		upZvals.swap(upZvalsNew);

// 		if (timingEnabled)
// 		{
// 			auto ms = [](const Clock::time_point& a, const Clock::time_point& b) {
// 				return std::chrono::duration_cast<std::chrono::milliseconds>(b - a).count();
// 			};
// 			std::cout << "[CWF_TIMING][ComplexLoopZuenko::Subdivide] level " << l << "/" << level
// 			          << " | BuildS0=" << ms(t0, t1) << "ms"
// 			          << " BuildS1=" << ms(t1, t2) << "ms"
// 			          << " S0*X=" << ms(t2, t3) << "ms"
// 			          << " S0*amp=" << ms(t3, t4) << "ms"
// 			          << " updateLoopedZvals=" << ms(t4, t5) << "ms"
// 			          << " S1*omega=" << ms(t5, t6) << "ms"
// 			          << " ConvertToVector3=" << ms(t6, t7) << "ms"
// 			          << " GetSubdividedEdges=" << ms(t7, t8) << "ms"
// 			          << " GetSubdividedFaces=" << ms(t8, t9) << "ms"
// 			          << " Populate=" << ms(t9, t10) << "ms"
// 			          << " iterTotal=" << ms(iterStart, t10) << "ms"
// 			          << std::endl;
// 		}
// 	}

// 	if (timingEnabled)
// 	{
// 		auto totalEnd = Clock::now();
// 		auto totalMs = std::chrono::duration_cast<std::chrono::milliseconds>(totalEnd - totalStart).count();
// 		int endVerts = _mesh.GetVertCount();
// 		int endEdges = _mesh.GetEdgeCount();
// 		int endFaces = _mesh.GetFaceCount();
// 		std::cout << "[CWF_TIMING][ComplexLoopZuenko::Subdivide] total=" << totalMs << "ms"
// 		          << " levels=" << level
// 		          << " V " << startVerts << "->" << endVerts
// 		          << " E " << startEdges << "->" << endEdges
// 		          << " F " << startFaces << "->" << endFaces
// 		          << std::endl;
// 	}
// }
