#include <argus/wrinkle_reconstruction/phase_reconstruction.hpp>
#include <complex>
#include <vector>
namespace argus
{
	WrinklePhase::WrinklePhase(Sheet& sheet)
	{
		initVBDPhaseReconstruction(sheet);
	}
	Real WrinklePhase::computePhaseReconstructionEnergy(argus::Sheet& sheet, const argus::VectorN& phiold, Real timestep)
	{
		const std::vector<trimesh::triangle_t>& faces = sheet.getTriangles();
		const Matrix3N& positions = sheet.getPositionsN3();
		const Matrix2N& dphis = sheet.getDphis1VectorPerFace();
		const VectorN& amps = sheet.getAmpVec();
		const VectorN& phis = sheet.getPhis();
		const VectorN& faceAreas = sheet.getFaceAreas();
		const VectorN& vertexMass = sheet.getInertia();
		Real beta = sheet.getLameBeta();
		Real alpha = sheet.getLameAlpha();
		Real thickness = sheet.getThickness();
		Eigen::Vector<Real, Eigen::Dynamic> faceEnergies, vertexEnergies;
		faceEnergies.resize(sheet.getFaceCount());
		vertexEnergies.resize(sheet.getVertexCount());
		// phase compatibility energy
		Real Ec = 0;
#ifdef USE_OMP
#pragma omp parallel for
#endif
		for (int faceId = 0; faceId < faces.size(); faceId++)
		{
			faceEnergies(faceId) = 0;
			Vector3 omegas = Vector3({ dphis(0, faceId), -dphis(1, faceId), -dphis(0, faceId) + dphis(1, faceId) });
			Real faceArea = faceAreas(faceId);
			for (int vid = 0; vid < 3; vid++)
			{
				int vidi = faces[faceId].v[vid];
				int vidj = faces[faceId].v[(vid + 1) % 3];
				Real ai = amps(vidi);
				Real aj = amps(vidj);
				Real avgamps = ai + aj;
				if (avgamps < 1e-6)
					continue;
				// edge vid, vid+1%3
				Real eijnormsq = (positions.col(vidj) - positions.col(vidi)).squaredNorm();
				Real lambdaij = beta * thickness * faceArea * (avgamps * avgamps) / (64 * eijnormsq);
				// original orientation
				Real omegaij = omegas(vid);
				Real phii = phis(vidi);
				Real phij = phis(vidj);
				std::complex<Real> eiwij = std::polar(1.0, omegaij);
				faceEnergies(faceId) += lambdaij * std::norm(std::polar(1.0, phij) - std::polar(1.0, omegaij) * std::polar(1.0, phii));
				faceEnergies(faceId) += lambdaij * std::norm(std::polar(1.0, phij) - std::polar(1.0, -omegaij) * std::polar(1.0, phii));
			}
		}
		// inertial term
		Real Ep = 0;
#ifdef USE_OMP
#pragma omp parallel for
#endif
		for (int i = 0; i < sheet.getVertexCount(); i++)
		{

			Real mui = vertexMass(i) * amps(i) * amps(i) / (2.0 * timestep * timestep);
			Real phii = phis(i);
			Real phiiold = phiold(i);
			vertexEnergies(i) = std::norm(std::polar(1.0, phii) - std::polar(1.0, phiiold));
		}
		// accumulate
		Ep = vertexEnergies.sum();
		Ec = faceEnergies.sum();
		return Ep + Ec;
	}
	void WrinklePhase::initVBDPhaseReconstruction(argus::Sheet& sheet) { sheet.requireVertexColoringInfo(); }
	std::complex<Real> WrinklePhase::getPhaseFromFaceEnergy(argus::Sheet& sheet, int vid, Real timestep, const Eigen::Vector<std::complex<Real>, Eigen::Dynamic>& x_current)
	{
		Real retval = 0;
		Real h = sheet.getThickness();
		Real beta = sheet.getLameBeta();
		const Eigen::Vector<Real, Eigen::Dynamic>& faceAreas = sheet.getFaceAreas();
		const Matrix3N& restPositiosn = sheet.getRestPositionsN3();

		const Eigen::VectorXi& faces = sheet.sheetInfo.vertexFacesSet.row(vid);
		const std::vector<trimesh::triangle_t>& triangles = sheet.getTriangles();
		const VectorN& amps = sheet.getAmpVec();
		const VectorN& inertia = sheet.getInertia();
		const Matrix2N& dphis = sheet.getDphis1VectorPerFace();
		const VectorN& phis = sheet.getPhis(); // should contain phin
		unsigned int nfaces = 0;
		std::complex<Real> zphi_inertial = x_current(vid);
		Real mi = inertia(vid);
		Real k1 = timestep > 0 ? mi / (2.0 * timestep * timestep) : 0;

		std::complex<Real> totalZ = k1 * zphi_inertial;
		Real totalWeight = k1;

		for (int i = 0; i < faces.rows(); i++)
		{
			if (faces(i) == -1)
				break;
			nfaces++;
		}
		// #ifdef USE_OMP
		// #pragma omp parallel for
		// #endif
		for (int i = 0; i < nfaces; i++)
		{
			int faceId = faces(i);
			Vector3 omegas = Vector3({ dphis.col(faceId)(0), dphis.col(faceId)(1) - dphis.col(faceId)(0), -dphis.col(faceId)(1) });
			int localVid = -1;
			for (int ii = 0; ii < 3; ii++)
			{
				int vidii = triangles.at(faceId).v[ii];
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
			// now since we have local vertex index,
			// the two incident edges on the vertex will be
			int globali = triangles.at(faceId).v[localVid];

			int locali, localj;
			locali = localVid;
			for (int localedge = 1; localedge < 3; localedge++)
			{
				localj = (locali + localedge) % 3;
				int globalj = triangles.at(faceId).v[localj];

				// now we have successfully identified the edge ij, get the contribution from omega ij
				Real Aij = (1.0 / 3.0) * faceAreas(faceId);
				Real ai = amps(globali);
				Real aj = amps(globalj);
				Real eijnormsquared = (restPositiosn.col(localj) - restPositiosn.col(locali)).squaredNorm();
				Real k2 = beta * h * Aij * (ai + aj) * (ai + aj) / (64 * eijnormsquared);
				Real omegaij = (localedge == 1) ? omegas(localVid) : omegas(localj);
				std::complex<Real> zphij = x_current(globalj);

				// Accumulate complex contribution
				totalZ += k2 * (zphij * std::polar(1.0, -omegaij) + zphij * std::polar(1.0, omegaij));
				totalWeight += 2.0 * k2;
			}
		}
		if (std::abs(totalZ) > 0)
			totalZ /= std::abs(totalZ);
		return totalZ;
	}

	void WrinklePhase::solvePhi(argus::Sheet& sheet, Real timestep)
	{
		Real youngsmodulus = sheet.getYoungsModulus();

		// start vbd like optimization
		const SheetInfo& vbdInfo = sheet.getSheetInfo();
		// get mass information
		VectorN& mass = sheet.getInertia();
		VectorN invInertia = sheet.getInvInertia();

		Eigen::Vector<std::complex<Real>, Eigen::Dynamic> xt = sheet.getzPhis();

		VectorN vphi = sheet.getPhiVelocities();
		VectorN phis = sheet.getPhis();
		Eigen::Vector<std::complex<Real>, Eigen::Dynamic> y;
		y.resize(xt.size());
		for (int i = 0; i < y.size(); i++)
		{
			y(i) = std::polar(1.0, phis(i) + timestep * vphi(i));
		}

		// 3. Initialize our working variable x for the VBD iterations
		Eigen::Vector<std::complex<Real>, Eigen::Dynamic> x = y;

		int niters = 100;
		int ncolors = vbdInfo.coloredSets.rows();

		for (int iter = 0; iter < niters; iter++)
		{
			// Inside the loop
			for (int color = 0; color < ncolors; color++)
			{
				argus::VectorNi vset;

				vset = vbdInfo.coloredSets.row(color);
				// #ifdef USE_OMP
				// #pragma omp parallel for schedule(static)
				// #endif
				for (int vtex = 0; vtex < vset.size(); vtex++)
				{
					int vid = vset(vtex);
					if (vid == -1)
						continue;
					if (invInertia(3 * vid) == 0)
						continue;
					std::complex<Real> xi, yi;
					xi = x(vid);
					yi = y(vid);
					argus::Real mvtex = mass(vid);

					x(vid) = getPhaseFromFaceEnergy(sheet, vid, timestep, x);
				}

				// Eigen::internal::set_is_malloc_allowed(true);
				// update positions
				sheet.getzPhis().noalias() = x;
			}
		}

		x.noalias() = sheet.getzPhis();
		for (int vid = 0; vid < x.size(); vid++)
		{
			Real phi_new = std::arg(x(vid));
			Real phi_old = std::arg(xt(vid));

			Real delta_phi = phi_new - phi_old;

			delta_phi = std::atan2(std::sin(delta_phi), std::cos(delta_phi));

			sheet.getPhiVelocities()(vid) = delta_phi / timestep;

			sheet.getPhis()(vid) = phi_new;
		}
		sheet.getzPhis() = x;
	}
}
