#include <argus/tfw_gravity.hpp>

namespace argus
{
#ifdef ARGUS_CHECKPOINT
	void TFWGravity::saveInfo()
	{
		nlohmann::json data;
		data["Energy"] = "TFW Gravity";
		data["Is external"] = isExternal();
		data["gravity"] = g;
		data["Amps affected"] = amp_enabled;
		simConf["Energies"].push_back(data);
	}
#endif
	Real TFWGravity::getEnergy(Sheet& sheet) const
	{
		Real energy = 0.0;
		int nvertices = sheet.getPositionsN3().cols();
		int nedges = sheet.getEdges().size();
		int nfaces = sheet.getTriangles().size();
		VectorN amplitudes = sheet.getAmplitudes();
		Matrix3N positions = sheet.getPositionsN3();
		std::vector<trimesh::triangle_t> faces = sheet.getTriangles();
		VectorN inertias = sheet.getInertia();
		// #ifdef USE_OMP
		// #pragma omp parallel for reduction(+ : energy)
		// #endif
		for (int i = 0; i < nvertices; i++)
			energy += inertias(i) * g.dot(sheet.getPositionsN3().col(i));
		return energy;
	}

	void TFWGravity::getEnergy(Sheet& sheet, Real& energy, VectorN* deriv, std::vector<Eigen::Triplet<Real> >* hess) const
	{
		energy = 0.0;
		int nvertices = sheet.getPositionsN3().cols();
		int nedges = sheet.getEdges().size();
		int nfaces = sheet.getTriangles().size();
		VectorN amplitudes = sheet.getAmplitudes();
		Matrix3N positions = sheet.getPositionsN3();
		std::vector<trimesh::triangle_t> faces = sheet.getTriangles();
		VectorN inertias = sheet.getInertia();
		// #ifdef USE_OMP
		// #pragma omp parallel for reduction(+ : energy)
		// #endif
		for (int i = 0; i < nvertices; i++)
			energy += inertias(i) * g.dot(sheet.getPositionsN3().col(i));
		if (deriv)
		{
			deriv->resize(nvertices + 2 * nfaces + 3 * nvertices);
			deriv->setZero();
			// #ifdef USE_OMP
			// #pragma omp parallel for
			// #endif
			for (int i = 0; i < nvertices; i++)
			{
				(*deriv)(Eigen::seq(nvertices + 2 * nfaces + 3 * i, nvertices + 2 * nfaces + 3 * i + 2)) = g * inertias(i);
			}
			// std::cout << "Deriv \n" << (*deriv) << std::endl;
			// exit(0);
		}
		// adding gravity to amplitudes in normal direction
		if (amp_enabled)
		{
			Real e_amps = 0.0;
			VectorN degrees = sheet.getVerticesDegree();
			for (int i = 0; i < nfaces; i++)
			{
				// gravity will force amplitudes to be small in horizontal triangles
				Vector3 vi = positions.col(faces[i].v[0]);
				Vector3 vj = positions.col(faces[i].v[1]);
				Vector3 vk = positions.col(faces[i].v[2]);
				Vector3 normal = (vj - vi).cross(vk - vi).normalized();
				for (int j = 0; j < 3; j++)
				{
					int vtex_degree = degrees(faces[i].v[j]);
					Real amp = amplitudes(faces[i].v[j]);
					e_amps += amp * (normal.dot(g)) / vtex_degree;
					if (deriv)
					{
						(*deriv)(faces[i].v[j]) += normal.dot(g) / vtex_degree;
					}
				}
			}
			energy += e_amps;
		}

		if (hess)
		{
			(*hess) = std::vector<Eigen::Triplet<Real> >(0);
		}
		// if(hess)
		// {
		//     hess->resize(4*nvertices + 2*nfaces, 4*nvertices + 2*nfaces);
		//     hess->setZero();
		// }
	}
	void TFWGravity::getWrinkleShellEnergyPerFace(Sheet& sheet, Real& energy, int faceId, VectorN* deriv, MatrixNN* hess) const
	{
		energy = 0;
		VectorN derivatives;
		VectorN inertias = sheet.getInertia();
		MatrixNN vertex_degrees = sheet.getVerticesDegree();
		int vi = sheet.getTriangles().at(faceId).v[0];
		int vj = sheet.getTriangles().at(faceId).v[1];
		int vk = sheet.getTriangles().at(faceId).v[2];
		energy += inertias(vi) * g.dot(sheet.getPositionsN3().col(vi)) / vertex_degrees(vi);
		energy += inertias(vj) * g.dot(sheet.getPositionsN3().col(vj)) / vertex_degrees(vj);
		energy += inertias(vk) * g.dot(sheet.getPositionsN3().col(vk)) / vertex_degrees(vk);
		// No need to handle deriv as dphis are not affected by gravity
		if (deriv)
		{
			deriv->resize(3 + 2 + 3 * 3 + 3 * 3);
			deriv->setZero();
			// todo add position derivatives divided by vertex degree
			(*deriv)(Eigen::seq(5, 7)) = inertias(vi) * g / vertex_degrees(vi);
			(*deriv)(Eigen::seq(8, 10)) += inertias(vj) * g / vertex_degrees(vj);
			(*deriv)(Eigen::seq(11, 13)) = inertias(vk) * g / vertex_degrees(vk);
		}
		if (hess)
		{
			hess->resize(23, 23);
			hess->setZero();
		}
	}

	void TFWGravity::addForces(Sheet& sheet,
	    Matrix3N& forces) const
	{
		int nvertices = sheet.getVertexCount();
		const VectorN& inertias = sheet.getInertia();
#ifdef USE_OMP
#pragma omp parallel for
#endif
		for (int i = 0; i < nvertices; i++)
		{
			// for each vertex
			// energy = inertias(i) * positions.col(i).dot(g)
			forces.col(i) += -inertias(i) * g;
		}
	}

	void TFWGravity::addForces(Sheet& sheet, VectorN& forces) const
	{
		int nvertices = sheet.getVertexCount();
		int nfaces = sheet.getFaceCount();
		const Matrix3N& positions = sheet.getPositionsN3();
		const VectorN& inertias = sheet.getInertia();
#ifdef USE_OMP
#pragma omp parallel for
#endif
		for (int i = 0; i < nvertices; i++)
		{
			// for each vertex
			// energy = inertias(i) * positions.col(i).dot(g)
			forces.segment(nvertices + 2 * nfaces + 3 * i, 3) += -inertias(i) * g;
		}
	}
	void TFWGravity::addForcesPerFace(Sheet& sheet, Eigen::Vector<Real, 23>& forces, int faceId) const
	{
		// only required for dphis
	}
}
