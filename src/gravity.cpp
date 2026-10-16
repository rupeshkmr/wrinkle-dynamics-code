#include <argus/gravity.hpp>

#include <argus/sheet.hpp>

namespace argus
{
#ifdef ARGUS_CHECKPOINT
	void Gravity::saveInfo()
	{
		nlohmann::json data;
		data["Energy"] = "TFW Gravity";
		data["Is external"] = isExternal();
		data["gravity"] = m_Gravity;
		simConf["Energies"].push_back(data);
	}
#endif
	void Gravity::getEnergy(Sheet& sheet, Real& energy, VectorN* deriv, std::vector<Eigen::Triplet<Real> >* hess) const
	{
		const trimesh::index_t vertexCount = sheet.getVertexCount();
		VectorN& inertia = sheet.getInertia();
#ifdef USE_OMP
#pragma omp parallel for reduction(+ : energy)
#endif
		for (trimesh::index_t v = 0; v < vertexCount; v++)
		{
			energy += inertia(v) * m_Gravity.dot(sheet.getPositions().col(v));
		}
		if (deriv)
		{
			(*deriv).resize(3 * vertexCount);
#ifdef USE_OMP
#pragma omp parallel for
#endif
			for (trimesh::index_t v = 0; v < vertexCount; v++)
			{
				(*deriv)(Eigen::seq(3 * v, 3 * v + 2)) = inertia(v) * m_Gravity;
			}
		}
	}

	Real Gravity::getEnergy(Sheet& sheet) const
	{

		Real energy = 0;
		const trimesh::index_t vertexCount = sheet.getVertexCount();
		const VectorN& inertia = sheet.getInertia();
		for (trimesh::index_t v = 0; v < vertexCount; v++)
		{
			energy -= inertia(v) * m_Gravity.dot(sheet.getPositionsN3().col(v));
		}

		return energy;
	}

	void Gravity::addForces(Sheet& sheet,
	    Matrix3N& forces) const
	{
		const trimesh::index_t vertexCount = sheet.getVertexCount();
		const VectorN& inertia = sheet.getInertia();

		for (trimesh::index_t v = 0; v < vertexCount; v++)
		{
			forces.col(v) += inertia(v) * m_Gravity;
		}
	}
}
