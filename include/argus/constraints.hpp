#ifndef CONSTRAINTS_HPP
#define CONSTRAINTS_HPP

#include "common.hpp"

namespace argus
{
	struct ConstraintMatrix
	{
		trimesh::index_t m_Index;
		Matrix33 m_Matrix;
	};
	struct FixedConstraint
	{
		trimesh::index_t m_Index;
	};
	struct PlaneConstraint
	{
		trimesh::index_t m_Index;
		Vector3 m_NormalizedNormal;
	};
	struct LineConstraint
	{
		trimesh::index_t m_Index;
		Vector3 m_NormalizedDirection;
	};
	// dynamic constraints
	struct DynamicConstraint
	{
		// pass in the dofs that need to be constrained
		std::vector<int> dofs;
		// the function should store thresholds and constraint information
		// input is all sheet dofs (active + inactive) and a std vector containing the active dofs
		// constraint will be applied onto the dofs that are active
		std::function<void(VectorN&, VectorN&, std::vector<int>, Real dt)> applyConstraints;
	};
}

#endif /* CONSTRAINTS_HPP */
