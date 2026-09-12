#ifndef REST_SHAPE_HPP
#define REST_SHAPE_HPP

#include "common.hpp"

namespace argus
{
	class Sheet;
	class RestShape
	{
		std::pair<bool, std::vector<Real>> m_RestLengths;
		std::pair<bool, std::vector<Real>> m_CrossLengths;
		std::pair<bool, std::vector<Real>> m_BridsonSineFactors;
		std::pair<bool, std::vector<Matrix32>> m_DeformedCoordsToDeformationGradient;
		std::pair<bool, std::vector<Matrix22>> m_uvDiffs;
		std::pair<bool, std::vector<Matrix22>> m_IToFTF;
		std::pair<bool, std::vector<Real>> m_Areas;
		// Geometry matrices
		// std::pair<bool, std::vector<Matrix33>> m_Normals; // per face rest normals
		std::pair<bool, std::vector<Matrix22>> m_Ibars; // stores per face rest first fundamental form
		std::pair<bool, std::vector<Matrix22>> m_IIbars;
		std::pair<bool, MatrixNN> m_Ts; // Transition matrices. Ts.block<2,2>(2*i,0) maps vectors from barycentric coordinates of face edgeFace(i,0) to barycentric coordinates of edgeFace(i,1). Ts.block<2,2>(2*i, 2) is opposite.
		std::pair<bool, MatrixNN> m_Js; // Js.block<2,2>(2*i,0) rotates vectors on face i (in face i's barycentric coordinates) to the perpendicular vector (under the metric abars[i])
		Sheet& m_Owner;

	public:
		RestShape(Sheet& sheet)
		    : m_Owner(sheet) { };
		void requireRestLengths();
		void requireCrossLengths();
		void requireBridsonSineFactors();
		void requireDeformedCoordsToDeformationGradient();
		void requireAreas();
		void requireIntrinsicGeometry();
		const std::vector<Real>& getRestLengths() const { return m_RestLengths.second; }
		const std::vector<Real>& getCrossLengths() const { return m_CrossLengths.second; }
		const std::vector<Real>& getBridsonSineFactors() { return m_BridsonSineFactors.second; }
		const std::vector<Matrix32>& getDeformedCoordsToDeformationGradient() const { return m_DeformedCoordsToDeformationGradient.second; }
		const std::vector<Matrix22>& getUVDiffs() const { return m_uvDiffs.second; }
		const std::vector<Real>& getAreas() const { return m_Areas.second; }
		const MatrixNN& getTs() const { return m_Ts.second; }
		const MatrixNN& getJs() const { return m_Js.second; }
		const std::pair<bool, std::vector<Matrix22>>& getIToFTF() { return m_IToFTF; }

		const std::vector<Matrix22>& getIbars() const
		{
			if (m_Ibars.first == false)
				throw " error in ibars";
			return m_Ibars.second;
		}
		const std::vector<Matrix22>& getIIbars() const { return m_IIbars.second; }
	};
}

#endif /* REST_SHAPE_HPP */
