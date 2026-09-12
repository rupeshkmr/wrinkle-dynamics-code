#include <argus/rest_shape.hpp>
#include <argus/sheet.hpp>
#include <argus/geometry_functions.hpp>
#include "../../vendor/SecondFundamentalForm/MidedgeAverageFormulation.h"
namespace argus
{
	void RestShape::requireIntrinsicGeometry()
	{
		const std::vector<trimesh::triangle_t>& faces = m_Owner.getTriangles();
		const Matrix3N& restPositions = m_Owner.getPositionsN3();
		const std::vector<trimesh::edge_t>& sheetEdges = m_Owner.getEdges();
		const trimesh::trimesh_t& sheetHalfEdges = m_Owner.getHalfEdges();
		const std::vector<trimesh::index_t>& edgeHalfEdges = sheetHalfEdges.getEdgeHalfEdges();
		int nedges = sheetEdges.size();
		if (!m_Ibars.first)
		{
			m_Ibars.first = true;
			m_Ibars.second.resize(faces.size());
			for (size_t i = 0; i < faces.size(); i++)
			{
				// Updating the first fundamental form for each face
				Matrix33 facePositions;
				// facePositions.col(0) = restPositions.block<3, 1>(0, faces[i].v[0]);
				// facePositions.col(1) = restPositions.block<3, 1>(0, faces[i].v[1]);
				// facePositions.col(2) = restPositions.block<3, 1>(0, faces[i].v[2]);

				facePositions.col(0) = restPositions.col(faces[i].v[0]);
				facePositions.col(1) = restPositions.col(faces[i].v[1]);
				facePositions.col(2) = restPositions.col(faces[i].v[2]);
				m_Ibars.second[i] = firstFundamentalForm(facePositions, false, false, nullptr, nullptr);
			}
		}
		if (!m_IIbars.first)
		{
			m_IIbars.first = true;
			m_IIbars.second.resize(faces.size());
			std::shared_ptr<SecondFundamentalFormDiscretization> sff;
			sff = std::make_shared<MidedgeAverageFormulation>(); // using midedgeavg no extra edge dofs required
			for (int faceId = 0; faceId < faces.size(); faceId++)
			{
				m_IIbars.second[faceId] = sff->secondFundamentalForm(m_Owner.getBaseMesh(), restPositions.transpose(), m_Owner.getEdgeDofs(), faceId, nullptr, nullptr);
			}
		}

		if (!m_Js.first)
		{
			m_Js.first = true;
			m_Js.second.resize(2 * faces.size(), 2);
			for (int i = 0; i < faces.size(); i++)
			{
				Matrix22 Jeuclid;
				Jeuclid << 0, -1,
				    1, 0;
				m_Js.second.block<2, 2>(2 * i, 0) = std::sqrt(m_Ibars.second[i].determinant()) * m_Ibars.second[i].inverse() * Jeuclid;
			}
		}

		if (!m_Ts.first)
		{
			m_Ts.first = true;
			// TODO Finish this
			m_Ts.second.resize(2 * nedges, 4);
			m_Ts.second.setZero();
			for (int i = 0; i < nedges; i++)
			{
				trimesh::index_t halfEdgeIndex = edgeHalfEdges[i];
				trimesh::trimesh_t::halfedge_t halfEdge = sheetHalfEdges.halfedge(halfEdgeIndex);
				// -1 face means boundary HalfEdge.
				if (halfEdge.face == -1 || sheetHalfEdges.halfedge(halfEdge.opposite_he).face == -1)
				{
					continue;
				}
				int face1 = halfEdge.face;
				int face2 = sheetHalfEdges.halfedge(halfEdge.opposite_he).face;

				int vert1 = sheetHalfEdges.halfedge(halfEdge.next_he).to_vertex;
				int vert2 = halfEdge.to_vertex;

				// write edge vert1->vert2 in barycentric coordinates on each face
				Vector2 barys[3] = { { 0, 0 }, { 1, 0 }, { 0, 1 } };
				Vector2 face1e(0, 0);
				Vector2 face2e(0, 0);
				for (int j = 0; j < 3; j++)
				{
					if (vert1 == faces[face1].v[j])
						face1e -= barys[j];
					else if (vert2 == faces[face1].v[j])
						face1e += barys[j];
					if (vert1 == faces[face2].v[j])
						face2e -= barys[j];
					else if (vert2 == faces[face2].v[j])
						face2e += barys[j];
				}

				Vector2 face1eperp = m_Js.second.block<2, 2>(2 * face1, 0) * face1e;
				Vector2 face2eperp = m_Js.second.block<2, 2>(2 * face2, 0) * face2e;
				Matrix22 face1basis;
				face1basis.col(0) = face1e;
				face1basis.col(1) = face1eperp;
				Matrix22 face2basis;
				face2basis.col(0) = face2e;
				face2basis.col(1) = face2eperp;
				m_Ts.second.block<2, 2>(2 * i, 0) = face2basis * face1basis.inverse();
				m_Ts.second.block<2, 2>(2 * i, 2) = face1basis * face2basis.inverse();
			}
		}
		if (!m_IToFTF.first)
		{
			m_IToFTF.first = true;
			const std::vector<trimesh::triangle_t>& triangles = m_Owner.getTriangles();
			const trimesh::index_t triangleCount = triangles.size();
			Matrix32 a {
				{ -1, -1 },
				{ 1, 0 },
				{ 0, 1 }
			};
			m_IToFTF.second.resize(triangleCount);
			Real totalError = 0;
			for (trimesh::index_t t = 0; t < triangleCount; t++)
			{
				trimesh::index_t v0 = triangles[t].v[0];
				trimesh::index_t v1 = triangles[t].v[1];
				trimesh::index_t v2 = triangles[t].v[2];
				Matrix33 x;
				x.col(0) = m_Owner.getPositionsN3().col(v0);
				x.col(1) = m_Owner.getPositionsN3().col(v1);
				x.col(2) = m_Owner.getPositionsN3().col(v2);

				Vector3 eij = m_Owner.getPositionsN3().col(v1) - m_Owner.getPositionsN3().col(v0);
				Vector3 eik = m_Owner.getPositionsN3().col(v2) - m_Owner.getPositionsN3().col(v0);
				Vector3 normal = (eij.cross(eik).normalized());
				Vector3 eijcap = eij.normalized();
				Vector3 bt = (normal.cross(eij)).normalized();
				m_IToFTF.second[t](0, 0) = eij.norm();
				m_IToFTF.second[t](0, 1) = eik.dot(eijcap);
				m_IToFTF.second[t](1, 0) = 0;
				m_IToFTF.second[t](1, 1) = eik.dot(bt);
			}
		}
	}

	void RestShape::requireRestLengths()
	{
		if (!m_RestLengths.first)
		{
			m_RestLengths.first = true;

			const std::vector<trimesh::edge_t>& edges = m_Owner.getEdges();
			const trimesh::index_t edgeCount = edges.size();
			Matrix3N& positions = m_Owner.getPositionsN3();

			m_RestLengths.second.resize(edgeCount);

			for (trimesh::index_t e = 0; e < edgeCount; e++)
			{
				m_RestLengths.second[e] = (positions.col(edges[e].v[0]) - positions.col(edges[e].v[1])).norm();
			}
		}
	}

	void RestShape::requireCrossLengths()
	{
		if (!m_CrossLengths.first)
		{
			m_CrossLengths.first = true;

			const trimesh::trimesh_t& halfEdges = m_Owner.getHalfEdges();
			Matrix3N& positions = m_Owner.getPositionsN3();
			const std::vector<trimesh::edge_t>& edges = m_Owner.getEdges();
			const trimesh::index_t edgeCount = edges.size();

			m_CrossLengths.second.resize(edgeCount);

			const std::vector<trimesh::index_t>& edgeHalfEdges = halfEdges.getEdgeHalfEdges();

			// Iterating over all Edges
			for (trimesh::index_t e = 0; e < edgeCount; e++)
			{
				trimesh::index_t halfEdgeIndex = edgeHalfEdges[e];
				trimesh::trimesh_t::halfedge_t halfEdge = halfEdges.halfedge(halfEdgeIndex);
				// -1 face means boundary HalfEdge.
				if (halfEdge.face == -1 || halfEdges.halfedge(halfEdge.opposite_he).face == -1)
				{
					continue;
				}
				const trimesh::index_t index1 = halfEdges.halfedge(halfEdge.next_he).to_vertex;
				const trimesh::index_t index2 = halfEdges.halfedge(halfEdges.halfedge(halfEdge.opposite_he).next_he).to_vertex;
				m_CrossLengths.second[halfEdge.edge] = (positions.col(index1) - positions.col(index2)).norm();
			}
		}
	}

	void RestShape::requireBridsonSineFactors()
	{
		if (!m_BridsonSineFactors.first)
		{
			m_BridsonSineFactors.first = true;

			const trimesh::trimesh_t& halfEdges = m_Owner.getHalfEdges();
			Matrix3N& positions = m_Owner.getPositionsN3();
			const std::vector<trimesh::edge_t>& edges = m_Owner.getEdges();
			const trimesh::index_t edgeCount = edges.size();

			m_BridsonSineFactors.second.resize(edgeCount);

			const std::vector<trimesh::index_t>& edgeHalfEdges = halfEdges.getEdgeHalfEdges();

			// Iterating over all Edges
			for (trimesh::index_t e = 0; e < edgeCount; e++)
			{
				trimesh::index_t halfEdgeIndex = edgeHalfEdges[e];
				trimesh::trimesh_t::halfedge_t halfEdge = halfEdges.halfedge(halfEdgeIndex);
				// -1 face means boundary HalfEdge.
				if (halfEdge.face == -1 || halfEdges.halfedge(halfEdge.opposite_he).face == -1)
				{
					continue;
				}
				const trimesh::index_t index1 = halfEdges.halfedge(halfEdge.next_he).to_vertex;
				const trimesh::index_t index2 = halfEdges.halfedge(halfEdges.halfedge(halfEdge.opposite_he).next_he).to_vertex;
				const trimesh::index_t index3 = edges[halfEdge.edge].v[0];
				const trimesh::index_t index4 = edges[halfEdge.edge].v[1];

				const Vector3 x1x3 = positions.col(index1) - positions.col(index3);
				const Vector3 x1x4 = positions.col(index1) - positions.col(index4);
				const Vector3 Normal1 = x1x3.cross(x1x4);

				const Vector3 x2x3 = positions.col(index2) - positions.col(index3);
				const Vector3 x2x4 = positions.col(index2) - positions.col(index4);
				const Vector3 Normal2 = x2x4.cross(x2x3);

				const Vector3 Edge = positions.col(index4) - positions.col(index3);

				Real sign = Normal1.normalized().cross(Normal2.normalized()).dot(Edge.normalized()) >= 0 ? 1 : -1;

				Real temp = (1 - Normal1.normalized().dot(Normal2.normalized())) / 2.0f;
				Real sineForceFactor = sign * sqrt(std::max(temp, (Real)0.0f));

				m_BridsonSineFactors.second[halfEdge.edge] = sineForceFactor;
			}
		}
	}
	void RestShape::requireDeformedCoordsToDeformationGradient()
	{
		if (!m_DeformedCoordsToDeformationGradient.first)
		{
			m_DeformedCoordsToDeformationGradient.first = true;
			const std::vector<trimesh::triangle_t>& triangles = m_Owner.getUVTriangles();
			const trimesh::index_t triangleCount = triangles.size();
			Matrix32 a {
				{ -1, -1 },
				{ 1, 0 },
				{ 0, 1 }
			};
			m_uvDiffs.first = true;
			argus::Matrix2N uvs = m_Owner.getUVs();

			m_DeformedCoordsToDeformationGradient.second.resize(triangleCount);
			m_uvDiffs.second.resize(triangleCount);

			// F [X1 - X0, X2 - X0] = [x1 - x0, x2 - x0]
			// [x1 - x0, x2 - x0] = a * [x1, x3, x0]

			for (trimesh::index_t t = 0; t < triangleCount; t++)
			{
				trimesh::index_t v0 = triangles[t].v[0];
				trimesh::index_t v1 = triangles[t].v[1];
				trimesh::index_t v2 = triangles[t].v[2];
				// [X1 - X0, X2 - X0]
				Matrix22 uvDiffs;
				uvDiffs << uvs.col(v1) - uvs.col(v0), uvs.col(v2) - uvs.col(v0);
				m_uvDiffs.second[t] = uvDiffs;

				// F = [x0, x1, x2] B
				m_DeformedCoordsToDeformationGradient.second[t] = a * uvDiffs.inverse();
			}
		}
	}

	/* void RestShape::requireDeformedCoordsToDeformationGradient()
	{
	    if (!m_DeformedCoordsToDeformationGradient.first)
	    {
	        m_DeformedCoordsToDeformationGradient.first = true;
	        m_uvDiffs.first = true;
	        const std::vector<trimesh::triangle_t>& triangles = m_Owner.getTriangles();
	        argus::Matrix2N uvs = m_Owner.getUVs();
	        const trimesh::index_t triangleCount = triangles.size();

	        m_DeformedCoordsToDeformationGradient.second.resize(triangleCount);
	        m_uvDiffs.second.resize(triangleCount);

	        // F [X1 - X0, X2 - X0] = [x1 - x0, x2 - x0]
	        // [x1 - x0, x2 - x0] = a * [x1, x3, x0]
	        Matrix32 a {
	            { -1, -1 },
	            { 1, 0 },
	            { 0, 1 }
	        };
	        for (trimesh::index_t t = 0; t < triangleCount; t++)
	        {
	            trimesh::index_t v0 = triangles[t].v[0];
	            trimesh::index_t v1 = triangles[t].v[1];
	            trimesh::index_t v2 = triangles[t].v[2];
	            // [X1 - X0, X2 - X0]
	            Matrix22 uvDiffs;
	            uvDiffs << uvs.col(v1) - uvs.col(v0), uvs.col(v2) - uvs.col(v0);
	            m_uvDiffs.second[t] = uvDiffs;

	            // F = [x0, x1, x2] B
	            m_DeformedCoordsToDeformationGradient.second[t] = a * uvDiffs.inverse();
	        }
	    }
	} */

	void RestShape::requireAreas()
	{
		if (!m_Areas.first)
		{
			m_Areas.first = true;

			const std::vector<trimesh::triangle_t>& triangles = m_Owner.getTriangles();
			Matrix3N& positions = m_Owner.getPositionsN3();
			const trimesh::index_t triangleCount = triangles.size();

			m_Areas.second.resize(triangleCount);

			for (trimesh::index_t t = 0; t < triangleCount; t++)
			{
				m_Areas.second[t] = 0.5 * (positions.col(triangles[t].v[1]) - positions.col(triangles[t].v[0])).cross(positions.col(triangles[t].v[2]) - positions.col(triangles[t].v[0])).norm();
			}
		}
	}
}
