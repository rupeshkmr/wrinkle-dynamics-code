#include <argus/sheet.hpp>
#include <argus/util.hpp>
#include <argus/geometry_functions.hpp>
#include <argus/graph_coloring.hpp>
#include <stdexcept>
#include <utility>
#include <cassert>
#include <algorithm>
#include <sys/stat.h>
#include <igl/adjacency_list.h>
#include <boost/filesystem.hpp>
#include "../../vendor/SecondFundamentalForm/SecondFundamentalFormDiscretization.h"
#include "../../vendor/SecondFundamentalForm/MidedgeAverageFormulation.h"

#define MAX_ADJ_LEN 32

namespace argus
{
	nlohmann::json simConf;
	std::shared_ptr<SecondFundamentalFormDiscretization> sffobject;

	Sheet::Sheet(const Material& material,
	    std::vector<trimesh::triangle_t>& triangles,
	    argus::Matrix3N& positions,
	    std::vector<trimesh::triangle_t>& uvtriangles,
	    argus::Matrix2N& uvs,
	    const SheetInfo& info)
	    : m_Material(material)
	    , m_Positions(std::move(positions))
	    , m_Triangles(std::move(triangles))
	    , m_UVTriangles(std::move(uvtriangles))
	    , m_UVs(std::move(uvs))
	    , m_UseWrinkleParameters(info.wrinkle_mesh)
	    , m_RestShape(*this)
	{
		sffobject = std::make_shared<MidedgeAverageFormulation>();
		curFrame = 0;
		this->sheetInfo = info;
#ifdef ARGUS_CHECKPOINT
		// for creating checkpoint dir
		time_t t = time(NULL);
		std::stringstream st;
		st << t;
		checkpointDirName = st.str();
#endif
		if (sheetInfo.seam_path != "")
		{
			std::string str;
			std::ifstream file(sheetInfo.seam_path);
			while (std::getline(file, str))
			{
				std::stringstream ss(str);
				std::string ch;
				std::vector<int> pairs;
				while (std::getline(ss, ch, ' '))
				{
					pairs.push_back(std::atoi(ch.c_str()));
				}
				m_SeamVertices.insert(std::pair<int, int>(pairs.at(0), pairs.at(1)));
			}
		}

		m_RestPositions = m_Positions;
		m_VertexCount = m_Positions.cols();
		m_FaceCount = m_Triangles.size();
		trimesh::unordered_edges_from_triangles(m_Triangles.size(), &m_Triangles[0], m_Edges);
		m_EdgeCount = m_Edges.size();
		m_HalfEdges.build(m_VertexCount, m_Triangles.size(), &m_Triangles[0], m_Edges.size(), &m_Edges[0]);
		m_Velocities.resize(3, m_VertexCount);
		m_Forces.resize(3, m_VertexCount);
		m_InvInertia.resize(3, m_VertexCount);
		m_Inertia.resize(m_VertexCount);
		m_Inertia.setZero();
		Real density = m_Material.getDensity();
		m_Velocities.setZero();

		// wrinkle stuff
		m_UseReferenceMesh = false;
		if (m_UseWrinkleParameters)
		{
			m_dPhisPerVertex.resize(3, m_VertexCount);
			// Setup wrinkle dynamics
			m_AmplitudeVelocities.resize(m_VertexCount);
			m_AmplitudeVelocities.setZero();
		}

		// required for quadratic bending model
		m_vertAreas.resize(m_VertexCount);
		std::fill(m_vertAreas.begin(), m_vertAreas.end(), 0.0);

		// fill vertex degrees matrix
		m_VerticesDegree.resize(m_VertexCount, 1);
		m_VerticesDegree.setZero();
		m_faceAreas.resize(m_FaceCount);
		int faceId = 0;
		for (int faceId = 0; faceId < m_FaceCount; faceId++)
		{
			trimesh::triangle_t t = m_Triangles.at(faceId);
			Real area = (m_Positions.col(t.j()) - m_Positions.col(t.i())).cross(m_Positions.col(t.k()) - m_Positions.col(t.i())).norm() / 2;
			m_faceAreas(faceId) = area;
			Real massContrib = area * density;
			m_Inertia(t.i()) += massContrib / 3;
			m_Inertia(t.j()) += massContrib / 3;
			m_Inertia(t.k()) += massContrib / 3;
			// for wrinkle stuff
			m_VerticesDegree(t.i()) += 1;
			m_VerticesDegree(t.j()) += 1;
			m_VerticesDegree(t.k()) += 1;
			m_vertAreas[t.i()] += area / 3.0;
			m_vertAreas[t.j()] += area / 3.0;
			m_vertAreas[t.k()] += area / 3.0;
		}

		for (trimesh::index_t i = 0; i < m_Inertia.size(); i++)
		{
			m_InvInertia.col(i) = Vector3::Constant(1.0 / m_Inertia(i));
		}

		setupQuadraturePoints();
		setupEdgeOppToVertexPerFaces();
		if (m_UseWrinkleParameters)
		{
			m_Amplitudes.resize(m_VertexCount);
			m_Amplitudes.setZero();
			m_ampsPerFace.resize(m_VertexCount);
			m_ampsPerFace.setZero();

			m_Phis.resize(m_VertexCount);
			m_Phis.setZero();

			m_zPhis.resize(m_VertexCount);
			m_zPhis.setZero();

			m_PhiVelocities.resize(m_VertexCount);
			m_PhiVelocities.setZero();

			m_dPhis.resize(m_Edges.size());
			m_dPhis.setZero();
			m_dPhis1VectorPerFace.resize(2, m_Triangles.size());
			m_dPhis1VectorPerFace.setZero();
			m_CornerPhis.resize(3, m_Triangles.size());
			m_CornerPhis.setZero();
		}
		// Filling base mesh information
		m_faceMatrix.resize(m_Triangles.size(), 3);
		for (int i = 0; i < m_Triangles.size(); i++)
		{
			m_faceMatrix(i, 0) = m_Triangles[i].v[0];
			m_faceMatrix(i, 1) = m_Triangles[i].v[1];
			m_faceMatrix(i, 2) = m_Triangles[i].v[2];
		}
		m_baseMesh = MeshConnectivity(m_faceMatrix);

		computeLaplacian();
		// for visualization
		m_strainValuesPerFace.resize(m_Triangles.size());
		m_strainValuesPerFace.setZero();
		// register material
		m_Material.registerRestShapeData(*this);

		if (!sheetInfo.clampedVertices.empty())
		{
			std::vector<std::vector<int>> formatted_clamp_data;

			// Iterate over the flat vector in chunks of 4 [idx, x, y, z]
			for (size_t i = 0; i < sheetInfo.clampedVertices.size(); i += 4)
			{
				if (i + 3 < sheetInfo.clampedVertices.size())
				{
					std::vector<int> nodeData;
					nodeData.push_back(sheetInfo.clampedVertices[i]); // Vertex Index
					nodeData.push_back(sheetInfo.clampedVertices[i + 1]); // Clamp X
					nodeData.push_back(sheetInfo.clampedVertices[i + 2]); // Clamp Y
					nodeData.push_back(sheetInfo.clampedVertices[i + 3]); // Clamp Z
					formatted_clamp_data.push_back(nodeData);
				}
			}

			// Pass the formatted data to the existing registration function
			if (!formatted_clamp_data.empty())
			{
				registerClampedDofs(formatted_clamp_data);
			}
		}
		m_ClampedAnim = sheetInfo.clampedVerticesMotion;

		// self-collisions
		m_adjacencyMat = util::buildAdjacencyMatrix(m_faceMatrix, m_VertexCount);
		m_averageEdgeLengths = 0;
		for (auto it = m_Edges.begin(); it != m_Edges.end(); it++)
		{
			m_averageEdgeLengths += (m_RestPositions.col(it->v[0]) - m_RestPositions.col(it->v[1])).norm();
		}
		m_averageEdgeLengths /= m_Edges.size();
		m_minimum_edge_lengths.resize(m_FaceCount);
		m_dphi_to_gradphis.resize(m_FaceCount);
		for (int i = 0; i < m_FaceCount; i++)
		{
			argus::Matrix32 M;
			argus::Matrix32 A;
			int vi = m_Triangles.at(i).v[0];
			int vj = m_Triangles.at(i).v[1];
			int vk = m_Triangles.at(i).v[2];
			A.col(0) = m_RestPositions.col(vj) - m_RestPositions.col(vi);
			A.col(1) = m_RestPositions.col(vk) - m_RestPositions.col(vi);
			M = A * (A.transpose() * A).inverse();
			m_dphi_to_gradphis[i] = M;
			// set delta values
			Vector2 delta;
			delta(0) = 1 / A.col(0).norm();
			delta(1) = 1 / A.col(1).norm();
			m_minimum_edge_lengths[i] = delta.norm();
		}
		m_cachedI.first = false;
		m_cachedI.second.resize(m_FaceCount);
		m_cachedII.first = false;
		m_cachedII.second.resize(m_FaceCount);
		m_cachedDphiDphiT.first = false;
		m_cachedDphiDphiT.second.resize(m_FaceCount);
		m_cachedDaDphiT.first = false;
		m_cachedDaDphiT.second.resize(m_FaceCount);
		m_cachedDaDaT.first = false;
		m_cachedDaDaT.second.resize(m_FaceCount);
		loadAdjacencyLists();
	}
	void Sheet::loadAdjacencyLists()
	{
		// 2. Initialize Adjacency Lists
		m_vertexFaceAdjacencyList.resize(m_VertexCount);
		m_vertexEdgeAdjacencyList.resize(m_VertexCount);
		// 3. Populate Face Adjacency
		for (int i = 0; i < m_FaceCount; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				int vIdx = m_Triangles[i].v[j];
				m_vertexFaceAdjacencyList[vIdx].push_back(i);
			}
		}

		// 4. Populate Edge Adjacency
		for (int i = 0; i < m_EdgeCount; ++i)
		{
			for (int j = 0; j < 2; ++j)
			{
				int vIdx = m_Edges[i].v[j];
				m_vertexEdgeAdjacencyList[vIdx].push_back(i);
			}
		}
	}

	void Sheet::updateBaseMesh()
	{
		m_baseMesh = MeshConnectivity(m_faceMatrix);
	}

	void Sheet::registerClampedDofs(std::vector<std::vector<int>> clampeddofs)
	{
		m_clampInfo.resize(clampeddofs.size(), 4);
		m_clampedVertices.resize(clampeddofs.size());
		for (int i = 0; i < clampeddofs.size(); i++)
		{
			int v = clampeddofs.at(i).at(0);
			m_ClampedDofs.push_back(v);

			m_clampedVertices(i) = v;
			m_InvInertia.col(v) = Vector3::Zero();
			// amplitude dof
			m_Amplitudes(v) = 0;
			m_clampInfo(i, 0) = v;
			for (int j = 1; j < clampeddofs.at(i).size(); j++)
			{
				m_clampInfo(i, j) = clampeddofs.at(i).at(j);
				if (clampeddofs.at(i).at(j) == 1)
				{
					m_ClampedDofs.push_back(m_VertexCount + 2 * m_FaceCount + 3 * v + j - 1);
				}
			}
		}
	}

	VectorN Sheet::getDofsDot(bool active)
	{
		// Note: Currently we are providing time derivative of dphis as zero
		// just for consistency
		VectorN filteredDofs;
		VectorN fullDofs;
		if (m_UseWrinkleParameters)
		{
			fullDofs.resize(4 * m_VertexCount + 2 * m_FaceCount);
			fullDofs(Eigen::seq(0, m_VertexCount - 1)) = m_AmplitudeVelocities;
			fullDofs(Eigen::seq(m_VertexCount, m_VertexCount + 2 * m_FaceCount - 1)) = VectorN::Zero(2 * m_FaceCount);
			fullDofs(Eigen::seq(m_VertexCount + 2 * m_FaceCount, Eigen::indexing::last)) = m_Velocities.reshaped();
		}
		else
		{
			fullDofs.resize(3 * m_VertexCount);
			fullDofs = m_Velocities.reshaped();
		}
		if (active)
		{
			filteredDofs.resize(m_ActiveDofs.size());
			filteredDofs = fullDofs(m_ActiveDofs);
			return filteredDofs;
		}
		else
			return fullDofs;
	}

	VectorN Sheet::getDofs(bool active)
	{
		VectorN filteredDofs;
		VectorN fullDofs;
		if (m_UseWrinkleParameters)
		{
			fullDofs.resize(4 * m_VertexCount + 2 * m_FaceCount);
			fullDofs(Eigen::seq(0, m_VertexCount - 1)) = m_Amplitudes;
			fullDofs(Eigen::seq(m_VertexCount, m_VertexCount + 2 * m_FaceCount - 1)) = m_dPhis1VectorPerFace.reshaped();
			fullDofs(Eigen::seq(m_VertexCount + 2 * m_FaceCount, Eigen::indexing::last)) = m_Positions.reshaped();
		}
		else
		{
			fullDofs.resize(3 * m_VertexCount);
			fullDofs = m_Positions.reshaped();
		}
		if (active)
		{
			filteredDofs.resize(m_ActiveDofs.size());
			filteredDofs = fullDofs(m_ActiveDofs);
			return filteredDofs;
		}
		else
			return fullDofs;
	}

	// removes clamped dofs from active dofs and inserts into m_ActiveDofs
	void Sheet::setActiveDofs(std::vector<int> activedofs)
	{
		m_ActiveDofs.clear();
		std::sort(activedofs.begin(), activedofs.end());
		std::sort(m_ClampedDofs.begin(), m_ClampedDofs.end());
		std::set_difference(activedofs.begin(), activedofs.end(), m_ClampedDofs.begin(), m_ClampedDofs.end(), std::back_inserter(m_ActiveDofs));
	}

	// TODO: Probelmatic remove this
	void Sheet::updateDofs(VectorN& dofs, VectorN dofsdot)
	{
		bool dofsdot_present = (dofsdot.size() != 0);
		if (m_UseWrinkleParameters)
		{
			int entry = 0;
			for (int dofIdx : m_ActiveDofs)
			{
				// update amps
				if (dofIdx < m_VertexCount)
				{
					m_Amplitudes(dofIdx) = dofs(entry);
					// update vel
					if (dofsdot_present)
						m_AmplitudeVelocities(dofIdx) = dofsdot(entry);
					entry++;
				}
				else if (dofIdx < m_VertexCount + 2 * m_FaceCount)
				{
					int faceId = int((dofIdx - m_VertexCount) / 2);
					int component = int((dofIdx - m_VertexCount) % 2);
					m_dPhis1VectorPerFace(component, faceId) = dofs(entry++);
				}
				else if (dofIdx < 4 * m_VertexCount + 2 * m_FaceCount)
				{
					int vid = int((dofIdx - m_VertexCount - 2 * m_FaceCount) / 3);
					int component = int((dofIdx - m_VertexCount - 2 * m_FaceCount) % 3);
					m_Positions(component, vid) = dofs(entry);
					if (dofsdot_present)
						m_Velocities(component, vid) = dofsdot(entry);
					entry++;
				}
				else
					throw "Error in updating dofs! Active dofs index overflow!";
			}
		}
		else
		{
			int entry = 0;
			for (int dofIdx : m_ActiveDofs)
			{
				if (dofIdx >= 3 * m_VertexCount)
					throw "Error in updating dofs! Active dofs index overflow!";

				int vid = int((dofIdx) / 3);
				int component = int((dofIdx) % 3);
				m_Positions(component, vid) = dofs(entry);
				if (dofsdot_present)
					m_Velocities(component, vid) = dofsdot(entry);
				entry++;
			}
		}
	}

	Real Sheet::getE(VectorN& X, VectorN& deriv, std::vector<Eigen::Triplet<Real>>* hess)
	{
		VectorN dofscopy = getDofs();
		updateDofs(X);
		deriv.resize(X.rows());
		VectorN tempderiv;

		tempderiv.resize(4 * m_VertexCount + 2 * m_FaceCount);
		tempderiv.setZero();
		Real energy = getEnergies(&tempderiv);
		deriv = tempderiv(m_ActiveDofs);
		updateDofs(dofscopy);
		return energy;
	}

	std::vector<argus::Vector3> Sheet::getFaceNormals()
	{
		m_faceNormals.resize(m_Triangles.size());
		for (int i = 0; i < m_Triangles.size(); i++)
		{
			argus::Vector3 vi = m_Positions.col(m_Triangles.at(i).v[0]);
			argus::Vector3 vj = m_Positions.col(m_Triangles.at(i).v[1]);
			argus::Vector3 vk = m_Positions.col(m_Triangles.at(i).v[2]);
			argus::Vector3 n = (vj - vi).cross(vk - vi);
			n = n.normalized();
			m_faceNormals.push_back(n);
		}
		return m_faceNormals;
	}

	void Sheet::computeLaplacian()
	{
		std::vector<Eigen::Vector3i> bnd_edges;
		Eigen::VectorXi newIndex;
		Eigen::MatrixXi bndV;
		Eigen::MatrixXi bndF;
		Eigen::MatrixXi OppBndV;
		Eigen::MatrixXi restF = m_faceMatrix;
		MatrixNN restV = m_Positions.transpose();
		int nverts = m_Positions.cols();
		igl::boundary_facets(restF, bndV, bndF, OppBndV);
		bnd_edges.clear();
		for (int i = 0; i < bndV.rows(); i++)
		{
			bnd_edges.push_back(Eigen::Vector3i(bndV(i, 0), bndV(i, 1), restF(bndF(i), OppBndV(i))));
		}

		newIndex.resize(restV.rows());
		for (int i = 0; i < newIndex.rows(); i++)
		{
			newIndex(i) = i;
		}

		laplacian.resize(nverts, nverts);
		Eigen::SparseMatrix<double> restLaplacian;
		igl::cotmatrix(restV, restF, restLaplacian);
		for (int i = 0; i < bnd_edges.size(); i++)
		{
			int vidx0 = bnd_edges[i](0);
			int vidx1 = bnd_edges[i](1);
			int vidx2 = bnd_edges[i](2);

			Eigen::Vector3d v0 = restV.row(vidx0); // edge vertex
			Eigen::Vector3d v1 = restV.row(vidx1); // edge vertex
			Eigen::Vector3d v2 = restV.row(vidx2); // edge opposite vertex
			double cotan12 = cotan(v0, v1, v2);
			double cotan02 = cotan(v1, v0, v2);
			restLaplacian.coeffRef(vidx0, vidx0) += 0.5 * cotan02;
			restLaplacian.coeffRef(vidx0, vidx1) += 0.5 * cotan12;
			restLaplacian.coeffRef(vidx0, vidx2) -= 0.5 * (cotan02 + cotan12);
			restLaplacian.coeffRef(vidx1, vidx0) += 0.5 * cotan02;
			restLaplacian.coeffRef(vidx1, vidx1) += 0.5 * cotan12;
			restLaplacian.coeffRef(vidx1, vidx2) -= 0.5 * (cotan02 + cotan12);
		}

		std::vector<Eigen::Triplet<double>> lapList;
		for (int k = 0; k < restLaplacian.outerSize(); k++)
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(restLaplacian, k); it; ++it)
			{
				lapList.push_back(Eigen::Triplet<double>(newIndex(it.row()), newIndex(it.col()), it.value()));
			}
		}

		laplacian.setFromTriplets(lapList.begin(), lapList.end());
	}
	void Sheet::buildDampingLaplacian()
	{
		Eigen::SparseMatrix<double> L(m_VertexCount * 3, m_VertexCount * 3);
		std::vector<Eigen::Triplet<double>> triplets;
		Matrix23 A;
		A << -1, 1, 0,
		    -1, 0, 1;
		m_DampingLaplacian.resize(3 * m_VertexCount, 3 * m_VertexCount);
		for (int faceId = 0; faceId < m_FaceCount; ++faceId)
		{
			int vidxi = m_Triangles[faceId].v[0];
			int vidxj = m_Triangles[faceId].v[1];
			int vidxk = m_Triangles[faceId].v[2];
			// compute orthogonal basis for each face based on the rest shape
			Vector3 eij = m_RestPositions.col(vidxj) - m_RestPositions.col(vidxi);
			Vector3 eik = m_RestPositions.col(vidxk) - m_RestPositions.col(vidxi);
			Vector3 n = eij.cross(eik);
			Vector3 u = eij.normalized();
			n.normalize();
			Vector3 bt = n.cross(u);
			bt.normalize();
			Matrix22 Dm;
			Dm << eij.dot(u), eik.dot(u), eij.dot(bt), eik.dot(bt);
			Eigen::Matrix<Real, 6, 9> G;
			Matrix23 DminvTA = Dm.inverse().transpose() * A;
			Matrix33 I = Matrix33::Identity();
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					G.block(i * 3, j * 3, 3, 3) = DminvTA(i, j) * I;
				}
			}
			Real w = getLameBeta() * getThickness() * (eij.cross(eik)).norm() / 4.0;
			MatrixNN Li = w * G.transpose() * G;
			// fill global laplacian matrix
			for (int i = 0; i < 3; i++)
			{
				int v_i = m_Triangles[faceId].v[i];
				for (int j = 0; j < 3; j++)
				{
					int v_j = m_Triangles[faceId].v[j];

					for (int r = 0; r < 3; r++)
					{
						for (int c = 0; c < 3; c++)
						{
							Real value = Li(3 * i + r, 3 * j + c);
							// Only add non-zero entries to keep the sparse matrix clean
							if (std::abs(value) > 1e-15)
							{
								triplets.push_back(Eigen::Triplet<Real>(3 * v_i + r, 3 * v_j + c, value));
							}
						}
					}
				}
			}
			// for (int i = 0; i < 3; i++)
			// {
			// 	int vidxi = m_Triangles[faceId].v[i];
			// 	for (int j = 0; j < 3; j++)
			// 	{
			// 		int vidxj = m_Triangles[faceId].v[j];

			// 		triplets.push_back(Eigen::Triplet(3 * vidxi, 3 * vidxj, Li(3 * i, 3 * j)));
			// 		triplets.push_back(Eigen::Triplet(3 * vidxi + 1, 3 * vidxj + 1, Li(3 * i + 1, 3 * j + 1)));
			// 		triplets.push_back(Eigen::Triplet(3 * vidxi + 2, 3 * vidxj + 2, Li(3 * i + 2, 3 * j + 2)));
			// 	}
			// }
		}
		m_DampingLaplacian.setFromTriplets(triplets.begin(), triplets.end());
		m_DampingLaplacian.makeCompressed();
	}
	void Sheet::setupEdgeOppToVertexPerFaces()
	{
		// std::cout << "Filling edge opposite to vertices for all triangles\n";
		std::vector<trimesh::index_t> edgeHalfEdges = m_HalfEdges.getEdgeHalfEdges();
		m_EdgeOppositeVertices.resize(m_Edges.size(), 2);
		m_EdgeOppositeVertices.setOnes();
		m_EdgeOppositeVertices *= -1;
		m_VtexOppositeEdgesPerFace.resize(m_Triangles.size());
		for (int i = 0; i < m_Triangles.size(); i++)
		{
			m_VtexOppositeEdgesPerFace[i].v[0] = -1;
			m_VtexOppositeEdgesPerFace[i].v[1] = -1;
			m_VtexOppositeEdgesPerFace[i].v[2] = -1;
		}
		for (int i = 0; i < m_Edges.size(); i++)
		{
			trimesh::index_t halfEdgeIndex = edgeHalfEdges[i];
			trimesh::trimesh_t::halfedge_t halfEdge = m_HalfEdges.halfedge(halfEdgeIndex);
			// GOTO face for this edge and store it for its corresponding opposite vertex
			int face = halfEdge.face;
			int oppFace = m_HalfEdges.halfedge(halfEdge.opposite_he).face;
			if (face != -1)
			{
				int oppvertex = m_HalfEdges.halfedge(halfEdge.next_he).to_vertex;
				int vid_face = -1;
				if (m_Triangles[face].v[0] == oppvertex)
					vid_face = 0;
				else if (m_Triangles[face].v[1] == oppvertex)
					vid_face = 1;
				else if (m_Triangles[face].v[2] == oppvertex)
					vid_face = 2;
				else
					std::cout << " Error in intializing vertex edge array !!!!\n";
				assert(m_VtexOppositeEdgesPerFace[face].v[vid_face] == -1);
				if (!(m_VtexOppositeEdgesPerFace[face].v[vid_face] == -1))
					throw std::runtime_error("Error in creating vertex opposite edges");
				m_VtexOppositeEdgesPerFace[face].v[vid_face] = i;
				if (m_EdgeOppositeVertices(i, 0) == -1)
					m_EdgeOppositeVertices(i, 0) = oppvertex;
				else if (m_EdgeOppositeVertices(i, 1) == -1)
					m_EdgeOppositeVertices(i, 1) = oppvertex;
				else
				{
					std::cout << "error in edge opposite vertices\n";
					exit(0);
				}
			}
			if (oppFace != -1)
			{
				halfEdge = m_HalfEdges.halfedge(halfEdge.opposite_he);
				int oppvertex = m_HalfEdges.halfedge(halfEdge.next_he).to_vertex;
				int vid_face = -1;
				if (m_Triangles[oppFace].v[0] == oppvertex)
					vid_face = 0;
				else if (m_Triangles[oppFace].v[1] == oppvertex)
					vid_face = 1;
				else if (m_Triangles[oppFace].v[2] == oppvertex)
					vid_face = 2;
				else
					std::cout << " Error in intializing vertex edge array !!!!\n";
				// assert(m_VtexOppositeEdgesPerFace[oppFace].v[vid_face] == -1);
				if (!(m_VtexOppositeEdgesPerFace[oppFace].v[vid_face] == -1))
					throw std::runtime_error("ERror while filling vertex opp face!");
				m_VtexOppositeEdgesPerFace[oppFace].v[vid_face] = i;
				if (m_EdgeOppositeVertices(i, 0) == -1)
					m_EdgeOppositeVertices(i, 0) = oppvertex;
				else if (m_EdgeOppositeVertices(i, 1) == -1)
					m_EdgeOppositeVertices(i, 1) = oppvertex;
				else
				{
					std::cout << "error in edge opposite vertices\n";
					exit(0);
				}
			}
		}
	}

	void Sheet::setupQuadraturePoints()
	{
		QuadraturePoints point;
		Real x = 1.0 / 6.0;
		Real weight = 1.0 / 3.0;
		Real pos[3] = { x, x, 1 - 2 * x };
		for (int i = 0; i < 3; i++)
		{
			point.u = pos[i];
			point.v = pos[(i + 1) % 3];
			point.weight = weight;
			m_quadPoints.push_back(point);
		}
	}

	Vector2 Sheet::computeDphi1VectorPerFace(int faceId, Matrix22* gradDphi)
	{
		if (gradDphi)
		{
			gradDphi->setIdentity();
		}
		return m_dPhis1VectorPerFace.col(faceId);
	}

	Vector2 Sheet::computeDphi(int faceId, Matrix23* gradDphi)
	{
		Eigen::Vector3i edgeIndices;
		for (int i = 0; i < 3; i++)
		{
			// Edge index of the edge opposite to the vertex i of face faceId
			// edgeIndices(i) = _state.baseMesh.faceEdge(faceId, i);
			edgeIndices(i) = m_VtexOppositeEdgesPerFace[faceId].v[i];
		}

		Real phiu = m_dPhis(edgeIndices(2));
		Real phiv = m_dPhis(edgeIndices(1));
		// std::cout << "phiu " << edgeIndices(2) << " phiv " << edgeIndices(1);
		// exit(0);
		int flagU = 1;
		int flagV = 1;
		// Convention that for each triangle, we are using dphis from edges that form the canonical basis for the face and for continuity, we want the
		// orientations to be same for all the faces
		if (m_Triangles[faceId].v[0] > m_Triangles[faceId].v[1])
		{
			flagU = -1;
		}
		if (m_Triangles[faceId].v[0] > m_Triangles[faceId].v[2])
		{
			flagV = -1;
		}

		phiu *= flagU;
		phiv *= flagV;

		if (gradDphi)
		{
			gradDphi->setZero();
			gradDphi->coeffRef(0, 2) = flagU;
			gradDphi->coeffRef(1, 1) = flagV;
		}
		// std::cout << "Phiu, Phiv" << phiu << " " << phiv << std::endl;
		return Vector2(phiu, phiv);
	}

	Matrix22 Sheet::computeDphiDphi1VectorPerFaceTensorNew(int faceId, int quadId, std::array<Matrix22, 2>* deriv, std::array<Matrix22, 4>* hess)
	{
		Vector2 dphi;
		dphi = computeDphi1VectorPerFace(faceId, NULL);
		if (deriv)
		{
			// Matrix22 grad_dphiT;
			// grad_dphiT(0, 0) = 2.0 * dphi(0);
			// grad_dphiT(0, 1) = dphi(1);
			// grad_dphiT(1, 0) = dphi(1);
			// grad_dphiT(1, 1) = 0.0;

			(*deriv)[0] << 2 * dphi(0), dphi(1), dphi(1), 0;
			// grad_dphiT(0, 0) = 0.0;
			// grad_dphiT(0, 1) = dphi(0);
			// grad_dphiT(1, 0) = dphi(0);
			// grad_dphiT(1, 1) = 2 * dphi(1);
			// deriv->push_back(grad_dphiT);

			(*deriv)[1] << 0, dphi(0), dphi(0), 2 * dphi(1);
		}
		if (hess)
		{
			// 0th entry
			(*hess)[0] << 2, 0, 0, 0;
			// 1st entry
			(*hess)[1] << 0, 1, 1, 0;
			// 2nd entry
			(*hess)[2] << 0, 1, 1, 0;
			// 4th entry
			(*hess)[3] << 0, 0, 0, 2;
		}
		return dphi * dphi.transpose();
	}
	Matrix22 Sheet::computeDphiDphi1VectorPerFaceTensor(int faceId, int quadId, std::vector<Matrix22>* deriv, std::vector<Matrix22>* hess)
	{
		Vector2 dphi;
		dphi = computeDphi1VectorPerFace(faceId, NULL);
		if (deriv)
		{
			Matrix22 grad_dphiT;
			grad_dphiT(0, 0) = 2.0 * dphi(0);
			grad_dphiT(0, 1) = dphi(1);
			grad_dphiT(1, 0) = dphi(1);
			grad_dphiT(1, 1) = 0.0;

			deriv->push_back(grad_dphiT);

			grad_dphiT(0, 0) = 0.0;
			grad_dphiT(0, 1) = dphi(0);
			grad_dphiT(1, 0) = dphi(0);
			grad_dphiT(1, 1) = 2 * dphi(1);
			deriv->push_back(grad_dphiT);
		}
		if (hess)
		{
			Matrix22 tH;
			// 0th entry
			tH << 2, 0, 0, 0;
			hess->push_back(tH);
			// 1st entry
			tH << 0, 1, 1, 0;
			hess->push_back(tH);
			// 2nd entry
			tH << 0, 1, 1, 0;
			hess->push_back(tH);
			// 4th entry
			tH << 0, 0, 0, 2;
			hess->push_back(tH);
		}
		return dphi * dphi.transpose();
	}

	Matrix22 Sheet::computeDphiDphiTensor(int faceId, int quadId, std::vector<Matrix22>* deriv, std::vector<Matrix33>* hessian)
	{
		Vector2 dphi;
		Matrix23 gradDphi;
		dphi = computeDphi(faceId, deriv ? &gradDphi : NULL);

		if (deriv)
		{
			for (int i = 0; i < 3; i++)
			{
				Matrix22 grad_dphidphiT;
				grad_dphidphiT(0, 0) = 2 * dphi(0) * gradDphi(0, i);
				grad_dphidphiT(0, 1) = dphi(0) * gradDphi(1, i) + dphi(1) * gradDphi(0, i);
				grad_dphidphiT(1, 0) = grad_dphidphiT(0, 1);
				grad_dphidphiT(1, 1) = 2 * dphi(1) * gradDphi(1, i);

				deriv->push_back(grad_dphidphiT);
			}
		}

		if (hessian)
		{
			hessian->resize(4);
			hessian->at(0) << 2 * gradDphi(0, 0) * gradDphi(0, 0), 2 * gradDphi(0, 1) * gradDphi(0, 0), 2 * gradDphi(0, 2) * gradDphi(0, 0),
			    2 * gradDphi(0, 0) * gradDphi(0, 1), 2 * gradDphi(0, 1) * gradDphi(0, 1), 2 * gradDphi(0, 2) * gradDphi(0, 1),
			    2 * gradDphi(0, 0) * gradDphi(0, 2), 2 * gradDphi(0, 1) * gradDphi(0, 2), 2 * gradDphi(0, 2) * gradDphi(0, 2);

			hessian->at(1) << gradDphi(0, 0) * gradDphi(1, 0) + gradDphi(1, 0) * gradDphi(0, 0), gradDphi(0, 1) * gradDphi(1, 0) + gradDphi(1, 1) * gradDphi(0, 0), gradDphi(0, 2) * gradDphi(1, 0) + gradDphi(1, 2) * gradDphi(0, 0),
			    gradDphi(0, 0) * gradDphi(1, 1) + gradDphi(1, 0) * gradDphi(0, 1), gradDphi(0, 1) * gradDphi(1, 1) + gradDphi(1, 1) * gradDphi(0, 1), gradDphi(0, 2) * gradDphi(1, 1) + gradDphi(1, 2) * gradDphi(0, 1),
			    gradDphi(0, 0) * gradDphi(1, 2) + gradDphi(1, 0) * gradDphi(0, 2), gradDphi(0, 1) * gradDphi(1, 2) + gradDphi(1, 1) * gradDphi(0, 2), gradDphi(0, 2) * gradDphi(1, 2) + gradDphi(1, 2) * gradDphi(0, 2);

			hessian->at(2) << gradDphi(0, 0) * gradDphi(1, 0) + gradDphi(1, 0) * gradDphi(0, 0), gradDphi(0, 1) * gradDphi(1, 0) + gradDphi(1, 1) * gradDphi(0, 0), gradDphi(0, 2) * gradDphi(1, 0) + gradDphi(1, 2) * gradDphi(0, 0),
			    gradDphi(0, 0) * gradDphi(1, 1) + gradDphi(1, 0) * gradDphi(0, 1), gradDphi(0, 1) * gradDphi(1, 1) + gradDphi(1, 1) * gradDphi(0, 1), gradDphi(0, 2) * gradDphi(1, 1) + gradDphi(1, 2) * gradDphi(0, 1),
			    gradDphi(0, 0) * gradDphi(1, 2) + gradDphi(1, 0) * gradDphi(0, 2), gradDphi(0, 1) * gradDphi(1, 2) + gradDphi(1, 1) * gradDphi(0, 2), gradDphi(0, 2) * gradDphi(1, 2) + gradDphi(1, 2) * gradDphi(0, 2);

			hessian->at(3) << 2 * gradDphi(1, 0) * gradDphi(1, 0), 2 * gradDphi(1, 1) * gradDphi(1, 0), 2 * gradDphi(1, 2) * gradDphi(1, 0),
			    2 * gradDphi(1, 0) * gradDphi(1, 1), 2 * gradDphi(1, 1) * gradDphi(1, 1), 2 * gradDphi(1, 2) * gradDphi(1, 1),
			    2 * gradDphi(1, 0) * gradDphi(1, 2), 2 * gradDphi(1, 1) * gradDphi(1, 2), 2 * gradDphi(1, 2) * gradDphi(1, 2);
		}

		return dphi * dphi.transpose();
	}

	Real Sheet::computeAmplitudesPerFaceFromQuad(int faceId, int quadId, Vector2* da, Real* gradA, Vector2* gradDA, Matrix33* hessianA, std::vector<Matrix33>* hessianDA)
	{
		Real amp = m_ampsPerFace(faceId);
		if (da)
		{
			(*da).setZero();
		}
		if (gradA)
		{
			*gradA = 1.0;
		}
		if (gradDA)
		{
			(*gradDA).setZero();
		}
		if (hessianA)
		{
			hessianA->setZero();
		}
		if (hessianDA)
		{
			hessianDA->resize(3);
			for (int i = 0; i < 3; i++)
			{
				hessianDA->at(i).setZero();
			}
		}
		return amp;
	}

	Real Sheet::computeAmplitudesFromQuadNew(int faceId, int quadId, Vector2* da, Vector3* gradA, Matrix23* gradDA, Matrix33* hessianA, std::array<Matrix33, 2>* hessianDA)
	{
		Vector3 amps;
		for (int i = 0; i < 3; i++)
		{
			// vertices(i) = _state.baseMesh.faceVertex(faceId, i);
			int vid = m_Triangles[faceId].v[i];
			amps(i) = m_Amplitudes(vid);
		}
		Real u = m_quadPoints[quadId].u;
		Real v = m_quadPoints[quadId].v;
		Real a = u * amps(1) + v * amps(2) + (1 - u - v) * amps(0);

		if (da)
		{
			da->coeffRef(0) = amps(1) - amps(0);
			da->coeffRef(1) = amps(2) - amps(0);
		}

		if (gradA)
		{
			gradA->coeffRef(0) = 1 - u - v;
			gradA->coeffRef(1) = u;
			gradA->coeffRef(2) = v;
		}

		if (hessianA)
		{
			hessianA->setZero();
		}

		if (gradDA)
		{
			gradDA->coeffRef(0, 0) = -1.0;
			gradDA->coeffRef(0, 1) = 1.0;
			gradDA->coeffRef(0, 2) = 0;

			gradDA->coeffRef(1, 0) = -1.0;
			gradDA->coeffRef(1, 1) = 0;
			gradDA->coeffRef(1, 2) = 1.0;
		}
		if (hessianDA)
		{
			(*hessianDA)[0].setZero();
			(*hessianDA)[1].setZero();
		}

		return a;
	}
	Vector2 Sheet::computeDaFromFace(int faceId, Matrix23* gradDA)
	{
		Vector3 amps;
		for (int i = 0; i < 3; i++)
		{
			// vertices(i) = _state.baseMesh.faceVertex(faceId, i);
			int vid = m_Triangles[faceId].v[i];
			amps(i) = m_Amplitudes(vid);
		}
		Vector2 da;
		da(0) = amps(1) - amps(0);
		da(1) = amps(2) - amps(0);
		if (gradDA)
		{
			gradDA->coeffRef(0, 0) = -1.0;
			gradDA->coeffRef(0, 1) = 1.0;
			gradDA->coeffRef(0, 2) = 0;

			gradDA->coeffRef(1, 0) = -1.0;
			gradDA->coeffRef(1, 1) = 0;
			gradDA->coeffRef(1, 2) = 1.0;
		}
		return da;
	}
	Real Sheet::computeAmplitudesFromQuad(int faceId, int quadId, Vector2* da, Vector3* gradA, Matrix23* gradDA, Matrix33* hessianA, std::vector<Matrix33>* hessianDA)
	{
		Vector3 amps;
		for (int i = 0; i < 3; i++)
		{
			// vertices(i) = _state.baseMesh.faceVertex(faceId, i);
			int vid = m_Triangles[faceId].v[i];
			amps(i) = m_Amplitudes(vid);
		}
		Real u = m_quadPoints[quadId].u;
		Real v = m_quadPoints[quadId].v;
		Real a = u * amps(1) + v * amps(2) + (1 - u - v) * amps(0);

		if (da)
		{
			da->coeffRef(0) = amps(1) - amps(0);
			da->coeffRef(1) = amps(2) - amps(0);
		}

		if (gradA)
		{
			gradA->coeffRef(0) = 1 - u - v;
			gradA->coeffRef(1) = u;
			gradA->coeffRef(2) = v;
		}

		if (hessianA)
		{
			hessianA->setZero();
		}

		if (gradDA)
		{
			gradDA->coeffRef(0, 0) = -1.0;
			gradDA->coeffRef(0, 1) = 1.0;
			gradDA->coeffRef(0, 2) = 0;

			gradDA->coeffRef(1, 0) = -1.0;
			gradDA->coeffRef(1, 1) = 0;
			gradDA->coeffRef(1, 2) = 1.0;
		}
		if (hessianDA)
		{
			hessianDA->resize(2);
			hessianDA->at(0).setZero();
			hessianDA->at(1).setZero();
		}

		return a;
	}
	Matrix22 Sheet::computeDaDaTensorNew(int faceId, int quadId, std::array<Matrix22, 3>* deriv, std::array<Matrix33, 4>* hessian)
	{
		Eigen::Vector2d da;
		Eigen::Vector3d gradA;
		Eigen::Matrix<Real, 2, 3> gradDA;
		Eigen::Matrix<Real, 3, 3> hessA;
		std::vector<Eigen::Matrix<Real, 3, 3>> hessDA;

		Real a = computeAmplitudesFromQuad(faceId, quadId, &da, NULL, deriv ? &gradDA : NULL, hessian ? &hessA : NULL, hessian ? &hessDA : NULL);

		if (deriv)
		{
			for (int i = 0; i < 3; i++)
			{
				Matrix22 grad_dadaT;
				grad_dadaT(0, 0) = 2 * da(0) * gradDA(0, i);
				grad_dadaT(0, 1) = da(0) * gradDA(1, i) + da(1) * gradDA(0, i);
				grad_dadaT(1, 0) = grad_dadaT(0, 1);
				grad_dadaT(1, 1) = 2 * da(1) * gradDA(1, i);

				(*deriv)[i].noalias() = grad_dadaT;
			}
		}

		if (hessian)
		{
			// M00
			(*hessian)[0] << 2, -2, 0,
			    -2, 2, 0,
			    0, 0, 0;
			(*hessian)[1] << 2, -1, -1,
			    -1, 0, 1,
			    -1, 1, 0;

			(*hessian)[2] << 2, -1, -1,
			    -1, 0, 1,
			    -1, 1, 0;

			(*hessian)[3] << 2, 0, -2,
			    0, 0, 0,
			    -2, 0, 2;

			/* hessian->at(0) <<
			    2 * gradDA(0, 0) * gradDA(0, 0) + 2.0 * da(0) * hessDA[0](0, 0),
			    2 * gradDA(0, 1) * gradDA(0, 0) + 2.0 * da(0) * hessDA[0](0, 1),
			    2 * gradDA(0, 2) * gradDA(0, 0) + 2.0 * da(0) * hessDA[0](0, 2),
			    2 * gradDA(0, 0) * gradDA(0, 1) + 2.0 * da(0) * hessDA[0](1, 0),
			    2 * gradDA(0, 1) * gradDA(0, 1) + 2.0 * da(0) * hessDA[0](1, 1),
			    2 * gradDA(0, 2) * gradDA(0, 1) + 2.0 * da(0) * hessDA[0](1, 2),
			    2 * gradDA(0, 0) * gradDA(0, 2) + 2.0 * da(0) * hessDA[0](2, 0),
			    2 * gradDA(0, 1) * gradDA(0, 2) + 2.0 * da(0) * hessDA[0](2, 1),
			    2 * gradDA(0, 2) * gradDA(0, 2) + 2.0 * da(0) * hessDA[0](2, 2);

			hessian->at(1) <<
			    gradDA(0, 0) * gradDA(1, 0) + gradDA(1, 0) * gradDA(0, 0) + da(0) * hessDA[1](0, 0) + da(1) * hessDA[0](0, 0), gradDA(0, 1)* gradDA(1, 0) + gradDA(1, 1) * gradDA(0, 0) + da(0) * hessDA[1](0, 1) + da(1) * hessDA[0](0, 1), gradDA(0, 2)* gradDA(1, 0) + gradDA(1, 2) * gradDA(0, 0) + da(0) * hessDA[1](0, 2) + da(1) * hessDA[0](0, 2),
			    gradDA(0, 0)* gradDA(1, 1) + gradDA(1, 0) * gradDA(0, 1) + da(0) * hessDA[1](1, 0) + da(1) * hessDA[0](1, 0), gradDA(0, 1)* gradDA(1, 1) + gradDA(1, 1) * gradDA(0, 1) + da(0) * hessDA[1](1, 1) + da(1) * hessDA[0](1, 1), gradDA(0, 2)* gradDA(1, 1) + gradDA(1, 2) * gradDA(0, 1) + da(0) * hessDA[1](1, 2) + da(1) * hessDA[0](1, 2),
			    gradDA(0, 0)* gradDA(1, 2) + gradDA(1, 0) * gradDA(0, 2) + da(0) * hessDA[1](2, 0) + da(1) * hessDA[0](2, 0), gradDA(0, 1)* gradDA(1, 2) + gradDA(1, 1) * gradDA(0, 2) + da(0) * hessDA[1](2, 1) + da(1) * hessDA[0](2, 1), gradDA(0, 2)* gradDA(1, 2) + gradDA(1, 2) * gradDA(0, 2) + da(0) * hessDA[1](2, 2) + da(1) * hessDA[0](2, 2);

			hessian->at(2) <<
			    gradDA(0, 0) * gradDA(1, 0) + gradDA(1, 0) * gradDA(0, 0) + da(0) * hessDA[1](0, 0) + da(1) * hessDA[0](0, 0), gradDA(0, 1)* gradDA(1, 0) + gradDA(1, 1) * gradDA(0, 0) + da(0) * hessDA[1](0, 1) + da(1) * hessDA[0](0, 1), gradDA(0, 2)* gradDA(1, 0) + gradDA(1, 2) * gradDA(0, 0) + da(0) * hessDA[1](0, 2) + da(1) * hessDA[0](0, 2),
			    gradDA(0, 0)* gradDA(1, 1) + gradDA(1, 0) * gradDA(0, 1) + da(0) * hessDA[1](1, 0) + da(1) * hessDA[0](1, 0), gradDA(0, 1)* gradDA(1, 1) + gradDA(1, 1) * gradDA(0, 1) + da(0) * hessDA[1](1, 1) + da(1) * hessDA[0](1, 1), gradDA(0, 2)* gradDA(1, 1) + gradDA(1, 2) * gradDA(0, 1) + da(0) * hessDA[1](1, 2) + da(1) * hessDA[0](1, 2),
			    gradDA(0, 0)* gradDA(1, 2) + gradDA(1, 0) * gradDA(0, 2) + da(0) * hessDA[1](2, 0) + da(1) * hessDA[0](2, 0), gradDA(0, 1)* gradDA(1, 2) + gradDA(1, 1) * gradDA(0, 2) + da(0) * hessDA[1](2, 1) + da(1) * hessDA[0](2, 1), gradDA(0, 2)* gradDA(1, 2) + gradDA(1, 2) * gradDA(0, 2) + da(0) * hessDA[1](2, 2) + da(1) * hessDA[0](2, 2);

			hessian->at(3) <<
			    2 * gradDA(1, 0) * gradDA(1, 0) + 2.0 * da(1) * hessDA[1](0, 0), 2 * gradDA(1, 1) * gradDA(1, 0) + 2.0 * da(1) * hessDA[1](0, 1), 2 * gradDA(1, 2) * gradDA(1, 0) + 2.0 * da(1) * hessDA[1](0, 2),
			    2 * gradDA(1, 0) * gradDA(1, 1) + 2.0 * da(1) * hessDA[1](1, 0), 2 * gradDA(1, 1) * gradDA(1, 1) + 2.0 * da(1) * hessDA[1](1, 1), 2 * gradDA(1, 2) * gradDA(1, 1) + 2.0 * da(1) * hessDA[1](1, 2),
			    2 * gradDA(1, 0) * gradDA(1, 2) + 2.0 * da(1) * hessDA[1](2, 0), 2 * gradDA(1, 1) * gradDA(1, 2) + 2.0 * da(1) * hessDA[1](2, 1), 2 * gradDA(1, 2) * gradDA(1, 2) + 2.0 * da(1) * hessDA[1](2, 2); */
		}

		return da * da.transpose();
	}
	Matrix22 Sheet::computeDaDaTensor(int faceId, int quadId, std::vector<Matrix22>* deriv, std::vector<Matrix33>* hessian)
	{
		Eigen::Vector2d da;
		Eigen::Vector3d gradA;
		Eigen::Matrix<Real, 2, 3> gradDA;
		Eigen::Matrix<Real, 3, 3> hessA;
		std::vector<Eigen::Matrix<Real, 3, 3>> hessDA;

		Real a = computeAmplitudesFromQuad(faceId, quadId, &da, NULL, deriv ? &gradDA : NULL, hessian ? &hessA : NULL, hessian ? &hessDA : NULL);

		if (deriv)
		{
			for (int i = 0; i < 3; i++)
			{
				Matrix22 grad_dadaT;
				grad_dadaT(0, 0) = 2 * da(0) * gradDA(0, i);
				grad_dadaT(0, 1) = da(0) * gradDA(1, i) + da(1) * gradDA(0, i);
				grad_dadaT(1, 0) = grad_dadaT(0, 1);
				grad_dadaT(1, 1) = 2 * da(1) * gradDA(1, i);

				deriv->push_back(grad_dadaT);
			}
		}

		if (hessian)
		{
			hessian->resize(4);
			// M00
			Matrix33 temp;
			temp << 2, -2, 0,
			    -2, 2, 0,
			    0, 0, 0;
			hessian->at(0) = temp;
			temp << 2, -1, -1,
			    -1, 0, 1,
			    -1, 1, 0;
			hessian->at(1) = temp;
			temp << 2, -1, -1,
			    -1, 0, 1,
			    -1, 1, 0;
			hessian->at(2) = temp;
			temp << 2, 0, -2,
			    0, 0, 0,
			    -2, 0, 2;
			hessian->at(3) = temp;
			/* hessian->at(0) <<
			    2 * gradDA(0, 0) * gradDA(0, 0) + 2.0 * da(0) * hessDA[0](0, 0),
			    2 * gradDA(0, 1) * gradDA(0, 0) + 2.0 * da(0) * hessDA[0](0, 1),
			    2 * gradDA(0, 2) * gradDA(0, 0) + 2.0 * da(0) * hessDA[0](0, 2),
			    2 * gradDA(0, 0) * gradDA(0, 1) + 2.0 * da(0) * hessDA[0](1, 0),
			    2 * gradDA(0, 1) * gradDA(0, 1) + 2.0 * da(0) * hessDA[0](1, 1),
			    2 * gradDA(0, 2) * gradDA(0, 1) + 2.0 * da(0) * hessDA[0](1, 2),
			    2 * gradDA(0, 0) * gradDA(0, 2) + 2.0 * da(0) * hessDA[0](2, 0),
			    2 * gradDA(0, 1) * gradDA(0, 2) + 2.0 * da(0) * hessDA[0](2, 1),
			    2 * gradDA(0, 2) * gradDA(0, 2) + 2.0 * da(0) * hessDA[0](2, 2);

			hessian->at(1) <<
			    gradDA(0, 0) * gradDA(1, 0) + gradDA(1, 0) * gradDA(0, 0) + da(0) * hessDA[1](0, 0) + da(1) * hessDA[0](0, 0), gradDA(0, 1)* gradDA(1, 0) + gradDA(1, 1) * gradDA(0, 0) + da(0) * hessDA[1](0, 1) + da(1) * hessDA[0](0, 1), gradDA(0, 2)* gradDA(1, 0) + gradDA(1, 2) * gradDA(0, 0) + da(0) * hessDA[1](0, 2) + da(1) * hessDA[0](0, 2),
			    gradDA(0, 0)* gradDA(1, 1) + gradDA(1, 0) * gradDA(0, 1) + da(0) * hessDA[1](1, 0) + da(1) * hessDA[0](1, 0), gradDA(0, 1)* gradDA(1, 1) + gradDA(1, 1) * gradDA(0, 1) + da(0) * hessDA[1](1, 1) + da(1) * hessDA[0](1, 1), gradDA(0, 2)* gradDA(1, 1) + gradDA(1, 2) * gradDA(0, 1) + da(0) * hessDA[1](1, 2) + da(1) * hessDA[0](1, 2),
			    gradDA(0, 0)* gradDA(1, 2) + gradDA(1, 0) * gradDA(0, 2) + da(0) * hessDA[1](2, 0) + da(1) * hessDA[0](2, 0), gradDA(0, 1)* gradDA(1, 2) + gradDA(1, 1) * gradDA(0, 2) + da(0) * hessDA[1](2, 1) + da(1) * hessDA[0](2, 1), gradDA(0, 2)* gradDA(1, 2) + gradDA(1, 2) * gradDA(0, 2) + da(0) * hessDA[1](2, 2) + da(1) * hessDA[0](2, 2);

			hessian->at(2) <<
			    gradDA(0, 0) * gradDA(1, 0) + gradDA(1, 0) * gradDA(0, 0) + da(0) * hessDA[1](0, 0) + da(1) * hessDA[0](0, 0), gradDA(0, 1)* gradDA(1, 0) + gradDA(1, 1) * gradDA(0, 0) + da(0) * hessDA[1](0, 1) + da(1) * hessDA[0](0, 1), gradDA(0, 2)* gradDA(1, 0) + gradDA(1, 2) * gradDA(0, 0) + da(0) * hessDA[1](0, 2) + da(1) * hessDA[0](0, 2),
			    gradDA(0, 0)* gradDA(1, 1) + gradDA(1, 0) * gradDA(0, 1) + da(0) * hessDA[1](1, 0) + da(1) * hessDA[0](1, 0), gradDA(0, 1)* gradDA(1, 1) + gradDA(1, 1) * gradDA(0, 1) + da(0) * hessDA[1](1, 1) + da(1) * hessDA[0](1, 1), gradDA(0, 2)* gradDA(1, 1) + gradDA(1, 2) * gradDA(0, 1) + da(0) * hessDA[1](1, 2) + da(1) * hessDA[0](1, 2),
			    gradDA(0, 0)* gradDA(1, 2) + gradDA(1, 0) * gradDA(0, 2) + da(0) * hessDA[1](2, 0) + da(1) * hessDA[0](2, 0), gradDA(0, 1)* gradDA(1, 2) + gradDA(1, 1) * gradDA(0, 2) + da(0) * hessDA[1](2, 1) + da(1) * hessDA[0](2, 1), gradDA(0, 2)* gradDA(1, 2) + gradDA(1, 2) * gradDA(0, 2) + da(0) * hessDA[1](2, 2) + da(1) * hessDA[0](2, 2);

			hessian->at(3) <<
			    2 * gradDA(1, 0) * gradDA(1, 0) + 2.0 * da(1) * hessDA[1](0, 0), 2 * gradDA(1, 1) * gradDA(1, 0) + 2.0 * da(1) * hessDA[1](0, 1), 2 * gradDA(1, 2) * gradDA(1, 0) + 2.0 * da(1) * hessDA[1](0, 2),
			    2 * gradDA(1, 0) * gradDA(1, 1) + 2.0 * da(1) * hessDA[1](1, 0), 2 * gradDA(1, 1) * gradDA(1, 1) + 2.0 * da(1) * hessDA[1](1, 1), 2 * gradDA(1, 2) * gradDA(1, 1) + 2.0 * da(1) * hessDA[1](1, 2),
			    2 * gradDA(1, 0) * gradDA(1, 2) + 2.0 * da(1) * hessDA[1](2, 0), 2 * gradDA(1, 1) * gradDA(1, 2) + 2.0 * da(1) * hessDA[1](2, 1), 2 * gradDA(1, 2) * gradDA(1, 2) + 2.0 * da(1) * hessDA[1](2, 2); */
		}

		return da * da.transpose();
	}

	// void Sheet::updateFirstFundamentalForms()
	// {
	// 	m_IValues.resize(m_Triangles.size());
	// 	for (size_t i = 0; i < m_Triangles.size(); i++)
	// 	{
	// 		// Updating the first fundamental form for each face
	// 		Matrix33 facePositions;
	// 		facePositions.col(0) = m_Positions.block<3, 1>(0,m_Triangles[i].v[0]);
	// 		facePositions.col(1) = m_Positions.block<3, 1>(0,m_Triangles[i].v[1]);
	// 		facePositions.col(2) = m_Positions.block<3, 1>(0,m_Triangles[i].v[2]);
	// 		m_IValues[i] = firstFundamentalForm(facePositions, false, false, nullptr, nullptr);
	// 	}
	// }
	// void Sheet::updateSecondFundamentalForms()
	// {
	// 	std::shared_ptr<SecondFundamentalFormDiscretization> sff;
	// 	sff = std::make_shared<MidedgeAverageFormulation>();
	// }

	Matrix22 Sheet::computeDaDphi1VectorPerFaceTensorNew(int faceId, int quadId, std::array<Matrix22, 5>* deriv, std::array<Eigen::Matrix<Real, 5, 5>, 4>* hess)
	{
		Vector2 da, dphi;
		Vector3 gradA;
		Matrix23 gradDA, gradDphi;
		Matrix33 hessA;
		std::vector<Matrix33> hessDA;

		Real a = computeAmplitudesFromQuadNew(faceId, quadId, &da, NULL, deriv ? &gradDA : NULL, NULL, NULL);
		dphi = computeDphi1VectorPerFace(faceId, NULL);
		// filling derivatives wrt a
		if (deriv)
		{

			for (int i = 0; i < 3; i++)
			{
				(*deriv)[i](0, 0) = gradDA(0, i) * dphi(0);
				(*deriv)[i](0, 1) = gradDA(0, i) * dphi(1);
				(*deriv)[i](1, 0) = gradDA(1, i) * dphi(0);
				(*deriv)[i](1, 1) = gradDA(1, i) * dphi(1);
			}
			// filling derivatives wrt dphi
			Matrix22 grad_dadphiT;
			(*deriv)[3](0, 0) = da(0);
			(*deriv)[3](0, 1) = 0.0;
			(*deriv)[3](1, 0) = da(1);
			(*deriv)[3](1, 1) = 0.0;
			(*deriv)[4](0, 1) = da(0);
			(*deriv)[4](0, 0) = 0.0;
			(*deriv)[4](1, 1) = da(1);
			(*deriv)[4](1, 0) = 0.0;
		}
		if (hess)
		{
			// Eigen::Matrix<Real, 5, 5> hM;
			// hM.setZero();
			// // 00
			// hM(0, 3) = -1;
			// hM(3, 0) = -1;
			// hM(1, 3) = 1;
			// hM(3, 1) = 1;
			// hess->push_back(hM);
			// hM.setZero();
			// // 01
			// hM(0, 4) = -1;
			// hM(4, 0) = -1;
			// hM(1, 4) = 1;
			// hM(4, 1) = 1;
			// hess->push_back(hM);
			// // 10
			// hM(0, 3) = -1;
			// hM(3, 0) = -1;
			// hM(2, 3) = 1;
			// hM(3, 2) = 1;
			// hess->push_back(hM);
			// // 11
			// hM(0, 4) = -1;
			// hM(4, 0) = -1;
			// hM(2, 4) = 1;
			// hM(4, 2) = 1;
			// hess->push_back(hM);
		}
		return da * dphi.transpose();
	}
	Matrix22 Sheet::computeDaDphi1VectorPerFaceTensor(int faceId, int quadId, std::vector<Matrix22>* deriv, std::vector<Eigen::Matrix<Real, 5, 5>>* hess)
	{
		Vector2 da, dphi;
		Vector3 gradA;
		Matrix23 gradDA, gradDphi;
		Matrix33 hessA;
		std::vector<Matrix33> hessDA;

		Real a = computeAmplitudesFromQuad(faceId, quadId, &da, NULL, deriv ? &gradDA : NULL, NULL, NULL);
		dphi = computeDphi1VectorPerFace(faceId, NULL);
		// filling derivatives wrt a
		if (deriv)
		{
			for (int i = 0; i < 3; i++)
			{
				Matrix22 grad_dadphiT;
				grad_dadphiT(0, 0) = gradDA(0, i) * dphi(0);
				grad_dadphiT(0, 1) = gradDA(0, i) * dphi(1);
				grad_dadphiT(1, 0) = gradDA(1, i) * dphi(0);
				grad_dadphiT(1, 1) = gradDA(1, i) * dphi(1);

				deriv->push_back(grad_dadphiT);
			}
			// filling derivatives wrt dphi
			Matrix22 grad_dadphiT;
			grad_dadphiT(0, 0) = da(0);
			grad_dadphiT(0, 1) = 0.0;
			grad_dadphiT(1, 0) = da(1);
			grad_dadphiT(1, 1) = 0.0;

			deriv->push_back(grad_dadphiT);
			grad_dadphiT(0, 1) = da(0);
			grad_dadphiT(0, 0) = 0.0;
			grad_dadphiT(1, 1) = da(1);
			grad_dadphiT(1, 0) = 0.0;

			deriv->push_back(grad_dadphiT);
		}
		if (hess)
		{
			Eigen::Matrix<Real, 5, 5> hM;
			hM.setZero();
			// 00
			hM(0, 3) = -1;
			hM(3, 0) = -1;
			hM(1, 3) = 1;
			hM(3, 1) = 1;
			hess->push_back(hM);
			hM.setZero();
			// 01
			hM(0, 4) = -1;
			hM(4, 0) = -1;
			hM(1, 4) = 1;
			hM(4, 1) = 1;
			hess->push_back(hM);
			// 10
			hM(0, 3) = -1;
			hM(3, 0) = -1;
			hM(2, 3) = 1;
			hM(3, 2) = 1;
			hess->push_back(hM);
			// 11
			hM(0, 4) = -1;
			hM(4, 0) = -1;
			hM(2, 4) = 1;
			hM(4, 2) = 1;
			hess->push_back(hM);
		}
		return da * dphi.transpose();
	}

	Matrix22 Sheet::computeDaDphiTensor(int faceId, int quadId, std::vector<Matrix22>* deriv, std::vector<Eigen::Matrix<Real, 6, 6>>* hessian)
	{
		Vector2 da, dphi;
		Vector3 gradA;
		Matrix23 gradDA, gradDphi;
		Matrix33 hessA;
		std::vector<Matrix33> hessDA;

		Real a = computeAmplitudesFromQuad(faceId, quadId, &da, NULL, deriv ? &gradDA : NULL, hessian ? &hessA : NULL, hessian ? &hessDA : NULL);
		dphi = computeDphi(faceId, deriv ? &gradDphi : NULL);

		if (deriv)
		{
			for (int i = 0; i < 3; i++)
			{
				Matrix22 grad_dadphiT;
				grad_dadphiT(0, 0) = gradDA(0, i) * dphi(0);
				grad_dadphiT(0, 1) = gradDA(0, i) * dphi(1);
				grad_dadphiT(1, 0) = gradDA(1, i) * dphi(0);
				grad_dadphiT(1, 1) = gradDA(1, i) * dphi(1);

				deriv->push_back(grad_dadphiT);
			}

			for (int i = 0; i < 3; i++)
			{
				Matrix22 grad_dadphiT;
				grad_dadphiT(0, 0) = gradDphi(0, i) * da(0);
				grad_dadphiT(0, 1) = gradDphi(1, i) * da(0);
				grad_dadphiT(1, 0) = gradDphi(0, i) * da(1);
				grad_dadphiT(1, 1) = gradDphi(1, i) * da(1);

				deriv->push_back(grad_dadphiT);
			}
		}

		if (hessian)
		{
			hessian->resize(4);

			for (int i = 0; i < 4; i++)
				hessian->at(i).setZero();

			hessian->at(0).block(0, 0, 3, 3) << hessDA[0](0, 0) * dphi(0), hessDA[0](0, 1) * dphi(0), hessDA[0](0, 2) * dphi(0),
			    hessDA[0](1, 0) * dphi(0), hessDA[0](1, 1) * dphi(0), hessDA[0](1, 2) * dphi(0),
			    hessDA[0](2, 0) * dphi(0), hessDA[0](2, 1) * dphi(0), hessDA[0](2, 2) * dphi(0);
			hessian->at(0).block(0, 3, 3, 3) << gradDA(0, 0) * gradDphi(0, 0), gradDA(0, 0) * gradDphi(0, 1), gradDA(0, 0) * gradDphi(0, 2),
			    gradDA(0, 1) * gradDphi(0, 0), gradDA(0, 1) * gradDphi(0, 1), gradDA(0, 1) * gradDphi(0, 2),
			    gradDA(0, 2) * gradDphi(0, 0), gradDA(0, 2) * gradDphi(0, 1), gradDA(0, 2) * gradDphi(0, 2);
			hessian->at(0).block(3, 0, 3, 3) = hessian->at(0).block(0, 3, 3, 3).transpose();

			hessian->at(1).block(0, 0, 3, 3) << hessDA[0](0, 0) * dphi(1), hessDA[0](0, 1) * dphi(1), hessDA[0](0, 2) * dphi(1),
			    hessDA[0](1, 0) * dphi(1), hessDA[0](1, 1) * dphi(1), hessDA[0](1, 2) * dphi(1),
			    hessDA[0](2, 0) * dphi(1), hessDA[0](2, 1) * dphi(1), hessDA[0](2, 2) * dphi(1);
			hessian->at(1).block(0, 3, 3, 3) << gradDA(0, 0) * gradDphi(1, 0), gradDA(0, 0) * gradDphi(1, 1), gradDA(0, 0) * gradDphi(1, 2),
			    gradDA(0, 1) * gradDphi(1, 0), gradDA(0, 1) * gradDphi(1, 1), gradDA(0, 1) * gradDphi(1, 2),
			    gradDA(0, 2) * gradDphi(1, 0), gradDA(0, 2) * gradDphi(1, 1), gradDA(0, 2) * gradDphi(1, 2);
			hessian->at(1).block(3, 0, 3, 3) = hessian->at(1).block(0, 3, 3, 3).transpose();

			hessian->at(2).block(0, 0, 3, 3) << hessDA[1](0, 0) * dphi(0), hessDA[1](0, 1) * dphi(0), hessDA[1](0, 2) * dphi(0),
			    hessDA[1](1, 0) * dphi(0), hessDA[1](1, 1) * dphi(0), hessDA[1](1, 2) * dphi(0),
			    hessDA[1](2, 0) * dphi(0), hessDA[1](2, 1) * dphi(0), hessDA[1](2, 2) * dphi(0);
			hessian->at(2).block(0, 3, 3, 3) << gradDA(1, 0) * gradDphi(0, 0), gradDA(1, 0) * gradDphi(0, 1), gradDA(1, 0) * gradDphi(0, 2),
			    gradDA(1, 1) * gradDphi(0, 0), gradDA(1, 1) * gradDphi(0, 1), gradDA(1, 1) * gradDphi(0, 2),
			    gradDA(1, 2) * gradDphi(0, 0), gradDA(1, 2) * gradDphi(0, 1), gradDA(1, 2) * gradDphi(0, 2);
			hessian->at(2).block(3, 0, 3, 3) = hessian->at(2).block(0, 3, 3, 3).transpose();

			hessian->at(3).block(0, 0, 3, 3) << hessDA[1](0, 0) * dphi(1), hessDA[1](0, 1) * dphi(1), hessDA[1](0, 2) * dphi(1),
			    hessDA[1](1, 0) * dphi(1), hessDA[1](1, 1) * dphi(1), hessDA[1](1, 2) * dphi(1),
			    hessDA[1](2, 0) * dphi(1), hessDA[1](2, 1) * dphi(1), hessDA[1](2, 2) * dphi(1);
			hessian->at(3).block(0, 3, 3, 3) << gradDA(1, 0) * gradDphi(1, 0), gradDA(1, 0) * gradDphi(1, 1), gradDA(1, 0) * gradDphi(1, 2),
			    gradDA(1, 1) * gradDphi(1, 0), gradDA(1, 1) * gradDphi(1, 1), gradDA(1, 1) * gradDphi(1, 2),
			    gradDA(1, 2) * gradDphi(1, 0), gradDA(1, 2) * gradDphi(1, 1), gradDA(1, 2) * gradDphi(1, 2);
			hessian->at(3).block(3, 0, 3, 3) = hessian->at(3).block(0, 3, 3, 3).transpose();
		}

		return da * dphi.transpose();
	}
	Matrix22 Sheet::computeDphi1VectorPerFaceDaTensorNew(int faceId, int quadId, std::array<Matrix22, 5>* deriv, std::array<Eigen::Matrix<Real, 5, 5>, 4>* hess)
	{
		std::array<Matrix22, 5> deriv1;
		std::array<Eigen::Matrix<Real, 5, 5>, 4> hess1;

		Matrix22 dadphTensor = computeDaDphi1VectorPerFaceTensorNew(faceId, quadId, deriv ? &deriv1 : NULL, hess ? &hess1 : NULL);

		if (deriv)
		{
			for (int i = 0; i < 5; i++)
				(*deriv)[i] = deriv1[i].transpose();
		}
		if (hess)
		{
			(*hess)[0] = hess1[0];
			(*hess)[1] = hess1[2];
			(*hess)[2] = hess1[1];
			(*hess)[3] = hess1[3];
		}
		return dadphTensor.transpose();
	}
	Matrix22 Sheet::computeDphi1VectorPerFaceDaTensor(int faceId, int quadId, std::vector<Matrix22>* deriv, std::vector<Eigen::Matrix<Real, 5, 5>>* hess)
	{
		std::vector<Matrix22> deriv1;
		std::vector<Eigen::Matrix<Real, 5, 5>> hess1;

		Matrix22 dadphTensor = computeDaDphi1VectorPerFaceTensor(faceId, quadId, deriv ? &deriv1 : NULL, hess ? &hess1 : NULL);

		if (deriv)
		{
			for (int i = 0; i < deriv1.size(); i++)
				deriv->push_back(deriv1[i].transpose());
		}
		if (hess)
		{
			hess->push_back(hess1.at(0));
			hess->push_back(hess1.at(2));
			hess->push_back(hess1.at(1));
			hess->push_back(hess1.at(3));
		}
		return dadphTensor.transpose();
	}

	Matrix22 Sheet::computeDphiDaTensor(int faceId, int quadId, std::vector<Matrix22>* deriv, std::vector<Eigen::Matrix<Real, 6, 6>>* hessian)
	{
		std::vector<Matrix22> deriv1;
		std::vector<Eigen::Matrix<Real, 6, 6>> hessian1;

		Matrix22 dadphTensor = computeDaDphiTensor(faceId, quadId, deriv ? &deriv1 : NULL, hessian ? &hessian1 : NULL);

		if (deriv)
		{
			for (int i = 0; i < deriv1.size(); i++)
				deriv->push_back(deriv1[i].transpose());
		}

		if (hessian)
		{
			hessian->resize(4);
			hessian->at(0) = hessian1[0];
			hessian->at(1) = hessian1[2];
			hessian->at(2) = hessian1[1];
			hessian->at(3) = hessian1[3];
		}

		return dadphTensor.transpose();
	}

	void Sheet::requireVertexColoringInfo()
	{
		std::cout << "\nBegin VBD info init########################" << std::endl;
		// fill vdbInfo here
		// graph coloring info
		// check if the info exists on disk
		struct stat buffer;

		std::string name = sheetInfo.sheetPath + "_vertexFacesSet.txt";
#ifndef ARGUS_QUIET
		std::cout << "Checking if color info file exists at: " << name << std::endl;
#endif
		if (stat(name.c_str(), &buffer) == 0)
		{
			argus::util::openEigenData(name, &sheetInfo.vertexFacesSet);
#ifndef ARGUS_QUIET
			std::cout << "Loaded vertexFacesSet \n";
#endif
		}
		else
		{
#ifndef ARGUS_QUIET
			std::cout << "vertexFacesSet not found! Creating.\n";
#endif
			// create vertex faces set
			sheetInfo.vertexFacesSet.resize(m_VertexCount, MAX_ADJ_LEN);
			sheetInfo.vertexFacesSet = -Eigen::MatrixXi::Ones(m_VertexCount, MAX_ADJ_LEN);
			VectorN idxMapping(m_VertexCount);
			idxMapping.setZero();
			for (int i = 0; i < m_FaceCount; i++)
			{
				trimesh::triangle_t triangle = m_Triangles.at(i);
				for (int j = 0; j < 3; j++)
				{
					unsigned int idx = triangle.v[j];
					unsigned int fillIdx = idxMapping(idx);
					if (fillIdx > MAX_ADJ_LEN)
					{
						std::cerr << "Error in creating the vertex neighbouring faces set! Total faces exceed the defined limit " << MAX_ADJ_LEN << std::endl;
						exit(0);
					}
					// std::cout << "Face " << i << " vertex " << idx << " fillIdx " << fillIdx << std::endl;
					sheetInfo.vertexFacesSet(idx, fillIdx) = i;
					idxMapping(idx)++;
				}
			}
			// exit(0);
			argus::util::writeEigenData(name, sheetInfo.vertexFacesSet);
		}
		name = sheetInfo.sheetPath + "_colorinfo.txt";
#ifdef ARGUS_DEBUG
		std::cout << "Checking if color info file exists at: " << name << std::endl;
#endif
		if (stat(name.c_str(), &buffer) == 0)
		{
			argus::util::openEigenData(name, &sheetInfo.coloredSets);
#ifdef ARGUS_DEBUG
			std::cout << "Loaded colored sets\n";
#endif
		}
		else
		{
#ifdef ARGUS_DEBUG
			std::cout << "colored sets not found! Creating.\n";
#endif
			unsigned int total_colors = 0;
			unsigned int max_vertices = 0;
			// load coloring info
			std::vector<std::vector<unsigned int>> temp;
			temp = graphColoringMcs(sheetInfo.vertexFacesSet, m_FaceCount);
			total_colors = temp.size();
			std::cout << total_colors << std::endl;
			for (unsigned int i = 0; i < total_colors; i++)
			{
				unsigned int cur = temp.at(i).size();
				if (cur > max_vertices)
					max_vertices = cur;
			}
			sheetInfo.coloredSets = -Eigen::MatrixXi::Ones(total_colors, max_vertices);
			for (unsigned int i = 0; i < total_colors; i++)
			{
				unsigned int cur = temp.at(i).size();
				for (unsigned int j = 0; j < cur; j++)
				{
					sheetInfo.coloredSets(i, j) = temp.at(i).at(j);
				}
			}

			argus::util::writeEigenData(name, sheetInfo.coloredSets);
		}

		sheetInfo.ndofs = m_VertexCount;
		sheetInfo.aat = argus::VectorN::Zero(m_VertexCount);
		sheetInfo.axt = argus::VectorN::Zero(3 * m_VertexCount);
		std::cout << "Successfully initialized VBD variables.####" << std::endl;
		sheetInfo.n_adjacentFaces.resize(m_VertexCount);
		for (int vid = 0; vid < m_VertexCount; vid++)
		{
			const Eigen::VectorXi& faces = sheetInfo.vertexFacesSet.row(vid);
			unsigned int nfaces = 0;
			for (int i = 0; i < faces.rows(); i++)
			{
				if (faces(i) == -1)
					break;
				nfaces++;
			}
			sheetInfo.n_adjacentFaces(vid) = nfaces;
		}
	}

	void Sheet::requireForcesJacobian()
	{
		const unsigned int dim = m_VertexCount * 3;
#ifdef BOGUS
		if (dim != m_ForcesPositionDeriv.rows() || dim != m_ForcesPositionDeriv.cols() || dim != m_ForcesVelocityDeriv.rows() || dim != m_ForcesVelocityDeriv.cols())
		{
			m_ForcesPositionDeriv.setCols(m_VertexCount);
			m_ForcesPositionDeriv.setRows(m_VertexCount);
			m_ForcesVelocityDeriv.setCols(m_VertexCount);
			m_ForcesVelocityDeriv.setRows(m_VertexCount);

			std::set<std::pair<size_t, size_t>> triplets = m_Material.getJacobianIndices(*this);
			m_ForcesPositionDeriv.reserve(triplets.size());
			m_ForcesVelocityDeriv.reserve(triplets.size());
			for (auto&& [i, j] : triplets)
			{
				m_ForcesPositionDeriv.insertBack(i, j);
				m_ForcesVelocityDeriv.insertBack(i, j);
			}
			m_ForcesPositionDeriv.finalize();
			m_ForcesVelocityDeriv.finalize();
		}
#endif
	}

	ForcesJacobian Sheet::getForcesJacobian()
	{
		requireForcesJacobian();
#ifdef BOGUS
		m_Material.getForcesJacobian(*this, m_ForcesPositionDeriv, m_ForcesVelocityDeriv);
		return { m_ForcesPositionDeriv, m_ForcesVelocityDeriv };
#endif
		ForcesJacobian d {};
		return d;
	}

	void Sheet::addConstraint(const LineConstraint& constraint)
	{
		m_OtherConstraints.push_back({ constraint.m_Index,
		    constraint.m_NormalizedDirection * constraint.m_NormalizedDirection.transpose() });
	}

	void Sheet::addConstraint(const PlaneConstraint& constraint)
	{
		m_OtherConstraints.push_back({ constraint.m_Index, Matrix33::Identity() - constraint.m_NormalizedNormal * constraint.m_NormalizedNormal.transpose() });
	}

	void Sheet::applyConstraints(Map<VectorN> vec)
	{
		Map<Matrix3N> vecN3 = GroupMatrixN3(vec);
		for (auto&& x : m_FixedConstraints)
		{
			vecN3.col(x.m_Index) = Vector3::Zero();
		}
		for (auto&& x : m_OtherConstraints)
		{
			vecN3.col(x.m_Index) = x.m_Matrix * vecN3.col(x.m_Index);
		}
		//  clamp constraints
		int nclampverts = m_clampedVertices.size();
		for (int i = 0; i < nclampverts; i++)
		{
			vecN3.col(m_clampedVertices(i)) *= 0;
		}
	}
	void Sheet::applyDynamicConstraints(Real dt)
	{
		if (m_ActiveDofs.size() == 0)
			throw "Error, active dofs not set!";
		// each dynamic constraint stores dofs and constraint info
		for (auto&& constraint : m_DynamicConstraints)
		{
			// pass sheet dofs (active + inactive)
			VectorN X = getDofs(false);
			VectorN Xdot = getDofsDot(false);
			// apply constraints on constraint dofs
			constraint.applyConstraints(X, Xdot, constraint.dofs, dt);
			// get active dofs
			X = X(m_ActiveDofs).eval();
			Xdot = Xdot(m_ActiveDofs).eval();
			// update sheet
			updateDofs(X, Xdot);
		}
		if (m_UseWrinkleParameters)
		{
			if (simConf.contains("wrinkle_self_coll"))
			{
				Real softness = 0;
				if (simConf["wrinkle_self_coll"].contains("softness"))
					softness = simConf["wrinkle_self_coll"]["softness"].get<Real>();
				applyWrinkleSelfCollisionConstraints(softness);
				applySmallDphiConstraints();
			}
#ifndef ARGUS_QUIET
			else
			{
				// std::cout << "Warning! Wrinkle self-collisions not enabled.";
			}
#endif
		}

		// // fix the motion of clamped vertices
		if (m_ClampedAnim.empty() == false)
		{
			if (m_ClampedAnim.find(curFrame) != m_ClampedAnim.end())
			{

				// 3. Get the map of tracked vertices for this frame
				const std::unordered_map<int, Vector3>& current_frame_data = m_ClampedAnim[curFrame];
				int nclampverts = m_clampedVertices.size();
				for (int i = 0; i < nclampverts; i++)
				{
					int target_vid = m_clampedVertices(i);
					if (current_frame_data.find(target_vid) != current_frame_data.end())
					{
						Vector3 pos = current_frame_data.at(target_vid);
						// std::cout << "Vertex " << target_vid << " pos update for frame " << curFrame << " = " << pos.transpose() << std::endl;
						m_Positions.col(target_vid) = pos;
						m_Velocities.col(target_vid) *= 0;
					}
				}
			}
		}
	}
	void Sheet::applyTopologicalConstraints()
	{
		// ensure seam vertices stay together
		for (auto pair = m_SeamVertices.begin(); pair != m_SeamVertices.end(); pair++)
		{
			unsigned int v1 = pair->first;
			unsigned int v2 = pair->second;
			Vector3 newPos = (m_Positions.col(v1) + m_Positions.col(v2)) / 2.0;
			m_Positions.col(v1) = m_Positions.col(v2) = newPos;
			Real newAmp = (m_Amplitudes(v1)); // + m_Amplitudes(v2)) / 2.0;
			m_Amplitudes(v1) = m_Amplitudes(v2) = newAmp;
		}
		for (auto it = m_clampedVertices.begin(); it != m_clampedVertices.end(); it++)
		{
			m_Velocities.col(*it) *= 0;
			m_Amplitudes(*it) = 0;
		}
	}

	bool Sheet::testTopologicalConstraints()
	{
		// ensure seam vertices stay together
		for (auto pair = m_SeamVertices.begin(); pair != m_SeamVertices.end(); pair++)
		{
			unsigned int v1 = pair->first;
			unsigned int v2 = pair->second;
			Vector3 oldpos1 = m_Positions.col(v1);
			Vector3 oldpos2 = m_Positions.col(v2);
			Real oldamp1 = m_Amplitudes(v1);
			Real oldamp2 = m_Amplitudes(v2);
			if ((oldpos1 - oldpos2).norm() > 0 || std::abs(oldamp1 - oldamp2) > 0)
			{
				std::cout << "Pair " << v1 << " and " << v2 << std::endl;
				return false;
			}
		}
		return true;
	}

	trimesh::index_t Sheet::materialPoint(const Vector3& point) const
	{
		// TODO: improve this brute force based approach
		trimesh::index_t closestIndex = -1;
		Real closestDistanceSq = std::numeric_limits<Real>::max();

		for (trimesh::index_t i = 0; i < m_Inertia.size(); i++)
		{
			Real distanceSq = (point - m_Positions.col(i)).norm();
			if (distanceSq < closestDistanceSq)
			{
				closestIndex = i;
				closestDistanceSq = distanceSq;
			}
		}
		return closestIndex;
	}

	// todo test parallel performance
	void Sheet::optimizeDphis1VectorPerFace(Real& energy, int max_iter, argus::Real tolerance)
	{
		int nfaces = m_Triangles.size();
		VectorN energies;
		energies.resize(nfaces);
		energies.setZero();
#ifdef USE_OMP
#pragma omp parallel for
#endif
		for (int face = 0; face < nfaces; face++)
		{
			// TODO: Add box solver here
			LBFGSpp::LBFGSParam<argus::Real> param;
			param.epsilon = tolerance;
			param.max_iterations = max_iter;
			LBFGSpp::LBFGSSolver solver(param);
			auto lbfgs_wrapper = [=](VectorN& x, VectorN& grad)
			{
				Vector2 init = m_dPhis1VectorPerFace.col(face);
				// tested that energies are independent
				m_dPhis1VectorPerFace.col(face) = x;
				Eigen::Vector<Real, 23> tempGrad;
				Real e = m_Material.getEnergyPerFace(*this, &tempGrad, face);
				grad.resize(2);
				grad = tempGrad(Eigen::seq(3, 4));
				m_dPhis1VectorPerFace.col(face) = init;
				return e;
			};
			// minimize
			VectorN X;
			X.resize(2);
			X = m_dPhis1VectorPerFace.col(face);
			Real energy_ = 0;
			try
			{
				solver.minimize(lbfgs_wrapper, X, energy_);
				energies(face) = energy_;
				m_dPhis1VectorPerFace.col(face) = X;
			}
			catch (const std::runtime_error& err)
			{
#ifndef ARGUS_QUIET
				std::cerr << "Error in optimizing wrinkle frequencies\n"
				          << err.what() << std::endl;
#endif
			}
			catch (...)
			{
#ifndef ARGUS_QUIET
				std::cerr << "Error in optimizing wrinkle frequencies\n"
				          << std::endl;
#endif
			}
		}
		energy = energies.sum();
	}

	void Sheet::clampDofs(VectorN& X)
	{
		unsigned int ndofs = X.size();
		// test if wrinkles are used or not
		if (m_UseWrinkleParameters)
		{
			if (ndofs == 4 * m_VertexCount)
			{
#ifdef USE_OMP
#pragma omp parallel for
#endif
				for (unsigned int i = 0; i < m_clampedVertices.size(); i++)
				{
					unsigned int vid = m_clampedVertices(i);
					// X(vid) = 0;
					X(Eigen::seq(m_VertexCount + 3 * vid, m_VertexCount + 3 * vid + 2)) = Vector3::Zero();
				}
			}
			else if (ndofs == m_VertexCount)
			{
				// clamp only amps
#ifdef USE_OMP
#pragma omp parallel for
#endif
				for (unsigned int i = 0; i < m_clampedVertices.size(); i++)
				{
					X(m_clampedVertices(i)) = 0;
				}
				// TODO: Use clampInfo
				/* for (unsigned int i = 0; i < m_clampInfo.rows(); i++)
				{
				    X(m_clampInfo(i, 0)) = 0;
				} */
			}
			else if (ndofs == 3 * m_VertexCount)
			{
				// clamp only pos
#ifdef USE_OMP
#pragma omp parallel for
#endif
				for (unsigned int i = 0; i < m_clampedVertices.size(); i++)
				{
					unsigned int vid = m_clampedVertices(i);
					X(Eigen::seq(3 * vid, 3 * vid + 2)) = Vector3::Zero();
				}
				// TODO: use clampInfo
				/* for (unsigned int i = 0; i < m_clampInfo.rows(); i++)
				{
				    Eigen::VectorXi clampV = m_clampInfo.row(i);
				    unsigned int cvid = clampV(0);
				    for (unsigned int j = 0; j < 3; j++)
				    {
				        if (clampV(j + 1))
				            X(3 * i + j) = 0;
				    }
				} */
			}
			else if (ndofs == 4 * m_VertexCount + 2 * m_FaceCount)
			{
				// clamp pos and amps
#ifdef USE_OMP
#pragma omp parallel for
#endif
				for (unsigned int i = 0; i < m_clampedVertices.size(); i++)
				{
					unsigned int vid = m_clampedVertices(i);
					X(vid) = 0;
					X(Eigen::seq(m_VertexCount + 2 * m_FaceCount + 3 * vid, m_VertexCount + 2 * m_FaceCount + 3 * vid + 2)) = Vector3::Zero();
				}
				/* for (unsigned int i = 0; i < m_clampInfo.rows(); i++)
				{
				    Eigen::VectorXi clampV = m_clampInfo.row(i);
				    unsigned int cvid = clampV(0);
				    X(cvid) = 0;
				    for (unsigned int j = 0; j < 3; j++)
				    {
				        if (clampV(j + 1))
				            X(m_VertexCount + 2 * m_FaceCount + 3 * i + j) = 0;
				    }
				} */
			}
			else
				throw "Error passed wrong vector for clamping dofs!";
		}
	}
	void Sheet::locatePotentialPureTensionFaces(std::set<int>& potentialPureTensionFaces)
	{
		// TODO remove, not needed already present in sheet class
		// Where?
		VectorN smallestEvals(m_FaceCount);

		std::vector<bool> isPureFaces(m_FaceCount, false);
		std::vector<bool> isPureVerts(m_VertexCount, false);

		for (int i = 0; i < m_FaceCount; i++)
		{
			Matrix22 Ibar = getRestShape().getIbars()[i];
			Matrix22 I;
			Matrix33 facePositions;
			int vidxi = m_Triangles[i].v[0];
			int vidxj = m_Triangles[i].v[1];
			int vidxk = m_Triangles[i].v[2];
			facePositions.col(0) = m_Positions.col(vidxi);
			facePositions.col(1) = m_Positions.col(vidxj);
			facePositions.col(2) = m_Positions.col(vidxk);
			I = firstFundamentalForm(facePositions, false, false); // correct till here
			Matrix22 diff = I - Ibar;
			Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::Matrix2d> solver(diff, Ibar);

			smallestEvals(i) = std::min(solver.eigenvalues()[0], solver.eigenvalues()[1]);
		}

		for (int i = 0; i < m_FaceCount; i++)
		{
			if (smallestEvals(i) >= 0)
			{
				isPureFaces[i] = true;
			}
		}

		for (int i = 0; i < m_FaceCount; i++)
		{
			if (isPureFaces[i])
			{
				for (int j = 0; j < 3; j++)
				{
					isPureVerts[m_Triangles.at(i).v[j]] = true;
				}
			}
		}

		//// if a face is not pure tension, but all its vertices are on the other pure faces, we think it is pure tension
		//// this is used to fill the mixed state "holes" in the pure tension region
		for (int i = 0; i < m_FaceCount; i++)
		{
			if (!isPureFaces[i])
			{
				bool isAllPureVerts = true;
				for (int j = 0; j < 3; j++)
				{
					int vid = m_Triangles.at(i).v[j];
					if (isPureVerts[vid] == false)
						isAllPureVerts = false;
				}
				isPureFaces[i] = isAllPureVerts;
			}
		}

		for (int i = 0; i < m_FaceCount; i++)
		{
			if (isPureFaces[i])
				potentialPureTensionFaces.insert(i);
		}
		// set this to empty means we never relax the integrability constriant.
		// potentialPureTensionFaces.clear();
	}

	void Sheet::estimateAmpOmegaFromStrain(Real amplitudeEstimate, VectorN& amp, Matrix2N& w, bool estimateAmps)
	{
		MatrixNN wguess(m_FaceCount, 2);
		wguess.setZero();
		std::set<int> tensionFaces;
		std::set<int> pureTensionVerts;
		Real length = m_Positions.row(0).maxCoeff() - m_Positions.row(0).minCoeff();
		Real width = m_Positions.row(1).maxCoeff() - m_Positions.row(1).minCoeff();
		Real height = m_Positions.row(2).maxCoeff() - m_Positions.row(2).minCoeff();
		Real bboxSize = std::max(length, width);
		bboxSize = std::max(bboxSize, height);
		locatePotentialPureTensionFaces(tensionFaces);
		amplitudeEstimate = 0.01 * (m_Positions.maxCoeff() - m_Positions.minCoeff());
		double newAmpEstimate = std::min(amplitudeEstimate, 0.1 * bboxSize);

		VectorN faceAmp(m_FaceCount);
		faceAmp.setZero();
		VectorN strainsPerFace(m_FaceCount);
		VectorN ampsPerFace(m_FaceCount);
		ampsPerFace.setZero();
		for (int i = 0; i < m_FaceCount; i++)
		{
			if (tensionFaces.find(i) != tensionFaces.end())
				continue;
			Matrix22 Ibar = getRestShape().getIbars()[i];
			Matrix22 I;
			Matrix33 facePositions;
			int vidxi = m_Triangles[i].v[0];
			int vidxj = m_Triangles[i].v[1];
			int vidxk = m_Triangles[i].v[2];
			facePositions.col(0) = m_Positions.col(vidxi);
			facePositions.col(1) = m_Positions.col(vidxj);
			facePositions.col(2) = m_Positions.col(vidxk);
			I = firstFundamentalForm(facePositions, false, false); // correct till here
			Matrix22 diff = I - Ibar;
			Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::Matrix2d> solver(diff, Ibar);
			Eigen::Vector2d evec = solver.eigenvectors().col(0);
			strainsPerFace(i) = std::fabs(solver.eigenvalues()[0]);
			faceAmp(i) = std::sqrt(2.0 * std::fabs(solver.eigenvalues()[0]));
			wguess.row(i) = (Ibar * evec).transpose(); // convert to one-form
			Real ai = m_Amplitudes(m_Triangles.at(i).v[0]);
			Real aj = m_Amplitudes(m_Triangles.at(i).v[1]);
			Real ak = m_Amplitudes(m_Triangles.at(i).v[2]);
			ampsPerFace(i) = (1.0 / 3.0) * (ai + aj + ak);
		}

		Real maxAmp = faceAmp.maxCoeff();
		// Real maxAmp = m_Amplitudes.maxCoeff() / 10.0;
		Real coeff = newAmpEstimate / (maxAmp + 1e-6);
		wguess = wguess / coeff;

		faceAmp = faceAmp * coeff;
		// if (false && !estimateAmps)
		// {
		// 	for (int i = 0; i < wguess.rows(); i++)
		// 	{
		// 		// get average face amplitude

		// 		if (ampsPerFace(i) >= 1e-2)
		// 		{
		// 			// wguess.row(i) = wguess.row(i).normalized();
		// 			// wguess.row(i) *= std::sqrt(2 * strainsPerFace(i)) / ((1.0 / 3.0) * (ai + aj + ak));
		// 			// if (i == 1225)
		// 			// {
		// 			// 	std::cout << coeff << " " << ampsPerFace(i) << " " << ampsPerFace.maxCoeff() << std::endl;
		// 			// 	// exit(0);
		// 			// }
		// 			Real ampratio = (ampsPerFace(i) / (ampsPerFace.maxCoeff()));
		// 			Real saferatio = std::max(0.1, ampratio);
		// 			wguess.row(i) /= saferatio;
		// 		}
		// 	}
		// }

		// std::cout << "wguess norm " << wguess.norm() << std::endl;
		// Ensures that the init dphis are aligned
		// ignoring for now
		// combField(F, abars, &faceAmp, wguess, combedW);

		w = wguess.transpose();
		if (estimateAmps)
		{
			amp.resize(m_VertexCount);
			amp.setZero();

			VectorN vdegree(m_VertexCount);
			vdegree.setZero();

			for (int i = 0; i < m_FaceCount; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					int vid = m_Triangles.at(i).v[j];
					vdegree(vid) += 1.0;
					amp(vid) += faceAmp(i);
				}
			}
			for (int i = 0; i < m_VertexCount; i++)
			{
				amp(i) /= vdegree(i);
			}

			for (auto& f : tensionFaces)
			{
				for (int j = 0; j < 3; j++)
				{
					int vid = m_Triangles.at(f).v[j];
					if (pureTensionVerts.find(vid) == pureTensionVerts.end())
						pureTensionVerts.insert(vid);
				}
			}

			for (auto& vid : pureTensionVerts)
				amp(vid) = 0;

			std::vector<std::vector<int>> vertNeis;
			igl::adjacency_list(m_faceMatrix, vertNeis);

			Eigen::SparseMatrix<Real> L(m_VertexCount, m_VertexCount);
			std::vector<Eigen::Triplet<Real>> Lcoeff;
			for (int i = 0; i < m_VertexCount; i++)
			{
				if (amp(i) > 0) // not pure tension vertices or clamped vertices
				{
					for (auto& neiV : vertNeis[i])
					{
						Lcoeff.push_back(Eigen::Triplet<double>(i, neiV, 1.0 / vertNeis[i].size()));
					}
				}
				else
					Lcoeff.push_back(Eigen::Triplet<double>(i, i, 1.0));
			}

			L.setFromTriplets(Lcoeff.begin(), Lcoeff.end());
			// apply twice by default
			amp = L * (L * amp);
			if (amp.norm() < 1e-4)
				amp = VectorN::Ones(amp.size()) * 1e-3;
		}
	}

	// wrinkle self collisions
	VectorN Sheet::getEdgeConstraints(Sheet& sheet)
	{
		std::vector<trimesh::edge_t> edges = sheet.getEdges();
		std::vector<trimesh::triangle_t> faces = sheet.getTriangles();
		std::vector<trimesh::triangle_t> VE = sheet.getVertexOppositeEdgesPerFace();
		argus::Matrix2N dphis = sheet.getDphis1VectorPerFace();
		argus::Real h = sheet.getThickness();
		// iterate through each edge to check for wrinkle self-collisions
		argus::VectorN edgeConstraints;
		Eigen::VectorXi edgeConstraintsFilled;
		edgeConstraints.resize(edges.size());
		edgeConstraintsFilled = -Eigen::VectorXi::Ones(edges.size());
		for (int i = 0; i < faces.size(); i++)
		{
			Eigen::VectorXi vids;
			vids.resize(3);
			vids[0] = faces.at(i).v[0];
			vids[1] = faces.at(i).v[1];
			vids[2] = faces.at(i).v[2];
			argus::Vector3 omegas;
			omegas(0) = dphis(0, i);
			omegas(2) = -dphis(1, i);
			omegas(1) = -omegas(0) - omegas(2);
			for (int j = 0; j < 3; j++)
			{
				int vid = vids[j];
				int eid = VE.at(i).v[j];
				if (edgeConstraintsFilled(eid) != -1)
					edgeConstraints(eid) = std::max(fabs(omegas((j + 1) % 3)), edgeConstraints(eid));
				else
				{
					edgeConstraints(eid) = fabs(omegas((j + 1) % 3));
					edgeConstraintsFilled(eid) = 1;
				}
			}
		}
		for (int i = 0; i < edges.size(); i++)
			if (edgeConstraintsFilled(i) != 1)
				throw "Error in edge constraints function";
		edgeConstraints = (edgeConstraints * h / PI).eval();
		return edgeConstraints;
	}

	std::vector<std::pair<int, Real>> Sheet::getSelfCollisionEdges(Sheet& sheet)
	{
		argus::VectorN edgeConstraints = getEdgeConstraints(sheet);
		double maxVal = edgeConstraints.maxCoeff();
		std::vector<std::pair<int, argus::Real>> problemEdges;
		std::mutex pedgemutex;
		std::vector<trimesh::edge_t> edges = sheet.getEdges();
		argus::Matrix3N positions = sheet.getPositionsN3();
		if (edges.size() != edgeConstraints.rows())
		{
			throw "Error, edge contraints not initialized!";
		}
		std::vector<std::pair<bool, Real>> ec(edges.size());

		// #ifdef USE_OMP
		// #pragma omp parallel for schedule(static)
		// #endif
		for (int i = 0; i < edges.size(); i++)
		{
			int vi = edges.at(i).v[0];
			int vj = edges.at(i).v[1];
			argus::Real edgeLen = (positions.col(vj) - positions.col(vi)).norm();
			// check for constraints
			if (edgeLen < edgeConstraints(i))
			{
				std::pair<bool, argus::Real> p(true, edgeConstraints(i));
				// pedgemutex.lock();
				problemEdges.push_back({ i, edgeConstraints(i) });
				// pedgemutex.unlock();
				ec[i] = p;
			}
			else
			{
				std::pair<bool, argus::Real> p(false, edgeConstraints(i));
				ec[i] = p;
			}
		}
		// for (int i = 0; i < edges.size(); i++)
		// {
		// 	if (ec[i].first)
		// 		problemEdges.push_back({ i, ec[i].second });
		// }
		return problemEdges;
	}
	void Sheet::getSelfCollisionEdgesNew(VectorN& collisionInfo)
	{
		// 		argus::VectorN edgeConstraints = getEdgeConstraints(sheet);
		// 		double maxVal = edgeConstraints.maxCoeff();
		// 		std::vector<std::pair<int, argus::Real>> problemEdges;
		// 		std::mutex pedgemutex;
		// 		std::vector<trimesh::edge_t> edges = sheet.getEdges();
		// 		argus::Matrix3N positions = sheet.getPositionsN3();
		// 		if (edges.size() != edgeConstraints.rows())
		// 		{
		// 			throw "Error, edge contraints not initialized!";
		// 		}
		// #ifdef USE_OMP
		// #pragma omp parallel for schedule(static)
		// #endif
		// 		for (int i = 0; i < edges.size(); i++)
		// 		{
		// 			int vi = edges.at(i).v[0];
		// 			int vj = edges.at(i).v[1];
		// 			argus::Real edgeLen = (positions.col(vj) - positions.col(vi)).norm();
		// 			// check for constraints
		// 			if (edgeLen < edgeConstraints(i))
		// 			{
		// 				std::pair<int, argus::Real> p(i, edgeConstraints(i));
		// 				pedgemutex.lock();
		// 				problemEdges.push_back(p);
		// 				pedgemutex.unlock();
		// 			}
		// 		}
		// return problemEdges;
	}
	void Sheet::applyDynamicCollisionConstraintsOnEdges(Sheet& sheet, Matrix3N& oldPositions, argus::Real dt)
	{
		std::vector<std::pair<int, argus::Real>> problemEdges = getSelfCollisionEdges(sheet);

		/* if(problemEdges.size() != 0)
		{
		    std::cout << "Problem Edges " << problemEdges.size() << std::endl;
		    exit(0);
		} */
		argus::Matrix3N& positions = sheet.getPositionsN3();
		argus::Matrix3N& velocities = sheet.getVelocitiesN3();
		std::vector<trimesh::edge_t>& edges = m_Edges;
		argus::Real h = sheet.getThickness();
		argus::VectorN massInv = sheet.getInvInertia();
		argus::Matrix3N updatedPositions = sheet.getPositionsN3();
		// continue from here
		auto getUpdatedConstraintPositions = [&sheet, problemEdges, positions, &updatedPositions, edges, h, massInv, dt](int index)
		{
			argus::Real delta = problemEdges.at(index).second;
			int eid = problemEdges.at(index).first;
			int viidx = edges.at(eid).v[0];
			int vjidx = edges.at(eid).v[1];
			argus::Vector3 xi = updatedPositions.col(viidx);
			argus::Vector3 xj = updatedPositions.col(vjidx);

			argus::Vector3 mInvi = massInv(Eigen::seq(3 * viidx, 3 * viidx + 2));
			argus::Vector3 mInvj = massInv(Eigen::seq(3 * vjidx, 3 * vjidx + 2));
			argus::Vector3 K = 2 * dt * dt * (xi - xj);
			argus::Vector3 Ki = mInvi.cwiseProduct(K);
			argus::Vector3 Kj = -mInvj.cwiseProduct(K);
			// setup quadratic coefficients
			argus::Real a = Kj.dot(Kj) - 2 * Kj.dot(Ki) + Ki.dot(Ki);
			argus::Real b = 2 * xj.dot(Kj) - 2 * xj.dot(Ki) - 2 * Kj.dot(xi) + 2 * xi.dot(Ki);
			argus::Real c = xj.dot(xj) + xi.dot(xi) - 2 * xj.dot(xi) - delta * delta;
			assert((b * b - 4 * a * c) >= 0 && "Cannot find suitable constraint projection force for wrinkle self-collisions");
			argus::Real l1 = (-b - sqrt(b * b - 4 * a * c)) / (2 * a);
			argus::Real l2 = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
			argus::Real l;

			if (l1 >= 0 && l2 >= 0)
				throw "Some error in computing constraint scale!";
			else if (l1 >= 0 && l2 <= 0)
				l = l1;
			else if (l2 >= 0 && l1 <= 0)
				l = l2;
			else
				l = 0;
			// Think when the constraint might get negative? It will be incorrect as it will push the vertices closer
			if (l < 0)
			{
#ifndef ARGUS_QUIET
				std::cerr << "Warning! Collision forces on edges pushing them inside. Results might be erroneous" << std::endl;
#endif
				l = 0;
			}
			argus::Vector3 gradCi = 2.0 * (xi - xj);
			argus::Vector3 gradCj = 2.0 * (xj - xi);
			argus::Vector3 cxi = xi + dt * dt * l * mInvi.cwiseProduct(gradCi);
			argus::Vector3 cxj = xj + dt * dt * l * mInvj.cwiseProduct(gradCj);
			// store updatedPosition
			// update positions
			updatedPositions.col(edges.at(eid).v[0]) = cxi;
			updatedPositions.col(edges.at(eid).v[1]) = cxj;
		};

		// get updated positions
		for (int i = 0; i < problemEdges.size(); i++)
		{
			getUpdatedConstraintPositions(i);
		}

		auto updateDynamicDofs = [&sheet, edges, problemEdges, updatedPositions, &positions, oldPositions, &velocities, dt](int index)
		{
			int eid = problemEdges.at(index).first;
			int vi = edges.at(eid).v[0];
			int vj = edges.at(eid).v[1];

			argus::Vector3 cxi = updatedPositions.col(vi);
			argus::Vector3 cxj = updatedPositions.col(vj);

			argus::Vector3 xi = positions.col(vi);
			argus::Vector3 xj = positions.col(vj);

			// the new velocity should be either 0 or in the direction of old velocity
			Vector3 oldxi = oldPositions.col(vi);
			Vector3 veliNew = (cxi - oldxi) / dt;
			Vector3 veliOld = velocities.col(edges.at(eid).v[0]);
			if (veliNew.norm() >= 1e-4 && veliNew.dot(veliOld) >= 0)
			{
				velocities.col(edges.at(eid).v[0]) = veliNew;
			}
			else if (veliNew.norm() >= 1e-4)
			{
				// velocities.col(edges.at(eid).v[0]) = veliNew- veliNew.dot(veliOld.normalized())*veliOld.normalized();
				velocities.col(edges.at(eid).v[0]) *= 0;
			}
			else
			{
				// no change in position
				velocities.col(edges.at(eid).v[0]) *= 0;
			}

			Vector3 oldxj = oldPositions.col(vj);
			Vector3 veljNew = (cxj - oldxj) / dt;
			Vector3 veljOld = velocities.col(edges.at(eid).v[1]);
			if (veljNew.norm() >= 1e-4 && veljNew.dot(veljOld) >= 0)
			{
				velocities.col(edges.at(eid).v[1]) = veljNew;
			}
			else if (veljNew.norm() >= 1e-4)
			{
				// velocities.col(edges.at(eid).v[1]) = veljNew - veljNew.dot(veljOld.normalized())*veljOld.normalized();
				velocities.col(edges.at(eid).v[1]) *= 0;
			}
			else
			{
				// all cases are covered
				velocities.col(edges.at(eid).v[1]) *= 0;
			}

			// fallback code if the new idea does not work
			// update velocities
			/* if((cxi-xi).norm() >= 1e-4)
			{
			    Vector3 oldx = oldPositions.col(vi);
			    velocities.col(edges.at(eid).v[0]) = (cxi - oldx)/dt;
			}

			if((cxj-xj).norm() >= 1e-4)
			{
			    Vector3 oldx = oldPositions.col(vj);
			    velocities.col(edges.at(eid).v[1]) = (cxj - oldx)/dt;
			} */
			// update positions
			positions.col(vi) = cxi;
			positions.col(vj) = cxj;
		};
		for (int i = 0; i < problemEdges.size(); i++)
		{
			updateDynamicDofs(i);
		}
	}
	// test if the triangle is squished in a way that edges do not capture inversive forces
	Real Sheet::getAlpha(Vector3 vi, Vector3 vj, Vector3 vk)
	{
		Real vivj = vi.dot(vj);
		Real vivk = vi.dot(vk);
		Real vjvj = vj.dot(vj);
		Real vjvk = vj.dot(vk);
		Real vkvk = vk.dot(vk);
		Real num = vjvk + vivj - vjvj - vivk;
		Real den = vjvk - vkvk - vjvj + vjvk;
		if (den == 0)
			throw "Error while checking for wrinkle self-collisions! The point is not on the line";
		return num / den;
	}

	std::vector<std::pair<Vector2, Eigen::Vector<int, 3>>> Sheet::getConstrainedFaceInfo(Sheet& sheet)
	{
		std::vector<std::pair<Vector2, Eigen::Vector<int, 3>>> problemInfo;
		const std::vector<trimesh::triangle_t>& faces = sheet.getTriangles();
		const Matrix3N& positions = sheet.getPositionsN3();
		const Matrix2N& dphis = sheet.getDphis1VectorPerFace();
		Real thickness = sheet.getThickness();
		int nfaces = sheet.getTriangles().size();
		int nverts = sheet.getPositionsN3().cols();
		// TODO: Parallelize
		//  std::ofstream op;
		//  op.open("debug.txt", std::ios::app);
		for (int i = 0; i < nfaces; i++)
		{
			/* op << std::endl;
			op << "Face " << i << std::endl;
			op << "dphis 0 " << dphis.col(i)(0) << " dphis 1 " << dphis.col(i)(1) << std::endl;
			op << "Vertex ids " << faces.at(i).v[0] << ", " << faces.at(i).v[1] << ", " << faces.at(i).v[2] << std::endl;
			op << std::endl; */
			Real omega01 = dphis.col(i)(0);
			Real omega20 = -dphis.col(i)(1);
			Real omega12 = -omega01 - omega20;
			for (int j = 0; j < 3; j++)
			{
				int idx0 = j;
				int idx1 = (j + 1) % 3;
				int idx2 = (j + 2) % 3;
				int vidxi = faces.at(i).v[idx0];
				int vidxj = faces.at(i).v[idx1];
				int vidxk = faces.at(i).v[idx2];
				// get positions
				Vector3 vi = positions.col(vidxi);
				Vector3 vj = positions.col(vidxj);
				Vector3 vk = positions.col(vidxk);
				Real alpha = getAlpha(vi, vj, vk);
				// get omegas
				Real omega;
				if (idx0 == 0)
				{
					omega = -omega01 - alpha * omega12;
				}

				if (idx0 == 1)
				{
					omega = -omega12 - alpha * omega20;
				}
				if (idx0 == 2)
				{
					omega = -omega20 - alpha * omega01;
				}
				Vector3 vijk = (1 - alpha) * vj + alpha * vk;
				/* op << "Inside Face " << i << "\nChosen vid " << vidxi << " j " << j << std::endl;
				op << "Orientation " << vidxi << ", " << vidxj << ", " << vidxk << std::endl;
				op << "x " << vi.transpose() << std::endl;
				op << "Alt Point " << vijk.transpose() << std::endl;
				op << "omega " << omega << std::endl;
				op << "thickness " << thickness << std::endl;
				op << "constraint " << fabs(omega * thickness/PI) << std::endl;
				op << "Altitude norm " << (vi-vijk).norm() << std::endl;
				op << "Alpha " << alpha << std::endl;  */
				// avoiding cases with zero area or poor deformations
				if (alpha > 0 && alpha < 1)
				{
					Vector3 altEdge = vijk - vi;
					Real constraint = fabs(omega * thickness / PI);
					if (altEdge.norm() < constraint)
					{
						std::pair<Vector2, Eigen::Vector<int, 3>> info;
						info.first = Vector2({ alpha, constraint });
						info.second = Eigen::Vector<int, 3>({ vidxi, vidxj, vidxk });
						problemInfo.push_back(info);
					}
				}
				// op << std::endl;
			}
		}
		// op.close();
		return problemInfo;
	}

	void Sheet::applyDynamicSelfCollisionConstraintsOnFaces(Sheet& sheet, Matrix3N& oldPositions, argus::Real dt)
	{
		int nfaces = sheet.getTriangles().size();
		int nverts = sheet.getPositionsN3().cols();
		std::vector<std::pair<Vector2, Eigen::Vector<int, 3>>> problemInfo = getConstrainedFaceInfo(sheet);

		argus::VectorN massInv = sheet.getInvInertia();
		Matrix3N updatedPositions = sheet.getPositionsN3();
		auto updatePositions = [&sheet, &massInv, dt](Matrix3N& updatedPositions, std::pair<Vector2, Eigen::Vector<int, 3>> info)
		{
			int vidxi, vidxj, vidxk;
			vidxi = info.second(0);
			vidxj = info.second(1);
			vidxk = info.second(2);

			Vector3 xi = updatedPositions.col(vidxi);
			Vector3 xj = updatedPositions.col(vidxj);
			Vector3 xk = updatedPositions.col(vidxk);

			argus::Vector3 mInvi = massInv(Eigen::seq(3 * vidxi, 3 * vidxi + 2));
			argus::Vector3 mInvj = massInv(Eigen::seq(3 * vidxj, 3 * vidxj + 2));
			argus::Vector3 mInvk = massInv(Eigen::seq(3 * vidxk, 3 * vidxk + 2));

			Real alpha = info.first(0);
			Real delta = info.first(1);
			Vector3 x_tilda = (1 - alpha) * xj + alpha * xk;

			Vector3 gradCi = (xi - x_tilda) / ((xi - x_tilda).norm());
			Vector3 Fci = dt * dt * mInvi.cwiseProduct(gradCi);
			Vector3 Fcj = dt * dt * mInvj.cwiseProduct(gradCi) * (alpha - 1);
			Vector3 Fck = dt * dt * mInvk.cwiseProduct(gradCi) * (-alpha);
			Vector3 t1 = xi - (1 - alpha) * xj - alpha * xk;
			Vector3 t2 = Fci - (1 - alpha) * Fcj - alpha * Fck;
			Real a = t2.dot(t2);
			Real b = 2 * t1.dot(t2);
			Real c = t1.dot(t1) - delta * delta;
			Real det = b * b - 4 * a * c;
			if (det < 0)
				throw "Error in finding correct forces for resolving wrinkle-self collisions\n";
			Real l1 = (-b + std::pow(det, 0.5)) / (2 * a);
			Real l2 = (-b - std::pow(det, 0.5)) / (2 * a);
			Real l = 0;
			if (l1 >= 0 && l2 >= 0)
				l = std::min(l1, l2);
			else if (l1 >= 0 && l2 < 0)
				l = l1;
			else if (l1 < 0 && l2 >= 0)
				l = l2;
			if (l < 0)
				throw "Scaling force cannot be negative while resolving face self-collisions! It is asking me to invert the triangle.\n";
			if ((t1 + l * t2).norm() - delta > 1e-4)
				throw "The solve is not satisfying the constraints while resolving wrinkle-self collisions in faces!\n";
			updatedPositions.col(vidxi) = xi + l * Fci;
			updatedPositions.col(vidxj) = xj + l * Fcj;
			updatedPositions.col(vidxk) = xk + l * Fck;
		};

		for (int i = 0; i < problemInfo.size(); i++)
		{
			try
			{
				updatePositions(updatedPositions, problemInfo.at(i));
			}
			catch (const char* error)
			{
				// std::cerr << "Error while resolving wrinkle self-collisions! " << error << std::endl;
				// throw "Error occurred while resolving wrinkle-self collisions for faces!";
			}
		}

		auto updateDynamicDofs = [&sheet, oldPositions, dt](Matrix3N& updatedPositions, std::pair<Vector2, Eigen::Vector<int, 3>> info)
		{
			Matrix3N& positions = sheet.getPositionsN3();
			Matrix3N& velocities = sheet.getVelocitiesN3();
			int vidxi, vidxj, vidxk;
			vidxi = info.second(0);
			vidxj = info.second(1);
			vidxk = info.second(2);

			Vector3 cxi = updatedPositions.col(vidxi);
			Vector3 cxj = updatedPositions.col(vidxj);
			Vector3 cxk = updatedPositions.col(vidxk);

			// the new velocity should be either 0 or in the direction of old velocity
			Vector3 oldxi = oldPositions.col(vidxi);
			Vector3 veliNew = (cxi - oldxi) / dt;
			Vector3 veliOld = velocities.col(vidxi);
			if (veliNew.norm() >= 1e-4 && veliNew.dot(veliOld) >= 0)
			{
				velocities.col(vidxi) = veliNew;
			}
			else if (veliNew.norm() >= 1e-4)
			{
				// velocities.col(edges.at(eid).v[0]) = veliNew- veliNew.dot(veliOld.normalized())*veliOld.normalized();
				velocities.col(vidxi) *= 0;
			}
			else
			{
				// no change in position
				velocities.col(vidxi) *= 0;
			}

			Vector3 oldxj = oldPositions.col(vidxj);
			Vector3 veljNew = (cxj - oldxj) / dt;
			Vector3 veljOld = velocities.col(vidxj);
			if (veljNew.norm() >= 1e-4 && veljNew.dot(veljOld) >= 0)
			{
				velocities.col(vidxj) = veljNew;
			}
			else if (veljNew.norm() >= 1e-4)
			{
				// velocities.col(edges.at(eid).v[1]) = veljNew - veljNew.dot(veljOld.normalized())*veljOld.normalized();
				velocities.col(vidxj) *= 0;
			}
			else
			{
				// all cases are covered
				velocities.col(vidxj) *= 0;
			}

			Vector3 oldxk = oldPositions.col(vidxk);
			Vector3 velkNew = (cxk - oldxk) / dt;
			Vector3 velkOld = velocities.col(vidxk);
			if (velkNew.norm() >= 1e-4 && velkNew.dot(velkOld) >= 0)
			{
				velocities.col(vidxk) = velkNew;
			}
			else if (velkNew.norm() >= 1e-4)
			{
				// velocities.col(edges.at(eid).v[1]) = veljNew - veljNew.dot(veljOld.normalized())*veljOld.normalized();
				velocities.col(vidxk) *= 0;
			}
			else
			{
				// all cases are covered
				velocities.col(vidxk) *= 0;
			}

			// fallback code if the new idea does not work
			// update velocities
			/* if((cxi-xi).norm() >= 1e-4)
			{
			    Vector3 oldx = oldPositions.col(vi);
			    velocities.col(edges.at(eid).v[0]) = (cxi - oldx)/dt;
			}

			if((cxj-xj).norm() >= 1e-4)
			{
			    Vector3 oldx = oldPositions.col(vj);
			    velocities.col(edges.at(eid).v[1]) = (cxj - oldx)/dt;
			} */
			// std::cout << "Problem face " << " alpha " << alpha << " constraint " << constraint << std::endl;
			//
			/* std::cout << "OldPositions " << oldxi.transpose() << std::endl;
			std::cout << "OldPositions " << oldxj.transpose() << std::endl;
			std::cout << "OldPositions " << oldxk.transpose() << std::endl;
			std::cout << "Positions " << cxi.transpose() << std::endl;
			std::cout << "Positions " << cxj.transpose() << std::endl;
			std::cout << "Positions " << cxk.transpose() << std::endl;
		exit(0); */
			// update positions
			positions.col(vidxi) = cxi;
			positions.col(vidxj) = cxj;
			positions.col(vidxk) = cxk;
		};

		for (int i = 0; i < problemInfo.size(); i++)
		{
			updateDynamicDofs(updatedPositions, problemInfo.at(i));
		}
	}

	void Sheet::resolveWrinkleSelfCollisions(VectorN& oldPositions, Real dt)
	{
		/*
		 * Prevent degenerate triangles by limiting wrinkle frequncies
		 */
#ifdef ARGUS_CHECKPOINT
		if (simConf.contains("Wrinkle self-Collision iterations") == false)
			simConf["Wrinkle self-Collision iterations"] = 100;
#endif
		// resolve wrinkle self-collisions
		Matrix3N oldPos = oldPositions.reshaped(3, m_VertexCount);
		// TODO: Remove collision resolution that modifies sheet positions directly and return the vector instead
		for (int i = 0; i < 5; i++)
		{
			applyDynamicCollisionConstraintsOnEdges(*this, oldPos, dt);
			applyDynamicSelfCollisionConstraintsOnFaces(*this, oldPos, dt);
		}
	}
	// Computes StVK energy
	void Sheet::computeGreenStrainPerFace()
	{
		int nfaces = m_FaceCount;
		for (int i = 0; i < nfaces; i++)
		{
			int faceId = i;
			Matrix33 facePositions;
			// facePositions.col(0) = m_Positions.block<3, 1>(0, m_Triangles[faceId].v[0]);
			// facePositions.col(1) = m_Positions.block<3, 1>(0, m_Triangles[faceId].v[1]);
			// facePositions.col(2) = m_Positions.block<3, 1>(0, m_Triangles[faceId].v[2]);
			facePositions.col(0) = m_Positions.col(m_Triangles[faceId].v[0]);
			facePositions.col(1) = m_Positions.col(m_Triangles[faceId].v[1]);
			facePositions.col(2) = m_Positions.col(m_Triangles[faceId].v[2]);
			Matrix22 I = firstFundamentalForm(facePositions, false, false, nullptr, nullptr);
			Matrix22 Ibar = m_RestShape.getIbars()[faceId];
			Real lameAlpha = getLameAlpha();
			Real lameBeta = getLameBeta();
			Matrix22 M0 = Ibar.inverse() * (I - Ibar);
			Real energy = SVNormSquared(lameAlpha, lameBeta, M0);
			m_strainValuesPerFace(faceId) = energy;
			/* if(faceId == 56)
			{
			    std::cout << I << std::endl;
			    std::cout << Ibar << std::endl;
			    std::cout << "Strain " << energy << std::endl;
			} */
		}
	}

	Real Sheet::getEhatEnergy(unsigned int vid, unsigned int fid, Real a, Vector2 dphi)
	{
		auto getAdjacentDphis = [this](int vid) -> VectorN
		{
			VectorN dphis;
			const Eigen::VectorXi& faces = sheetInfo.vertexFacesSet.row(vid);
			std::vector<int> adjfaces;
			for (int i = 0; i < faces.size(); i++)
			{
				if (faces(i) != -1)
					adjfaces.push_back(faces(i));
			}
			int nfaces = adjfaces.size();
			dphis.resize(2 * nfaces);
			for (int i = 0; i < nfaces; i++)
			{
				int fid = adjfaces.at(i);
				dphis(Eigen::seq(2 * i, 2 * i + 1)) = m_dPhis1VectorPerFace.col(fid);
			}
			return dphis;
		};
		auto getNewDphis = [this, getAdjacentDphis](int vid, Real anew)
		{
			VectorN oldDphis = getAdjacentDphis(vid);
			Real aold = m_Amplitudes(vid);
			std::vector<int> fids;
			VectorN deriv;
			MatrixNN hess;
			getVertexEnergy1(vid, fids, NULL, NULL);
			int nadjFaces = sheetInfo.n_adjacentFaces(vid);
			MatrixNN Eff;
			VectorN Ef, Eaf, dphi, dphinew, graddphinew;
			Ef = deriv(Eigen::seq(1, 2 * nadjFaces));
			Eaf = hess.col(0)(Eigen::seq(1, Eigen::last));
			Eff = hess.block(1, 1, 2 * nadjFaces, 2 * nadjFaces);
			dphi = oldDphis;
			// compatiblity check
			int Effrows, Effcols, Efsize, Eafsize, dphisize;
			dphisize = dphi.size();
			Effrows = Eff.rows();
			Effcols = Eff.cols();
			Efsize = Ef.size();
			Eafsize = Eaf.size();
			if (Effcols != Efsize)
				throw "Eff incompatible with Ef!";
			if (Effcols != Eafsize)
				throw "Eff incompatible with Eaf!";
			if (Effrows != dphisize)
				throw "Eff incompatible with dphi!";
			dphinew = dphi - Eff.inverse() * (Ef + Eaf * (anew - aold));
			return dphinew;
		};

		Real abak = m_Amplitudes(vid);
		Vector2 dphisbak = m_dPhis1VectorPerFace.col(fid);

		m_dPhis1VectorPerFace.col(fid) = dphi;
		m_Amplitudes(vid) = a;
		std::vector<int> fids;
		Real energy = getVertexEnergy1(vid, fids, NULL, NULL);
		m_Amplitudes(vid) = abak;
		m_dPhis1VectorPerFace.col(fid) = dphisbak;
		return energy;
	}
	void Sheet::writeEnergyContour(unsigned int vid, unsigned int fid, unsigned int subframe, MatrixNN stateMat, std::vector<std::pair<VectorN, Matrix2N>> dvals)
	{
		int niters = dvals.size();
		Real baka = m_Amplitudes(vid);
		Real bakh = m_dPhis1VectorPerFace.col(fid)(0);
		Real bakh1 = m_dPhis1VectorPerFace.col(fid)(1);
		Real a0 = stateMat(0, 0);
		Real a0_1 = stateMat(1, 0);
		Real h00 = stateMat(2, 0);
		Real h01 = stateMat(3, 0);
		Real a1 = stateMat(0, 1);
		Real h10 = stateMat(2, 1);
		Real h11 = stateMat(3, 1);

		Real mina = std::min(a0, a0_1);
		mina = std::min(mina, a1);
		Real maxa = std::max(a0, a0_1);
		maxa = std::max(maxa, a1);
		Real minh0 = std::min(h00, h10); // - 0.5 * std::fabs(std::min(h0, h1));
		Real maxh0 = std::max(h00, h10); // + 0.5 * std::fabs(std::max(h0, h1));

		Real minh1 = std::min(h10, h11); // - 0.5 * std::fabs(std::min(h0, h1));
		Real maxh1 = std::max(h10, h11); // + 0.5 * std::fabs(std::max(h0, h1));

		VectorN avals, hvals0, hvals1;
		avals.resize(niters);
		hvals0.resize(niters);
		hvals1.resize(niters);
		for (int i = 0; i < niters; i++)
		{
			avals(i) = dvals.at(i).first(vid);
			hvals0(i) = dvals.at(i).second.col(fid)(0);
			hvals1(i) = dvals.at(i).second.col(fid)(1);
		}
		std::string writeFolder = argus::project_root + "checkpoints/" + checkpointDirName + "/energycontour/";
		boost::filesystem::create_directories(writeFolder);
		std::string savepath;
		if (niters > 0)
		{
			savepath = writeFolder + "avals.csv." + std::to_string(curFrame) + "_" + std::to_string(subframe);
			util::writeEigenData(savepath, avals);

			savepath = writeFolder + "hvals0.csv." + std::to_string(curFrame) + "_" + std::to_string(subframe);
			util::writeEigenData(savepath, hvals0);

			savepath = writeFolder + "hvals1.csv." + std::to_string(curFrame) + "_" + std::to_string(subframe);
			util::writeEigenData(savepath, hvals1);

			Real newmina = avals.minCoeff();
			Real newmaxa = avals.maxCoeff();
			mina = std::min(mina, newmina);
			maxa = std::max(maxa, newmaxa);
			Real newminh0 = hvals0.minCoeff();
			Real newmaxh0 = hvals0.maxCoeff();
			minh0 = std::min(minh0, newminh0);
			maxh0 = std::max(maxh0, newmaxh0);
			Real newminh1 = hvals1.minCoeff();
			Real newmaxh1 = hvals1.maxCoeff();
			minh1 = std::min(minh1, newminh1);
			maxh1 = std::max(maxh1, newmaxh1);
		}

		// setup
		mina = mina - 1 * std::fabs(mina);
		maxa = maxa + 1 * std::fabs(maxa);
		minh0 = minh0 - 1 * std::fabs(minh0);
		maxh0 = maxh0 + 1 * std::fabs(maxh0);
		minh1 = minh1 - 1 * std::fabs(minh1);
		maxh1 = maxh1 + 1 * std::fabs(maxh1);

		Real samples = 100;
		MatrixNN agrid, hgrid0, hgrid1;
		Real astep, hstep0, hstep1;
		astep = (maxa - mina) / samples;
		hstep0 = (maxh0 - minh0) / samples;
		hstep1 = (maxh1 - minh1) / samples;
		int points = samples + 1;
		agrid.resize(points, points);
		hgrid0.resize(points, points);
		hgrid1.resize(points, points);
		// std::cout << " a0 " << a0 << " a1 " << a1 << " min a " << mina << " maxa " << maxa << std::endl;
		// std::cout << " h0 " << h0 << " h1 " << h1 << " min h " << minh << " maxh " << maxh << std::endl;
		MatrixNN energyGrid0, energyGrid1, ehatEnergyGrid0, ehatEnergyGrid1;
		energyGrid0.resize(points, points);
		energyGrid1.resize(points, points);
		ehatEnergyGrid0.resize(points, points);
		ehatEnergyGrid1.resize(points, points);
		for (int i = 0; i < points; i++)
		{
			Real agridval = mina + astep * Real(i);
			m_Amplitudes(vid) = agridval;
			for (int j = 0; j < points; j++)
			{
				agrid(i, j) = agridval;
				hgrid0(i, j) = minh0 + hstep0 * Real(j);
				m_dPhis1VectorPerFace.col(fid)(0) = hgrid0(i, j);
				energyGrid0(i, j) = getVertexEnergy1(vid);
				ehatEnergyGrid0(i, j) = getEhatEnergy(vid, fid, agridval, m_dPhis1VectorPerFace.col(fid));
			}
			// reset h
			m_dPhis1VectorPerFace.col(fid)(0) = bakh;
			for (int j = 0; j < points; j++)
			{
				hgrid1(i, j) = minh1 + hstep1 * Real(j);
				m_dPhis1VectorPerFace.col(fid)(1) = hgrid1(i, j);
				energyGrid1(i, j) = getVertexEnergy1(vid);
				ehatEnergyGrid1(i, j) = getEhatEnergy(vid, fid, agridval, m_dPhis1VectorPerFace.col(fid));
			}
			// reset h
			m_dPhis1VectorPerFace.col(fid)(1) = bakh1;
		}
		energyGrid1 = (energyGrid1.transpose()).eval();
		energyGrid0 = (energyGrid0.transpose()).eval();
		// write energyGrid
		savepath = writeFolder + "energies0.csv." + std::to_string(curFrame) + "_" + std::to_string(subframe);
		util::writeEigenData(savepath, energyGrid0);
		savepath = writeFolder + "energies1.csv." + std::to_string(curFrame) + "_" + std::to_string(subframe);
		util::writeEigenData(savepath, energyGrid1);

		savepath = writeFolder + "ehats0.csv." + std::to_string(curFrame) + "_" + std::to_string(subframe);
		util::writeEigenData(savepath, ehatEnergyGrid0);
		savepath = writeFolder + "ehats1.csv." + std::to_string(curFrame) + "_" + std::to_string(subframe);
		util::writeEigenData(savepath, ehatEnergyGrid1);
		// write a0,h0 and a1, h1
		savepath = writeFolder + "state.csv." + std::to_string(curFrame) + "_" + std::to_string(subframe);
		util::writeEigenData(savepath, stateMat);
		agrid = (agrid.transpose()).eval();
		savepath = writeFolder + "X.csv." + std::to_string(curFrame) + "_" + std::to_string(subframe);
		util::writeEigenData(savepath, agrid);
		hgrid0 = (hgrid0.transpose()).eval();
		savepath = writeFolder + "Y0.csv." + std::to_string(curFrame) + "_" + std::to_string(subframe);
		util::writeEigenData(savepath, hgrid0);
		hgrid1 = (hgrid1.transpose()).eval();
		savepath = writeFolder + "Y1.csv." + std::to_string(curFrame) + "_" + std::to_string(subframe);
		util::writeEigenData(savepath, hgrid1);

		m_Amplitudes(vid) = a0;
		m_dPhis1VectorPerFace.col(fid)(0) = h00;
		m_dPhis1VectorPerFace.col(fid)(1) = h01;
		Real energybefore = getVertexEnergy1(vid);
		m_Amplitudes(vid) = a1;
		m_dPhis1VectorPerFace.col(fid)(0) = h10;
		m_dPhis1VectorPerFace.col(fid)(1) = h11;
		Real energyafter = getVertexEnergy1(vid);
		Vector2 energies = Vector2({ energybefore, energyafter });
		savepath = writeFolder + "energies_ba.csv." + std::to_string(curFrame) + "_" + std::to_string(subframe);
		util::writeEigenData(savepath, energies);
		m_Amplitudes(vid) = baka;
		m_dPhis1VectorPerFace.col(fid)(0) = bakh;
		m_dPhis1VectorPerFace.col(fid)(1) = bakh1;
	}
	void Sheet::writeEnergyContourBak(unsigned int vid, unsigned int fid, unsigned int dphiid, unsigned int subframe, Real a0, Real a0_1, Real a1, Real h0, Real h1, std::vector<std::pair<VectorN, Matrix2N>> dvals)
	{
		int niters = dvals.size();
		Real baka = m_Amplitudes(vid);
		Real bakh = m_dPhis1VectorPerFace.col(fid)(0);
		Real bakh1 = m_dPhis1VectorPerFace.col(fid)(1);
		Real mina = std::min(a0, a0_1);
		mina = std::min(mina, a1);
		Real maxa = std::max(a0, a0_1);
		maxa = std::max(maxa, a1);
		Real minh = std::min(h0, h1); // - 0.5 * std::fabs(std::min(h0, h1));
		Real maxh = std::max(h0, h1); // + 0.5 * std::fabs(std::max(h0, h1));

		VectorN avals, hvals;
		avals.resize(niters);
		hvals.resize(niters);
		for (int i = 0; i < niters; i++)
		{
			avals(i) = dvals.at(i).first(vid);
			hvals(i) = dvals.at(i).second.col(fid)(dphiid);
		}
		std::string writeFolder = argus::project_root + "checkpoints/" + checkpointDirName + "/energycontour/";
		boost::filesystem::create_directories(writeFolder);
		std::string savepath;
		if (niters > 0)
		{
			savepath = writeFolder + "avals.csv." + std::to_string(curFrame) + "_" + std::to_string(subframe);
			util::writeEigenData(savepath, avals);

			savepath = writeFolder + "hvals.csv." + std::to_string(curFrame) + "_" + std::to_string(subframe);
			util::writeEigenData(savepath, hvals);

			Real newmina = avals.minCoeff();
			Real newmaxa = avals.maxCoeff();
			mina = std::min(mina, newmina);
			maxa = std::max(maxa, newmaxa);
			Real newminh = hvals.minCoeff();
			Real newmaxh = hvals.maxCoeff();
			minh = std::min(minh, newminh);
			maxh = std::max(maxh, newmaxh);
		}

		// setup
		mina = mina - 1 * std::fabs(mina);
		maxa = maxa + 1 * std::fabs(maxa);
		minh = minh - 1 * std::fabs(minh);
		maxh = maxh + 1 * std::fabs(maxh);
		Real samples = 100;
		MatrixNN agrid, hgrid;
		Real astep, hstep;
		astep = (maxa - mina) / samples;
		hstep = (maxh - minh) / samples;
		int points = samples + 1;
		agrid.resize(points, points);
		hgrid.resize(points, points);
		// std::cout << " a0 " << a0 << " a1 " << a1 << " min a " << mina << " maxa " << maxa << std::endl;
		// std::cout << " h0 " << h0 << " h1 " << h1 << " min h " << minh << " maxh " << maxh << std::endl;
		MatrixNN energyGrid, energyGrid1;
		energyGrid.resize(points, points);
		energyGrid1.resize(points, points);
		for (int i = 0; i < points; i++)
		{
			Real agridval = mina + astep * Real(i);
			m_Amplitudes(vid) = agridval;
			for (int j = 0; j < points; j++)
			{
				agrid(i, j) = agridval;
				hgrid(i, j) = minh + hstep * Real(j);
				m_dPhis1VectorPerFace.col(fid)(dphiid) = hgrid(i, j);
				energyGrid(i, j) = getVertexEnergy1(vid);
			}
		}
		energyGrid = (energyGrid.transpose()).eval();
		savepath = writeFolder + "energies.csv." + std::to_string(curFrame) + "_" + std::to_string(subframe);
		// write energyGrid
		util::writeEigenData(savepath, energyGrid);
		// write a0,h0 and a1, h1
		VectorN state;
		state.resize(5);
		state(0) = a0;
		state(1) = h0;
		state(2) = a1;
		state(3) = h1;
		state(4) = a0_1;

		savepath = writeFolder + "state.csv." + std::to_string(curFrame) + "_" + std::to_string(subframe);
		util::writeEigenData(savepath, state);
		agrid = (agrid.transpose()).eval();
		savepath = writeFolder + "X.csv." + std::to_string(curFrame) + "_" + std::to_string(subframe);
		util::writeEigenData(savepath, agrid);
		hgrid = (hgrid.transpose()).eval();
		savepath = writeFolder + "Y.csv." + std::to_string(curFrame) + "_" + std::to_string(subframe);
		util::writeEigenData(savepath, hgrid);
		m_Amplitudes(vid) = a0;
		m_dPhis1VectorPerFace.col(fid)(dphiid) = h0;
		Real energybefore = getVertexEnergy1(vid);
		m_Amplitudes(vid) = a1;
		m_dPhis1VectorPerFace.col(fid)(dphiid) = h1;
		Real energyafter = getVertexEnergy1(vid);
		Vector2 energies = Vector2({ energybefore, energyafter });
		savepath = writeFolder + "energies_ba.csv." + std::to_string(curFrame) + "_" + std::to_string(subframe);
		util::writeEigenData(savepath, energies);
		m_Amplitudes(vid) = baka;
		m_dPhis1VectorPerFace.col(fid)(dphiid) = bakh;
	}
	void Sheet::cacheParameters(std::array<bool, 5> parameters)
	{
		/*
		 * Cache parameters before running energies
		 * Warning! Do not run in parallel loops
		 */
		if (parameters[0])
		{
			// cache I
			m_cachedI.first = true;
#ifdef USE_OMP
#pragma omp parallel for schedule(static)
#endif
			for (unsigned int i = 0; i < m_FaceCount; i++)
			{
				Matrix33 facePositions;
				facePositions.col(0) = m_Positions.col(m_Triangles[i].v[0]);
				facePositions.col(1) = m_Positions.col(m_Triangles[i].v[1]);
				facePositions.col(2) = m_Positions.col(m_Triangles[i].v[2]);
				m_cachedI.second[i] = firstFundamentalForm(facePositions, false, false, NULL, NULL); // correct till here
			}
		}
		else
		{
			m_cachedI.first = false;
		}
		if (parameters[1])
		{
			// cache II
			m_cachedII.first = true;

#ifdef USE_OMP
#pragma omp parallel for schedule(static)
#endif
			for (unsigned int i = 0; i < m_FaceCount; i++)
			{
				// std::shared_ptr<SecondFundamentalFormDiscretization> sff;
				// sff = std::make_shared<MidedgeAverageFormulation>();
				m_cachedII.second[i] = sffobject->secondFundamentalFormNew(m_baseMesh, m_Positions, m_edgeDofs, i, NULL, NULL);
			}
		}
		else
		{
			m_cachedII.first = false;
		}
		if (parameters[2])
		{
			// cache DaDaT
			m_cachedDaDaT.first = true;

#ifdef USE_OMP
#pragma omp parallel for schedule(static)
#endif
			for (unsigned int i = 0; i < m_FaceCount; i++)
			{
				m_cachedDaDaT.second[i] = computeDaDaTensorNew(i, 0, NULL, NULL);
			}
		}
		else
		{
			m_cachedDaDaT.first = false;
		}
		if (parameters[3])
		{
			// cache DphiDphiT
			m_cachedDphiDphiT.first = true;

#ifdef USE_OMP
#pragma omp parallel for schedule(static)
#endif
			for (unsigned int i = 0; i < m_FaceCount; i++)
			{
				m_cachedDphiDphiT.second[i] = computeDphiDphi1VectorPerFaceTensorNew(i, 0, NULL, NULL);
			}
		}
		else
		{
			m_cachedDphiDphiT.first = false;
		}
		if (parameters[4])
		{
			// cache DaDphiT
			m_cachedDaDphiT.first = true;

#ifdef USE_OMP
#pragma omp parallel for schedule(static)
#endif
			for (unsigned int i = 0; i < m_FaceCount; i++)
			{
				m_cachedDaDphiT.second[i] = computeDaDphi1VectorPerFaceTensorNew(i, 0, NULL, NULL);
			}
		}
		else
		{
			m_cachedDaDphiT.first = false;
		}
	}
	void Sheet::clearCache()
	{
		m_cachedDaDaT.first = false;
		m_cachedDaDphiT.first = false;
		m_cachedDphiDphiT.first = false;
		m_cachedI.first = false;
		m_cachedII.first = false;
	}
	// 	void Sheet::applyWrinkleSelfCollisionConstraints()
	// 	{
	// 		// Math
	// 		// For each face, we have a deformation gradient F and wrinkle frequencies stored as one forms on two edges in reference mesh (w0,w1)
	// 		// Now, we need to ensure that wrinkles over an edge do not expand outside the base mesh length along the wrinkle:
	// 		// omega * h / PI <= ||edge|| - eq 1
	// 		// We intend to derive the constraint relation to prevent this and ensure constraints are satisfied vid pbd
	// 		// Solution: Let's choose an arbitrary vector v, such that it's deformed state would be F v
	// 		// Now, since the energy model assumes a single wrinkling direction in a face, we can safely assume that
	// 		// if ||F tildav|| = sigma2 ||tildav|| (i.e., tildav is the maximum compression direction along wrinkle) - eq 2
	// 		// from eq 1 and eq 2, we can say that
	// 		// |gradphi . v| * h/PI <= ||F v||, square both sides
	// 		//   v^T F^TF v -v^T(h*h/(PI*PI)) gradphi grapdphi^T v >=0
	// 		// v^T (F^T F - (h*h/(PI*PI) gradphi * gradphi^T) v >=0
	// 		// let h*h/(PI*PI) = deltasq
	// 		// this requires us the ensure that matrix (F^T F - deltasq gradphi gradphi^T) is positive semidefinite
	// 		// we know F^T F = Y^TIY is symmetric, and we know that outer product is symmetric (here Y=X^-1, see notes from fem)
	// 		// so the matrix A = Y^TIY - deltasq * gradphi * gradphi^T is symmetric
	// 		// then we just need to ensure that the sum of eigen values of A and the product of eigen values of A are positive
	// 		// det(A) >=0 and trace(A) >=0
	// 		// First assume I = [[l,m],[m,c]]
	// 		// And Y = [[tildaa, tildab],[tildac, tildad]]
	// 		// Then, Y^TIY =[[a,b],[b,c]]
	// 		// use det trick to reduce trace(A) constraint and get constraints
	// 		// C1 := a - deltasq * w0*w0 > 0 (i.e., A(0,0)>0)
	// 		// C2 := det(A)>0

	// 		Real deltasq = (getThickness() * getThickness() / (PI * PI));
	// 		// filter out faces that do not fail constraints
	// 		std::vector<bool> c1fails(m_FaceCount);
	// 		std::vector<bool> c2fails(m_FaceCount);
	// 		// std::vector<bool> c3fails(m_FaceCount);
	// 		int niters = 5; // TODO: What is the optimal number of iterations
	// 		if (!m_RestShape.getIToFTF().first)
	// 			throw "Error, matrix not initialized to process wrinkle self-collisions!";
	// 		for (int iter = 0; iter < niters; iter++)
	// 		{
	// #ifdef USE_OMP
	// #pragma omp parallel for schedule(static)
	// #endif
	// 			for (int faceId = 0; faceId < m_FaceCount; faceId++)
	// 			{
	// 				int vi = m_Triangles[faceId].v[0];
	// 				int vj = m_Triangles[faceId].v[1];
	// 				int vk = m_Triangles[faceId].v[2];
	// 				Matrix22 Y = m_RestShape.getIToFTF().second[faceId].inverse();
	// 				Matrix22 X = m_RestShape.getIToFTF().second[faceId];
	// 				Vector3 eij = m_Positions.col(vj) - m_Positions.col(vi);
	// 				Vector3 eik = m_Positions.col(vk) - m_Positions.col(vi);
	// 				Real w0oldBasis = m_dPhis1VectorPerFace.col(faceId)(0);
	// 				Real w1oldBasis = m_dPhis1VectorPerFace.col(faceId)(1);
	// 				Real w0 = w0oldBasis / X(0, 0);
	// 				Real w1 = (w1oldBasis - w0oldBasis * (X(0, 1) / X(0, 0))) / X(1, 1);
	// 				Real l = eij.dot(eij);
	// 				Real m = eij.dot(eik);
	// 				Real n = eik.dot(eik);
	// 				Matrix22 I = Matrix22({ { l, m }, { m, n } });
	// 				Matrix22 YTIY = Y.transpose() * I * Y;
	// 				if ((YTIY - YTIY.transpose()).norm() >= 1e-4 * YTIY.norm())
	// 					throw "Error FTF is not symmetric!";
	// 				Real a = YTIY(0, 0); // tildaa * (tildaa * l + tildac * m) + tildac * (tildaa * m + tildac * n);
	// 				Real b = YTIY(0, 1); // tildab * (tildaa * l + tildac * m) + tildad * (tildaa * m + tildac * n);
	// 				Real c = YTIY(1, 1); // tildab * (tildab * l + tildad * m) + tildad * (tildab * m + tildad * n);
	// 				// first constraint
	// 				Real C1 = a - deltasq * w0 * w0;
	// 				if (C1 < 0)
	// 					c1fails[faceId] = true;
	// 				else
	// 					c1fails[faceId] = false;
	// 				Real C2 = a * c - deltasq * (a * w1 * w1 + c * w0 * w0 - 2 * b * w0 * w1) - b * b;
	// 				if (C2 < 0)
	// 					c2fails[faceId] = true;
	// 				else
	// 					c2fails[faceId] = false;
	// 				/* Real C3 = c - deltasq * w1 * w1;
	// 				if (C3 < 0)
	// 				    c3fails[faceId] = true;
	// 				else
	// 				    c3fails[faceId] = false; */
	// 			}
	// 			for (int faceId = 0; faceId < m_FaceCount; faceId++)
	// 			{

	// 				if (!(c1fails[faceId] || c2fails[faceId]))
	// 					continue;
	// 				int vi = m_Triangles[faceId].v[0];
	// 				int vj = m_Triangles[faceId].v[1];
	// 				int vk = m_Triangles[faceId].v[2];
	// 				Real minvi = 1.0 / m_Inertia(vi);
	// 				Real minvj = 1.0 / m_Inertia(vj);
	// 				Real minvk = 1.0 / m_Inertia(vk);
	// 				if (minvi == 0 && minvj == 0 && minvk == 0)
	// 					continue;
	// 				Matrix22 Y = m_RestShape.getIToFTF().second[faceId].inverse();
	// 				Vector3 eij = m_Positions.col(vj) - m_Positions.col(vi);
	// 				Vector3 eik = m_Positions.col(vk) - m_Positions.col(vi);
	// 				// Real w0 = m_dPhis1VectorPerFace.col(faceId)(0);
	// 				// Real w1 = m_dPhis1VectorPerFace.col(faceId)(1);
	// 				Real w0oldBasis = m_dPhis1VectorPerFace.col(faceId)(0);
	// 				Real w1oldBasis = m_dPhis1VectorPerFace.col(faceId)(1);
	// 				Matrix22 X = m_RestShape.getIToFTF().second[faceId];
	// 				Real w0 = w0oldBasis / X(0, 0);
	// 				Real w1 = (w1oldBasis - w0oldBasis * (X(0, 1) / X(0, 0))) / X(1, 1);

	// 				Real l = eij.dot(eij);
	// 				Real m = eij.dot(eik);
	// 				Real n = eik.dot(eik);
	// 				Matrix22 I = Matrix22({ { l, m }, { m, n } });
	// 				Matrix22 YTIY = Y.transpose() * I * Y;
	// 				if ((YTIY - YTIY.transpose()).norm() >= 1e-4 * YTIY.norm())
	// 					throw "Error FTF is not symmetric!";
	// 				Real a = YTIY(0, 0); // tildaa * (tildaa * l + tildac * m) + tildac * (tildaa * m + tildac * n);
	// 				Real b = YTIY(0, 1); // tildab * (tildaa * l + tildac * m) + tildad * (tildaa * m + tildac * n);
	// 				Real c = YTIY(1, 1); // tildab * (tildab * l + tildad * m) + tildad * (tildab * m + tildad * n);
	// 				Real tildaa = Y(0, 0);
	// 				Real tildab = Y(0, 1);
	// 				Real tildac = Y(1, 0);
	// 				Real tildad = Y(1, 1);

	// 				Vector3 dldxi = -2 * eij;
	// 				Vector3 dmdxi = -eij - eik;
	// 				Vector3 dndxi = -2 * eik;
	// 				Vector3 dldxj = 2 * eij;
	// 				Vector3 dmdxj = eik;
	// 				Vector3 dndxj = Vector3::Zero();
	// 				Vector3 dldxk = Vector3::Zero();
	// 				Vector3 dmdxk = eij;
	// 				Vector3 dndxk = 2 * eik;
	// 				if (c1fails[faceId])
	// 				{
	// 					Vector3 dadxi = tildaa * (tildaa * dldxi + tildac * dmdxi) + tildac * (tildaa * dmdxi + tildac * dndxi);
	// 					Vector3 dadxj = tildaa * (tildaa * dldxj + tildac * dmdxj) + tildac * (tildaa * dmdxj + tildac * dndxj);
	// 					// constraint violated
	// 					// first constraint
	// 					Real C1 = a - deltasq * w0 * w0;
	// 					Vector3 gradCi = dadxi;
	// 					Vector3 gradCj = dadxj;
	// 					Real s = -C1 / (minvi * gradCi.dot(gradCi) + minvj * gradCj.dot(gradCj));
	// 					// update i
	// 					if (m_InvInertia.col(vi)(0) != 0)
	// 					{
	// 						m_Positions.col(vi) += s * gradCi;
	// 					}
	// 					// update j
	// 					if (m_InvInertia.col(vj)(0) != 0)
	// 					{
	// 						m_Positions.col(vj) += s * gradCj;
	// 					}
	// 				}
	// 				if (c2fails[faceId])
	// 				{
	// 					Real C2 = a * c - deltasq * (a * w1 * w1 + c * w0 * w0 - 2 * b * w0 * w1) - b * b;
	// 					Vector3 dadxi = tildaa * (tildaa * dldxi + tildac * dmdxi) + tildac * (tildaa * dmdxi + tildac * dndxi);
	// 					Vector3 dadxj = tildaa * (tildaa * dldxj + tildac * dmdxj) + tildac * (tildaa * dmdxj + tildac * dndxj);
	// 					Vector3 dadxk = tildaa * (tildaa * dldxk + tildac * dmdxk) + tildac * (tildaa * dmdxk + tildac * dndxk);
	// 					Vector3 dbdxi = tildaa * (tildab * dldxi + tildad * dmdxi) + tildac * (tildab * dmdxi + tildad * dndxi);
	// 					Vector3 dbdxj = tildaa * (tildab * dldxj + tildad * dmdxj) + tildac * (tildab * dmdxj + tildad * dndxj);
	// 					Vector3 dbdxk = tildaa * (tildab * dldxk + tildad * dmdxk) + tildac * (tildab * dmdxk + tildad * dndxk);
	// 					Vector3 dcdxi = tildab * (tildab * dldxi + tildad * dmdxi) + tildad * (tildab * dmdxi + tildad * dndxi);
	// 					Vector3 dcdxj = tildab * (tildab * dldxj + tildad * dmdxj) + tildad * (tildab * dmdxj + tildad * dndxj);
	// 					Vector3 dcdxk = tildab * (tildab * dldxk + tildad * dmdxk) + tildad * (tildab * dmdxk + tildad * dndxk);
	// 					Real dC2dxa = c - deltasq * w1 * w1;
	// 					Real dC2dxb = 2 * (deltasq * w0 * w1 - b);
	// 					Real dC2dxc = a - deltasq * w0 * w0;
	// 					Vector3 gradCi = dC2dxa * dadxi + dC2dxb * dbdxi + dC2dxc * dcdxi;
	// 					Vector3 gradCj = dC2dxa * dadxj + dC2dxb * dbdxj + dC2dxc * dcdxj;
	// 					Vector3 gradCk = dC2dxa * dadxk + dC2dxb * dbdxk + dC2dxc * dcdxk;
	// 					Real s = -C2 / (minvi * gradCi.dot(gradCi) + minvj * gradCj.dot(gradCj) + minvk * gradCk.dot(gradCk));
	// 					// update i
	// 					if (m_InvInertia.col(vi)(0) != 0)
	// 					{
	// 						m_Positions.col(vi) += s * gradCi;
	// 					}
	// 					// update j
	// 					if (m_InvInertia.col(vj)(0) != 0)
	// 					{
	// 						m_Positions.col(vj) += s * gradCj;
	// 					}
	// 					// update k
	// 					if (m_InvInertia.col(vk)(0) != 0)
	// 					{
	// 						m_Positions.col(vk) += s * gradCk;
	// 					}
	// 				}
	// 			}
	// 		}
	// 	}
	void Sheet::applySmallDphiConstraints()
	{
#ifdef USE_OMP
#pragma omp parallel for
#endif
		for (int faceId = 0; faceId < m_FaceCount; faceId++)
		{
			argus::Vector2 dphi = m_dPhis1VectorPerFace.col(faceId);
			const argus::Matrix32& M = m_dphi_to_gradphis.at(faceId);
			Real delta = m_minimum_edge_lengths.at(faceId) * PI / 2.0;
			argus::Vector3 gradphi = M * dphi;
			if (gradphi.norm() < delta)
			{
				m_dPhis1VectorPerFace.col(faceId) *= 1e-4;
			}
		}
	}
	void Sheet::applyWrinkleSelfCollisionConstraints(Real softness)
	{
		// Math
		// For each face, we have a deformation gradient F and wrinkle frequencies stored as one forms on two edges in reference mesh (w0,w1)
		// Now, we need to ensure that wrinkles over an edge do not expand outside the base mesh length along the wrinkle:
		// omega * h / PI <= ||edge|| - eq 1
		// We intend to derive the constraint relation to prevent this and ensure constraints are satisfied vid pbd
		// Solution: Let's choose an arbitrary vector v, such that it's deformed state would be F v
		// Now, since the energy model assumes a single wrinkling direction in a face, we can safely assume that
		// if ||F tildav|| = sigma2 ||tildav|| (i.e., tildav is the maximum compression direction along wrinkle) - eq 2
		// from eq 1 and eq 2, we can say that
		// |gradphi . v| * h/PI <= ||F v||, square both sides
		//   v^T F^TF v -v^T(h*h/(PI*PI)) gradphi grapdphi^T v >=0
		// v^T (F^T F - (h*h/(PI*PI) gradphi * gradphi^T) v >=0
		// let h*h/(PI*PI) = deltasq
		// this requires us the ensure that matrix (F^T F - deltasq gradphi gradphi^T) is positive semidefinite
		// we know F^T F = I is symmetric, and we know that outer product is symmetric
		// so the matrix A = I - deltasq * gradphi * gradphi^T is symmetric
		// then we just need to ensure that the sum of eigen values of A and the product of eigen values of A are positive
		// det(A) >=0 and trace(A) >=0
		// First assume I = [[a,b],[b,c]]
		// use det trick to reduce trace(A) constraint and get constraints
		// C1 := a - deltasq * w0*w0 > 0 (i.e., A(0,0)>0)
		// C2 := det(A)>0

		Real deltasq = (getThickness() * getThickness() / (PI * PI));
		// filter out faces that do not fail constraints
		std::vector<bool> c1fails(m_FaceCount);
		std::vector<bool> c2fails(m_FaceCount);
		int niters = 5; // TODO: What is the optimal number of iterations
		for (int iter = 0; iter < niters; iter++)
		{
#ifdef USE_OMP
#pragma omp parallel for schedule(static)
#endif
			for (int faceId = 0; faceId < m_FaceCount; faceId++)
			{
				int vi = m_Triangles[faceId].v[0];
				int vj = m_Triangles[faceId].v[1];
				int vk = m_Triangles[faceId].v[2];
				Vector3 eij = m_Positions.col(vj) - m_Positions.col(vi);
				Vector3 eik = m_Positions.col(vk) - m_Positions.col(vi);
				Real w0 = m_dPhis1VectorPerFace.col(faceId)(0);
				Real w1 = m_dPhis1VectorPerFace.col(faceId)(1);
				Real a = eij.dot(eij);
				Real b = eij.dot(eik);
				Real c = eik.dot(eik);
				// first constraint
				Real C1 = a - deltasq * w0 * w0;
				if (C1 < 0)
					c1fails[faceId] = true;
				else
					c1fails[faceId] = false;
				Real C2 = a * c - deltasq * (a * w1 * w1 + c * w0 * w0 - 2 * b * w0 * w1) - b * b;
				if (C2 < 0)
					c2fails[faceId] = true;
				else
					c2fails[faceId] = false;
			}
			for (int faceId = 0; faceId < m_FaceCount; faceId++)
			{

				if (!(c1fails[faceId] || c2fails[faceId]))
					continue;
				int vi = m_Triangles[faceId].v[0];
				int vj = m_Triangles[faceId].v[1];
				int vk = m_Triangles[faceId].v[2];
				Real minvi = 1.0 / m_Inertia(vi);
				Real minvj = 1.0 / m_Inertia(vj);
				Real minvk = 1.0 / m_Inertia(vk);
				if (minvi == 0 && minvj == 0 && minvk == 0)
					continue;
				Vector3 eij = m_Positions.col(vj) - m_Positions.col(vi);
				Vector3 eik = m_Positions.col(vk) - m_Positions.col(vi);
				Real w0 = m_dPhis1VectorPerFace.col(faceId)(0);
				Real w1 = m_dPhis1VectorPerFace.col(faceId)(1);
				Real a = eij.dot(eij);
				Real b = eij.dot(eik);
				Real c = eik.dot(eik);
				if (c1fails[faceId])
				{
					// constraint violated
					// first constraint
					Real C1 = a - deltasq * w0 * w0;
					Vector3 gradCi = -2 * (eij);
					Vector3 gradCj = 2 * eij;
					Real s = -C1 / (softness + minvi * gradCi.dot(gradCi) + minvj * gradCj.dot(gradCj));
					// update i
					if (m_InvInertia.col(vi)(0) != 0)
					{
						m_Positions.col(vi) += s * gradCi;
					}
					// update j
					if (m_InvInertia.col(vj)(0) != 0)
					{
						m_Positions.col(vj) += s * gradCj;
					}
				}
				if (c2fails[faceId])
				{
					Real C2 = a * c - deltasq * (a * w1 * w1 + c * w0 * w0 - 2 * b * w0 * w1) - b * b;
					Vector3 dadxi = -2 * eij;
					Vector3 dbdxi = -eij - eik;
					Vector3 dcdxi = -2 * eik;
					Vector3 dadxj = 2 * eij;
					Vector3 dbdxj = eik;
					Vector3 dcdxj = Vector3::Zero();
					Vector3 dadxk = Vector3::Zero();
					Vector3 dbdxk = eij;
					Vector3 dcdxk = 2 * eik;
					Real dC2dxa = c - deltasq * w1 * w1;
					Real dC2dxb = 2 * (deltasq * w0 * w1 - b);
					Real dC2dxc = a - deltasq * w0 * w0;
					Vector3 gradCi = dC2dxa * dadxi + dC2dxb * dbdxi + dC2dxc * dcdxi;
					Vector3 gradCj = dC2dxa * dadxj + dC2dxb * dbdxj + dC2dxc * dcdxj;
					Vector3 gradCk = dC2dxa * dadxk + dC2dxb * dbdxk + dC2dxc * dcdxk;
					Real s = -C2 / (softness + minvi * gradCi.dot(gradCi) + minvj * gradCj.dot(gradCj) + minvk * gradCk.dot(gradCk));
					// update i
					if (m_InvInertia.col(vi)(0) != 0)
					{
						m_Positions.col(vi) += s * gradCi;
					}
					// update j
					if (m_InvInertia.col(vj)(0) != 0)
					{
						m_Positions.col(vj) += s * gradCj;
					}
					// update k
					if (m_InvInertia.col(vk)(0) != 0)
					{
						m_Positions.col(vk) += s * gradCk;
					}
				}
			}
		}
	}
}
