#ifndef SHEET_HPP
#define SHEET_HPP

#include "common.hpp"
#include <mutex>
#include <memory>
#include <igl/boundary_facets.h>
#include <igl/cotmatrix.h>
#include "material.hpp"
#include "dynamics_object.hpp"
#include "constraints.hpp"
#include "rest_shape.hpp"
#include <LBFGS.h>
#include <LBFGSB.h>
#include "../../vendor/MeshLib/MeshConnectivity.h"
#include "../vendor/MeshLib/IntrinsicGeometry.h"

namespace argus
{
	void saveCheckpoint(std::vector<std::string> data_to_save, std::string parent_directory, int frame, argus::Sheet& sheet, nlohmann::json state);
	class Sheet : public DynamicsObject
	{
		Material m_Material;
		RestShape m_RestShape;
		std::mutex m_mutex;
		argus::Matrix3N m_Positions;
		std::vector<trimesh::triangle_t> m_UVTriangles; // refer to different vertices for seams
		std::vector<trimesh::triangle_t> m_Triangles;
		std::vector<trimesh::edge_t> m_Edges;
		trimesh::trimesh_t m_HalfEdges;
		argus::Matrix2N m_UVs;
		argus::MatrixNN m_VerticesDegree;
		std::vector<FixedConstraint> m_FixedConstraints;
		std::vector<ConstraintMatrix> m_OtherConstraints;
		std::vector<DynamicConstraint> m_DynamicConstraints;
		VectorN m_Inertia;
		Matrix3N m_InvInertia;
		Matrix3N m_Velocities;
		Matrix3N m_Forces;
		VectorN m_WrinkleForces;

#ifdef BOGUS
		BlockSparseMatrix m_ForcesPositionDeriv;
		BlockSparseMatrix m_ForcesVelocityDeriv;
#endif

		trimesh::index_t m_VertexCount;
		trimesh::index_t m_EdgeCount;
		int m_FaceCount;

		// stored clamped dofs
		// TODO : remove
		std::vector<int> m_ClampedDofs;
		std::vector<int> m_ActiveDofs;

		// For Wrinkle stuff
		// is wrinkle mesh?
		bool m_UseWrinkleParameters;
		argus::Matrix3N m_WrinkleMeshPositions;
		Eigen::MatrixXi m_WrinkleMeshFaces;
		bool m_UseReferenceMesh;
		argus::Matrix3N m_ReferenceMeshPositions; // reference mesh
		std::vector<trimesh::triangle_t> m_ReferenceMeshTriangles; // reference faces: the faces at seam have different vertices than in the simulation mesh
		std::map<int, int> m_SeamVertices; // contain the mapping seam vertex in refernce space -> seam vertex in rest space
		argus::Matrix3N m_RestPositions;
		argus::VectorN m_Amplitudes;
		argus::VectorN m_ampsPerFace;
		argus::VectorN m_dPhis;
		argus::VectorN m_Phis;
		Eigen::Vector<std::complex<Real>, Eigen::Dynamic> m_zPhis;
		argus::VectorN m_PhiVelocities;
		argus::Matrix3N m_CornerPhis; // store the phi values per corner of a face
		argus::Matrix2N m_dPhis1VectorPerFace;
		argus::Matrix3N m_dPhisPerVertex; // Dphis in 3D World coordinates
		Eigen::MatrixXi m_faceMatrix;
		Eigen::SparseMatrix<int> m_adjacencyMat;
		// store faces adjacent to each vertex
		std::vector<std::vector<int>> m_vertexFaceAdjacencyList;
		// store edges adjacent to each vertex
		std::vector<std::vector<int>> m_vertexEdgeAdjacencyList;
		std::vector<trimesh::triangle_t> m_WrinkledTriangles;
		std::vector<argus::Vector3> m_faceNormals;
		std::vector<trimesh::triangle_t> m_VtexOppositeEdgesPerFace; // m_VtexOppositeEdgesPerFace[faceId].v[i] = edgeId opposite to vertex i in face faceId
		Eigen::MatrixXi m_EdgeOppositeVertices; // m_EdgeOppositeVertices(edgeId, i) = vertexId opposite to edgeId i=0 for first vertex and i=1 for the second vertex
		Eigen::SparseMatrix<argus::Real> laplacian;
		Eigen::SparseMatrix<argus::Real> m_DampingLaplacian;
		std::vector<argus::Matrix22> m_IValues;
		std::vector<QuadraturePoints> m_quadPoints;
		std::vector<argus::Real> m_vertAreas;
		Eigen::Vector<Real, Eigen::Dynamic> m_faceAreas;
		// TODO: merge with the existing constraints
		// std::vector<bool> m_freeDOFMap;
		//
		// std::vector<int> m_idxOfDOFs;
		// std::vector<int> m_freeDOFInvIdxMap;
		VectorN m_edgeDofs;
		MeshConnectivity m_baseMesh;
		std::string checkpointDirName; // a string containing the number of seconds passed since 1970 Jan 1
		std::map<std::string, int> paramsToSave;
		nlohmann::json checkpointParams;
		// Store the mesh optimization energies
		VectorN m_elasticEnergies;
		VectorN m_wrinkleEnergies;
		// store clamped vertices info
		Eigen::VectorXi m_clampedVertices; // to store the vertex indices
		Eigen::MatrixXi m_clampInfo; // stores individual clamped dimensions
		// Store per face strain
		VectorN m_strainValuesPerFace;
		// Velocities of amplitudes
		VectorN m_AmplitudeVelocities;
		std::vector<unsigned int> m_collisionVertices;
		// for self collisions
		Real m_averageEdgeLengths;
		// frame
		unsigned int curFrame;
		// for dphi penalty energy
		std::vector<Real> m_minimum_edge_lengths;
		std::vector<Matrix32> m_dphi_to_gradphis;
		std::pair<bool, std::vector<Matrix22>> m_cachedI;
		std::pair<bool, std::vector<Matrix22>> m_cachedII;
		std::pair<bool, std::vector<Matrix22>> m_cachedDaDaT;
		std::pair<bool, std::vector<Matrix22>> m_cachedDphiDphiT;
		std::pair<bool, std::vector<Matrix22>> m_cachedDaDphiT;
		std::map<int, std::unordered_map<int, Vector3>> m_ClampedAnim;

	public:
		// Sheet takes ownership of positions, but triangles is
		// used to compute half edge structure so can be reused
		// TODO: Check if UVs should be mandatory
		Sheet(const Material& material,
		    std::vector<trimesh::triangle_t>& triangles,
		    argus::Matrix3N& positions,
		    std::vector<trimesh::triangle_t>& uvtriangles, // possible for meshes with seams // TODO: not required remove this
		    argus::Matrix2N& uvs,
		    const SheetInfo& info);

		trimesh::index_t materialPoint(const Vector3& point) const;

		const trimesh::trimesh_t& getHalfEdges() const { return m_HalfEdges; }
		const std::vector<trimesh::edge_t>& getEdges() const override { return m_Edges; }
		const std::vector<trimesh::triangle_t>& getTriangles() const override { return m_Triangles; }
		const std::vector<trimesh::triangle_t>& getUVTriangles() const override { return m_UVTriangles; }
		std::vector<trimesh::triangle_t>& getMTriangles() { return m_Triangles; }
		const argus::Matrix2N& getUVs() const { return m_UVs; }

		void addConstraint(const FixedConstraint& constraint) { m_FixedConstraints.push_back(constraint); }
		void addConstraint(const LineConstraint& constraint);
		void addConstraint(const PlaneConstraint& constraint);
		void addConstraint(const DynamicConstraint& constraint) { m_DynamicConstraints.push_back(constraint); }

		Map<VectorN> getInvInertia() override { return FlattenMatrixN3(m_InvInertia); }
		Map<VectorN> getVelocities() override { return FlattenMatrixN3(m_Velocities); }
		Map<VectorN> getAmpVec() override { return Map<VectorN>(m_Amplitudes.data(), m_VertexCount); }
		Map<Matrix2N> getDphis() override { return Map<Matrix2N>(m_dPhis1VectorPerFace.data(), 2, m_FaceCount); }
		Map<VectorN> getAmpVelocities() override { return Map<VectorN>(m_AmplitudeVelocities.data(), m_VertexCount); }
		Matrix3N& getVelocitiesN3() { return m_Velocities; }
		Map<VectorN> getPositions() override { return FlattenMatrixN3(m_Positions); }
		Map<VectorN> getRestPositions() override { return FlattenMatrixN3(m_RestPositions); }

		Matrix3N& getPositionsN3() { return m_Positions; }
		const Matrix3N& getPositionsN3() const { return m_Positions; }
		VectorN& getInertia() override { return m_Inertia; }
		trimesh::index_t getVertexCount() override { return m_VertexCount; }
		trimesh::index_t getEdgeCount() override { return m_EdgeCount; }
		trimesh::index_t getFaceCount() override { return m_FaceCount; }

		Map<VectorN> getForces() override { return FlattenMatrixN3(m_Material.getForces(*this, m_Forces)); }
		VectorN getWrinkleForces() override
		{
			m_WrinkleForces.resize(4 * m_VertexCount + 2 * m_FaceCount);
			m_WrinkleForces.setZero();
			return m_Material.getForces(*this, m_WrinkleForces);
		}
		Real getWrinkleEnergiesPerFace(int faceId, Eigen::Vector<Real, 23>* deriv = NULL) override
		{
			return m_Material.getEnergyPerFace(*this, deriv, faceId);
		}

		void requireForcesJacobian() override;
		void requireVertexColoringInfo() override;
		ForcesJacobian getForcesJacobian() override;
		void applyConstraints(Map<VectorN> vec) override;
		// dynamic constraints
		void applyDynamicConstraints(Real dt) override;
		void applyTopologicalConstraints() override;
		bool testTopologicalConstraints();
		RestShape& getRestShape() { return m_RestShape; }
		// register clamped vertices
		void registerClampedDofs(std::vector<std::vector<int>> clampeddofs);
		// get clamped dofs
		// std::vector<int>& getClampedDofsIdx() { return m_ClampedDofs; }
		std::vector<int>& getActiveDofsIdx() { return m_ActiveDofs; }
		// get dofs as a vector
		// TODO: Changing this, check if any error is there
		VectorN getDofs(bool active = false) override;
		VectorN getDofsDot(bool active = false) override;
		// update sheet
		void updateDofs(VectorN& dofs, VectorN dofsdot = VectorN()) override;
		// For Wrinkle stuff
		// TODO: Delete used for setting up the system with constraints
		// std::vector<int>& getFreeDOFsIdx() { return m_idxOfDOFs; };
		// std::vector<bool>& getFreeDOFsMap() { return m_freeDOFMap; };
		// std::vector<int>& getFreeDOFInvIdxMap() { return m_freeDOFInvIdxMap; };
		// get vertex degree matrix
		MatrixNN& getVerticesDegree() { return m_VerticesDegree; }
		// TODO: remove as managed by saveCheckpoint
		nlohmann::json& getCheckpointParams() { return checkpointParams; } // contains the data on which parameters were present in the checkpoint info
		// unique checkpoint directory name based on current time
		std::string& getCheckpointDirName() { return checkpointDirName; };
		// nfaces x 3 matrix containing faceid, vi id, vj id, vk id as each row
		const Eigen::MatrixXi& getFaceMatrix() const { return m_faceMatrix; }
		// returns mesh topology used for wrinkle phase reconstruction
		const MeshConnectivity& getBaseMesh() const { return m_baseMesh; }
		// returns the edge id opposite to each vertex of a face in the mesh
		const std::vector<trimesh::triangle_t>& getVertexOppositeEdgesPerFace() const { return m_VtexOppositeEdgesPerFace; } // edgeid opp to vertex per face
		// returns the vertex id opposite to each edge of each face in the mesh
		const Eigen::MatrixXi& getEdgeOppositeVertices() const { return m_EdgeOppositeVertices; }
		// upsampled mesh triangles
		std::vector<trimesh::triangle_t>& getWrinkledMeshTriangles() { return m_WrinkledTriangles; }
		// positions from the upsampled wrinkle mesh
		argus::Matrix3N& getWrinkleMeshPositions() { return m_WrinkleMeshPositions; }
		// upsampled mesh faces matrix
		Eigen::MatrixXi& getWrinkleMeshFaces() { return m_WrinkleMeshFaces; }
		// rest positions
		Matrix3N& getRestPositionsN3() { return m_RestPositions; }
		// amplitudes per vertex
		VectorN& getAmplitudes() { return m_Amplitudes; }
		// amplitudes per face
		VectorN& getAmplitudesPerFace() { return m_ampsPerFace; }
		// dphi one forms per edge
		VectorN& getDPhis() { return m_dPhis; }
		// free degrees of freedom per edge: 0 for mid edge average normal
		VectorN& getEdgeDofs() { return m_edgeDofs; }
		// wrinkle phase per vertex
		VectorN& getPhis() { return m_Phis; }
		Eigen::Vector<std::complex<Real>, Eigen::Dynamic>& getzPhis() { return m_zPhis; }
		VectorN& getPhiVelocities() { return m_PhiVelocities; }
		// phase recovered for each face vertex (nfaces x 3)
		argus::Matrix3N& getCornerPhis() { return m_CornerPhis; }
		// dphis stored as 1 vectors per face
		Matrix2N& getDphis1VectorPerFace() { return m_dPhis1VectorPerFace; }
		// dphis averaged per vertex in reference frame
		Matrix3N& getDphisPerVertex() { return m_dPhisPerVertex; }
		// required for quadratic bending energy model
		std::vector<argus::Real>& getVertAreas() { return m_vertAreas; }
		// quad points for wrinkle energies
		std::vector<QuadraturePoints>& getQuadPoints() { return m_quadPoints; }
		// get amplitude at a point in a triangle face from amplitudes stored per vertex
		Real computeAmplitudesFromQuad(int faceId, int quadId, Vector2* da, Vector3* gradA, Matrix23* gradDA, Matrix33* hessianA, std::vector<Matrix33>* hessianDA);
		// get amplitude at a point in a triangle face from amplitudes stored per face
		Real computeAmplitudesPerFaceFromQuad(int faceId, int quadId, Vector2* da, Real* gradA, Vector2* gradDA, Matrix33* hessianA, std::vector<Matrix33>* hessianDA);
		// get dphi from dphi one form stored per face
		Vector2 computeDphi(int faceId, Matrix23* gradDphi);
		// get dphi stored as one vector per face
		Vector2 computeDphi1VectorPerFace(int faceId, Matrix22* gradDphi);
		// dphi dphiT dphi stored per face
		Matrix22 computeDphiDphi1VectorPerFaceTensor(int faceId, int quadId, std::vector<Matrix22>* deriv, std::vector<Matrix22>* hess = NULL);
		// da dphiT dphi stored per face, amp stored per vertex
		Matrix22 computeDaDphi1VectorPerFaceTensor(int faceId, int quadId, std::vector<Matrix22>* deriv, std::vector<Eigen::Matrix<Real, 5, 5>>* hess = NULL); // dphi daT dphi stored per face, amp stored per vertex
		Matrix22 computeDphi1VectorPerFaceDaTensor(int faceId, int quadId, std::vector<Matrix22>* deriv, std::vector<Eigen::Matrix<Real, 5, 5>>* hess = NULL);
		// da daT amp stored per vertex
		Matrix22 computeDaDaTensor(int faceId, int quadId, std::vector<Matrix22>* deriv, std::vector<Matrix33>* hessian);
		// dphi dphiT dphi stored as one forms per edge
		Matrix22 computeDphiDphiTensor(int faceId, int quadId, std::vector<Matrix22>* deriv, std::vector<Matrix33>* hessian);
		// da dphiT dphib stored as one forms per edge, amp as scalar per vertex
		Matrix22 computeDaDphiTensor(int faceId, int quadId, std::vector<Matrix22>* deriv, std::vector<Eigen::Matrix<Real, 6, 6>>* hessian);
		// dphi daT dphi stored as one forms per edge, amp stored per vertex
		Matrix22 computeDphiDaTensor(int faceId, int quadId, std::vector<Matrix22>* deriv, std::vector<Eigen::Matrix<Real, 6, 6>>* hessian);
		// define the quadrature points (three point quadrature)
		void setupQuadraturePoints();
		// defines edge opposite to vertices of each face
		void setupEdgeOppToVertexPerFaces();
		// name explains itself
		std::vector<argus::Vector3> getFaceNormals();
		// returns first fundamental matrix for each face as a vector of matrices
		std::vector<Matrix22>& getIValues() { return m_IValues; }
		// used for amplitude initialization and wrinkle phase field reconstruction
		void computeLaplacian();
		Eigen::SparseMatrix<argus::Real>& getLaplacian() { return laplacian; }
		// return lame alpha parameter
		Real getLameAlpha() { return m_Material.getLameAlpha(); } // original 3.2967 * 1e4
		Real getYoungsModulus() { return m_Material.getY(); }
		// return lame beta parameter
		Real getLameBeta() { return m_Material.getLameBeta(); } // original 3.85 * 1e4
		// return sheet thickness
		Real getThickness() { return m_Material.getThickness(); } // original 0.1
		// clear sheet forces
		void clearSheetForces() { m_Material.clearSheetForces(); }
		// add sheet forces
		void addSheetForces(std::vector<std::shared_ptr<SheetForces>> forces) { m_Material.addSheetForces(forces); };
		// get sheet density
		Real getDensity() { return m_Material.getDensity(); }
		// get material
		Material& getMaterial() { return m_Material; }
		// reference mesh specific functions
		Matrix3N& getReferenceMeshPositions() { return m_ReferenceMeshPositions; }
		std::vector<trimesh::triangle_t>& getReferenceMeshTriangles() { return m_ReferenceMeshTriangles; }
		void registerReferenceMesh() { m_UseReferenceMesh = true; }
		bool useReferenceMesh() { return m_UseReferenceMesh; }
		// useful for cut mesh setup
		std::map<int, int>& getSeamVertices() { return m_SeamVertices; }
		// checkpointing specific functions
		// need to remove can be done using save_checkpoint
		void addParamToSave(std::string key) { paramsToSave[key] = 1; }
		void removeParamToSave(std::string key) { paramsToSave.erase(key); }
		std::map<std::string, int>& getParamsToSave() { return paramsToSave; }
		// Contains the energy per step, required for saving the checkpoints for visualization
		argus::VectorN& getElasticEnergies() { return m_elasticEnergies; }
		argus::VectorN& getWrinkleEnergies() { return m_wrinkleEnergies; }
		// Contains the eigen value of the strain per face, used for visualization
		VectorN& getStrainValuesPerFace() { return m_strainValuesPerFace; }
		// Returns rate of change of amplitude wrt time
		VectorN& getAmpsVelocities() { return m_AmplitudeVelocities; }
		// returns sheet energies, derivatives and hessians for sheet with no wrinkle parameters
		Real getEnergies(VectorN* deriv = NULL, std::vector<Eigen::Triplet<Real>>* Hess = NULL, bool vbd = false) override { return m_Material.getEnergies(*this, deriv, Hess, vbd); };

		Real getEnergies(VectorN& dofs, VectorN* deriv = NULL, std::vector<Eigen::Triplet<Real>>* Hess = NULL, bool vbd = false) override { return m_Material.getEnergies(*this, dofs, deriv, Hess, vbd); }
		Real getVertexEnergy(int vid, Vector3* deriv = NULL, Matrix33* hess = NULL) override { return m_Material.getVertexEnergy(*this, vid, deriv, hess); }
		Real getVertexEnergyNew(int vid, Vector3* deriv = NULL, Matrix33* hess = NULL) { return m_Material.getVertexEnergyNew(*this, vid, deriv, hess); }
		Real getVertexEnergy4New(int vid, Vector4* deriv = NULL, Matrix44* hess = NULL) { return m_Material.getVertexEnergy4New(*this, vid, deriv, hess); }

		Real getVertexEnergy4(int vid, Vector4* deriv = NULL, Matrix44* hess = NULL) override { return m_Material.getVertexEnergy4(*this, vid, deriv, hess); }
		Real getVertexEnergy1(int vid, Real* deriv = NULL, Real* hess = NULL) override { return m_Material.getVertexEnergy1(*this, vid, deriv, hess); }
		Real getVertexEnergy1New(int vid, Real* deriv = NULL, Real* hess = NULL) { return m_Material.getVertexEnergy1New(*this, vid, deriv, hess); }

		Real getVertexEnergy1(int vid, std::vector<int>& fids, VectorN* deriv = NULL, MatrixNN* hess = NULL) override { return m_Material.getVertexEnergy1(*this, vid, fids, deriv, hess); }
		Real getEnergy() override { return m_Material.getEnergy(*this); }
		Eigen::MatrixXi getFaces() override { return getFaceMatrix(); }
		Real getAverageEdgeLengths() override { return m_averageEdgeLengths; }
		Eigen::SparseMatrix<int> getAdjacencyMatrix() override { return m_adjacencyMat; }
		// fill active dofs
		void setActiveDofs(std::vector<int> activedofs) override;
		// for lbfgs
		Real getE(VectorN& X, VectorN& deriv, std::vector<Eigen::Triplet<Real>>* hess = NULL) override;

		// optimize dphis 1 vector per face parallely
		// void optimizeDphis1VectorPerFace(std::vector<Real>& energies, int max_iter= 1000);
		void optimizeDphis1VectorPerFace(Real& energy, int max_iter, Real tolerance);
		Real getWrinkleShellEnergiesPerFace(int i, VectorN& grad, MatrixNN* hess = NULL) override { return m_Material.getWrinkleShellEnergiesPerFace(*this, i, &grad, hess); }
		Eigen::VectorXi& getClampedVertices() override { return m_clampedVertices; }
		Eigen::MatrixXi& getClampInfo() { return m_clampInfo; }
		void getWrinkleInitializationParams(VectorN& initAmps, VectorN& initDphis);
		const SheetInfo& getSheetInfo() const override
		{
			return sheetInfo;
		}
		void clampDofs(VectorN& X) override;
		// info for vbd
		SheetInfo sheetInfo {};
		// struct
		// {
		// 	std::string sheetPath;
		// 	Eigen::MatrixXi vertexFacesSet;
		// 	Eigen::MatrixXi coloredSets;
		// 	VectorN axt;
		// 	VectorN aat;
		// 	int ndofs;
		// } vbdInfo {};

		void estimateAmpOmegaFromStrain(Real amplitudeEstimate, VectorN& amps, Matrix2N& w, bool estimateAmps);
		void locatePotentialPureTensionFaces(std::set<int>& potentialPureTensionFaces);
		void getInitialWrinkleParameters(Real amplitudeEstimate, VectorN& amp, Matrix2N& w, bool estimateAmps) override { estimateAmpOmegaFromStrain(amplitudeEstimate, amp, w, estimateAmps); }
		void writeCheckpoint(std::vector<std::string> data_to_save, std::string parent_directory, int frame, nlohmann::json state) override { saveCheckpoint(data_to_save, parent_directory, frame, *this, state); }
		void optimizeParallelEnergies() override
		{
			Real e = 0;
			optimizeDphis1VectorPerFace(e, 100, 1e-6);
		}

		VectorN getEdgeConstraints(Sheet& sheet);
		void applyDynamicCollisionConstraintsOnEdges(Sheet& sheet, Matrix3N& oldPositions, Real dt);
		std::vector<std::pair<int, Real>> getSelfCollisionEdges(Sheet& sheet);
		void applyDynamicSelfCollisionConstraintsOnFaces(Sheet& sheet, Matrix3N& oldPositions, Real dt);
		std::vector<std::pair<Vector2, Eigen::Vector<int, 3>>> getConstrainedFaceInfo(Sheet& sheet);
		Real getAlpha(Vector3 vi, Vector3 vj, Vector3 vk);
		void resolveWrinkleSelfCollisions(VectorN& oldPositions, Real dt) override;
		void lock() { m_mutex.lock(); }
		void unlock() { m_mutex.unlock(); }
		std::vector<unsigned int> getCollisionVertices() { return m_collisionVertices; }
		void setCollisionVertices(std::vector<unsigned int> collverts) override { m_collisionVertices = collverts; }
		void computeGreenStrainPerFace();
		Eigen::MatrixXi getVertexFacesSet() override { return sheetInfo.vertexFacesSet; }
		void writeEnergyContour(unsigned int vid, unsigned int fid, unsigned int subframe, MatrixNN stateMat, std::vector<std::pair<VectorN, Matrix2N>> dvals) override;
		void writeEnergyContourBak(unsigned int vid, unsigned int fid, unsigned int dphiid, unsigned int subframe, Real a0, Real a0_1, Real a1, Real h0, Real h1, std::vector<std::pair<VectorN, Matrix2N>> dvals);
		void updateCurFrame(unsigned int f) { curFrame = f; }
		// Real writeEhatEnergy(std::string basedir, unsigned int vid, unsigned int fid, Real a);
		Real getEhatEnergy(unsigned int vid, unsigned int fid, Real a, Vector2 dphi);
		std::vector<Real> getMinimumEdgeLengths() { return m_minimum_edge_lengths; }
		const std::vector<Matrix32>& getDphiToGradPhis() const { return m_dphi_to_gradphis; }
		void updateBaseMesh();
		Matrix22 computeDphiDphi1VectorPerFaceTensorNew(int faceId, int quadId, std::array<Matrix22, 2>* deriv = NULL, std::array<Matrix22, 4>* hess = NULL);
		Matrix22 computeDaDaTensorNew(int faceId, int quadId, std::array<Matrix22, 3>* deriv = NULL, std::array<Matrix33, 4>* hessian = NULL);
		Real computeAmplitudesFromQuadNew(int faceId, int quadId, Vector2* da = NULL, Vector3* gradA = NULL, Matrix23* gradDA = NULL, Matrix33* hessianA = NULL, std::array<Matrix33, 2>* hessianDA = NULL);
		Matrix22 computeDaDphi1VectorPerFaceTensorNew(int faceId, int quadId, std::array<Matrix22, 5>* deriv = NULL, std::array<Eigen::Matrix<Real, 5, 5>, 4>* hess = NULL);
		Matrix22 computeDphi1VectorPerFaceDaTensorNew(int faceId, int quadId, std::array<Matrix22, 5>* deriv = NULL, std::array<Eigen::Matrix<Real, 5, 5>, 4>* hess = NULL);
		const std::pair<bool, std::vector<Matrix22>>& getCachedIs() { return m_cachedI; }
		const std::pair<bool, std::vector<Matrix22>>& getCachedIIs() { return m_cachedII; }
		const std::pair<bool, std::vector<Matrix22>>& getCachedDaDaT() { return m_cachedDaDaT; }
		const std::pair<bool, std::vector<Matrix22>>& getCachedDphiDpiT() { return m_cachedDphiDphiT; }
		const std::pair<bool, std::vector<Matrix22>>& getCachedDaDpiT() { return m_cachedDaDphiT; }
		bool cachedIsValid() { return m_cachedI.first; }
		bool cachedIIsValid() { return m_cachedII.first; }
		bool cachedDaDaTValid() { return m_cachedDaDaT.first; }
		bool cachedDphiDphiTValid() { return m_cachedDphiDphiT.first; }
		bool cachedDaDphiTValid() { return m_cachedDaDphiT.first; }
		void cacheParameters(std::array<bool, 5> parameters) override;
		Vector2 computeDaFromFace(int faceId, Matrix23* gradDA = NULL);
		void clearCache() override;
		void getSelfCollisionEdgesNew(VectorN& collisionInfo);
		void applyWrinkleSelfCollisionConstraints(Real softness = 0);
		void loadAdjacencyLists();
		const std::vector<std::vector<int>>& getVertexFaceAdjacencyList() const override { return m_vertexFaceAdjacencyList; }
		const std::vector<std::vector<int>>& getVertedEdgeAdjacencyList() const override { return m_vertexEdgeAdjacencyList; }
		bool& setWrinkleEnergy() { return m_UseWrinkleParameters; }
		// for laplacian damping
		void buildDampingLaplacian();
		const Eigen::SparseMatrix<Real>& getDampingLaplacian() const override { return m_DampingLaplacian; }
		// delete after usage
		void printDebInfo(int faceId, int quadId) override { m_Material.printDebInfo(*this, faceId, quadId); }
		const Eigen::Vector<Real, Eigen::Dynamic>& getFaceAreas() { return m_faceAreas; }
		void applySmallDphiConstraints();
	};
}

#endif /* SHEET_HPP */
