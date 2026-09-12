#ifndef DYNAMICS_OBJECT_HPP
#define DYNAMICS_OBJECT_HPP

#include "common.hpp"

namespace argus
{
	struct ForcesJacobian
	{
#ifdef BOGUS
		BlockSparseMatrix& forcesPositionDeriv;
		BlockSparseMatrix& forcesVelocityDeriv;
#endif
	};
	class DynamicsObject
	{
	public:
		// TODO: Think whether we should allow copying objects
		DynamicsObject(const DynamicsObject&) = delete;
		DynamicsObject() { };
		virtual ~DynamicsObject() = default;
		// TODO: Remove this as we can treat it as dofs and dofsdot
		virtual Map<VectorN> getPositions() = 0;
		virtual Map<VectorN> getRestPositions() = 0;
		virtual Map<VectorN> getVelocities() = 0;
		virtual Eigen::MatrixXi getFaces() = 0;
		virtual Eigen::SparseMatrix<int> getAdjacencyMatrix() = 0;
		virtual Real getAverageEdgeLengths() = 0;
		virtual Map<VectorN> getAmpVec() = 0;
		virtual Map<VectorN> getAmpVelocities() = 0;
		virtual Map<Matrix2N> getDphis() = 0;
		virtual VectorN& getInertia() = 0;
		virtual Map<VectorN> getInvInertia() = 0;
		virtual const SheetInfo& getSheetInfo() const = 0;
		// todo remove not needed
		virtual Real getEnergies(VectorN* deriv = NULL, std::vector<Eigen::Triplet<Real>>* Hess = NULL, bool vbd = false) = 0;
		virtual Real getEnergies(VectorN& dofs, VectorN* deriv = NULL, std::vector<Eigen::Triplet<Real>>* Hess = NULL, bool vbd = false) = 0;
		virtual Real getE(VectorN& X, VectorN& deriv, std::vector<Eigen::Triplet<Real>>* hess = NULL) = 0;
		virtual Real getVertexEnergy(int vid, Vector3* deriv = NULL, Matrix33* hess = NULL) = 0;
		virtual Real getVertexEnergy4(int vid, Vector4* deriv = NULL, Matrix44* hess = NULL) = 0;
		virtual Real getVertexEnergy1(int vid, Real* deriv = NULL, Real* hess = NULL) = 0;
		virtual Real getVertexEnergy1(int vid, std::vector<int>& fids, VectorN* deriv = NULL, MatrixNN* hess = NULL) = 0;
		virtual Real getWrinkleShellEnergiesPerFace(int i, VectorN& grad, MatrixNN* hess = NULL) = 0;
		virtual Real getEnergy() = 0;
		virtual Map<VectorN> getForces() = 0;
		virtual void requireForcesJacobian() = 0;
		virtual void requireVertexColoringInfo() = 0;
		virtual ForcesJacobian getForcesJacobian() = 0;
		virtual void applyConstraints(Map<VectorN> vec) = 0;
		virtual void applyDynamicConstraints(Real dt) = 0;
		virtual void applyTopologicalConstraints() = 0;
		virtual void setActiveDofs(std::vector<int> activedofs) = 0;
		virtual VectorN getDofs(bool active = false) = 0;
		virtual VectorN getDofsDot(bool active = false) = 0;
		virtual void updateDofs(VectorN& dofs, VectorN dofsdot = VectorN()) = 0;
		virtual trimesh::index_t getVertexCount() = 0;
		virtual trimesh::index_t getFaceCount() = 0;
		virtual trimesh::index_t getEdgeCount() = 0;
		virtual void clampDofs(VectorN& X) = 0;
		virtual void getInitialWrinkleParameters(Real amplitudeEstimate, VectorN& amp, Matrix2N& w, bool estimateAmps) = 0;
		virtual void writeCheckpoint(std::vector<std::string> data_to_save, std::string parent_directory, int frame, nlohmann::json state) = 0;
		virtual void writeEnergyContour(unsigned int vid, unsigned int fid, unsigned int subframe, MatrixNN state, std::vector<std::pair<VectorN, Matrix2N>> dvals) = 0;
		virtual void optimizeParallelEnergies() = 0;
		virtual void resolveWrinkleSelfCollisions(VectorN& oldPositions, Real dt) = 0;
		virtual void setCollisionVertices(std::vector<unsigned int> collverts) = 0;
		virtual Eigen::MatrixXi getVertexFacesSet() = 0;
		virtual void cacheParameters(std::array<bool, 5>) = 0;
		virtual void clearCache() = 0;
		virtual const std::vector<trimesh::triangle_t>& getTriangles() const = 0;
		virtual const std::vector<trimesh::triangle_t>& getUVTriangles() const = 0;
		virtual const std::vector<trimesh::edge_t>& getEdges() const = 0;
		virtual Eigen::VectorXi& getClampedVertices() = 0;
		virtual const std::vector<std::vector<int>>& getVertexFaceAdjacencyList() const = 0;
		virtual const std::vector<std::vector<int>>& getVertedEdgeAdjacencyList() const = 0;
		// delete after use
		virtual void printDebInfo(int faceId, int quadId) = 0;
		virtual VectorN getWrinkleForces() = 0;
		virtual Real getWrinkleEnergiesPerFace(int faceId, Eigen::Vector<Real, 23>* deriv) = 0;
		virtual const Eigen::SparseMatrix<Real>& getDampingLaplacian() const = 0;
	};

}

#endif /* DYNAMICS_OBJECT_HPP */
