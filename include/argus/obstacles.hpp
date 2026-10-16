#ifndef ARGUS_OBSTACLES_HPP
#define ARGUS_OBSTACLES_HPP

#include "common.hpp"
#include "sheet.hpp"
#include "obstacle_helpers.hpp"
#include <nlohmann/json.hpp>
#include <Discregrid/All>
#include <future>
#include <memory>
#include <map>
#include "nanoflann.hpp"

struct SMPLPointCloud
{
	const Eigen::MatrixXd* pts;

	SMPLPointCloud(const Eigen::MatrixXd* pts_)
	    : pts(pts_)
	{
	}

	inline size_t kdtree_get_point_count() const
	{
		return pts->rows();
	}

	inline double kdtree_get_pt(const size_t idx, int dim) const
	{
		return (*pts)(idx, dim);
	}

	template <class BBOX>
	bool kdtree_get_bbox(BBOX&) const { return false; }
};
namespace argus
{
	enum class ObstacleType
	{
		PLANE,
		SPHERE,
		SDF,
		CYLINDER
	};

	static const std::map<std::string, ObstacleType> obstacleMap = {
		{ "PLANE", ObstacleType::PLANE },
		{ "SPHERE", ObstacleType::SPHERE },
		{ "CYLINDER", ObstacleType::CYLINDER },
		{ "SDF", ObstacleType::SDF }
	};

	class Obstacle
	{
	public:
		virtual ~Obstacle() = default;

		// use maybe in future
		// virtual Real const getPenaltyEnergy(Sheet& sheet, VectorN* deriv, std::vector<Eigen::Triplet<Real>>* hess) const { return 0; };
		virtual void project(VectorN& X, VectorN& V, const VectorN& Xold, Real dt, unsigned int frame = 0, std::vector<int>* contactVertices = NULL, unsigned int substep = 0)
		{
			solveCollisionVertex(this, X, V, Xold, dt, getThreshold(), getFrictionK(), substep, contactVertices);
		}
		// virtual void project(Real dt, VectorN& X, VectorN& V, const VectorN& Xold, const VectorN& W, const std::vector<trimesh::triangle_t>& faces, unsigned int frame = 0, std::vector<int>* contactVertices = NULL, nlohmann::json* config = NULL) const = 0;
		virtual void project(Real dt, VectorN& X, VectorN& V, const VectorN& Xold, const VectorN& W, const std::vector<trimesh::triangle_t>& faces, unsigned int frame = 0, std::vector<int>* contactVertices = NULL, nlohmann::json* config = NULL) const
		{
			int probe_points = 10;
			if (config && config->contains("probe_points"))
			{
				probe_points = (*config)["probe_points"];
			}
			probeWeights weights;
			switch (probe_points)
			{
			case 10:
			{
				weights = weights10;
				break;
			}
			case 15:
			{
				weights = weights15;
				break;
			}
			case 25:
			{
				weights = weights25;
				break;
			}
			default:
			{
				weights = weights10;
			}
			}

			solveCollisionQuadrature(this, dt, X, V, Xold, W, faces, getThreshold(), const_cast<Obstacle*>(this)->getFrictionK(), weights, contactVertices, probe_points);
		}
		virtual Real getSDF(const Vector3& p, Vector3* gradient = nullptr) const = 0;
		virtual bool hasExplicitVelocity() const { return false; }
		virtual Vector3 getVelocity(const Vector3& p) const { return Vector3::Zero(); }

		virtual Real getFrictionK() = 0;
		virtual void setTransform(const Vector3& t, const Eigen::Quaternion<Real>& r, const Vector3& s, Real dt) = 0;
		virtual void updateFrame() = 0;
		virtual void startAsyncLoad(int frameIdx) = 0;
		virtual void updateToActive() = 0;
		virtual std::string getName() = 0;
		virtual std::vector<int> getHitVertices() const { return {}; }
		virtual void exportToVTK(const std::string& path) const { }
		virtual void freezeMesh(bool freeze) { }
		virtual const Real getThreshold() const = 0;
		virtual bool inBounds(const Vector3& p, Real threshold) const { return true; }
		virtual bool meshAvailable() = 0;

		virtual std::tuple<Eigen::Matrix<Real, Eigen::Dynamic, 3>, std::vector<std::vector<trimesh::index_t>>> getMesh() const { throw "Not implemented!"; }
	};

	class Plane : public Obstacle
	{
		struct
		{
			Vector3 center;
			Vector3 axes[3];
			Vector3 extents;
		} OBB;
		Vector3 m_baseExtents;
		Vector3 m_v0, m_v1, m_v2, m_v3;
		Real m_threshold, m_Kpenalty, m_thickness;

		Real m_Kfriction;

		std::string m_name;
		Vector3 m_linearVelocity = Vector3::Zero();
		Vector3 m_angularVelocity = Vector3::Zero();
		Vector3 m_center_prev;
		bool m_hasExplicitVelocity;
		Eigen::Quaternion<Real> m_rotation;

	public:
		Plane(nlohmann::json& data) { build(data); }
		Real getSDF(const Vector3& p, Vector3* gradient = nullptr) const override;
		bool hasExplicitVelocity() const override { return true; }
		Vector3 getVelocity(const Vector3& p) const override;

		Real getFrictionK() override { return m_Kfriction; }
		bool meshAvailable() override { return true; }
		void setTransform(const Vector3& t, const Eigen::Quaternion<Real>& r, const Vector3& s, Real dt) override;
		void updateFrame() override;

		void startAsyncLoad(int frameIdx) override { }
		void updateToActive() override { }
		std::string getName() override { return m_name; }
		const Real getThreshold() const override { return m_threshold; }
		bool inBounds(const Vector3& p, Real threshold) const override;
		std::tuple<Eigen::Matrix<Real, Eigen::Dynamic, 3>, std::vector<std::vector<trimesh::index_t>>> getMesh() const override;

	private:
		void build(nlohmann::json& data);
		void updateOBB();
		bool isInsideOBB(Vector3 pos) const;
	};

	class Sphere : public Obstacle
	{
		Vector3 m_center, m_center_prev;
		Real m_radius, m_threshold, m_Kpenalty;
		Real x_min, y_min, z_min, x_max, y_max, z_max;

		Real m_Kfriction;

		std::string m_name;
		Vector3 m_linearVelocity = Vector3::Zero();
		Vector3 m_angularVelocity = Vector3::Zero();
		bool m_hasExplicitVelocity = false;
		Eigen::Quaternion<Real> m_rotation;

	public:
		Sphere(nlohmann::json& data) { build(data); }
		bool hasExplicitVelocity() const override { return true; }
		Vector3 getVelocity(const Vector3& p) const override;

		Real getFrictionK() override { return m_Kfriction; }

		Real getSDF(const Vector3& p, Vector3* gradient = nullptr) const override;
		bool meshAvailable() override { return true; }

		void setTransform(const Vector3& t, const Eigen::Quaternion<Real>& r, const Vector3& s, Real dt) override;
		void updateFrame() override;

		void startAsyncLoad(int frameIdx) override { }
		void updateToActive() override { }
		std::string getName() override { return m_name; }
		const Real getThreshold() const override { return m_threshold; }
		bool inBounds(const Vector3& p, Real threshold) const override;
		std::tuple<Eigen::Matrix<Real, Eigen::Dynamic, 3>, std::vector<std::vector<trimesh::index_t>>> getMesh() const override;

	private:
		void build(nlohmann::json& data);
		void updateAABB();
		bool isInsideAABB(Vector3 pos) const;
	};

	class Cylinder : public Obstacle
	{
		Vector3 m_center, m_axis, m_center_prev, m_axis_prev, m_center_init, m_axis_init, m_local_offset;
		Real m_radius, m_threshold, m_height, m_Kpenalty;

		Real m_Kfriction;

		std::string m_name;
		Vector3 m_linearVelocity = Vector3::Zero();
		Vector3 m_angularVelocity = Vector3::Zero();
		bool m_hasExplicitVelocity = false;
		Real x_min, x_max, y_min, y_max, z_min, z_max;
		Eigen::Quaternion<Real> m_rotation;

	public:
		Cylinder(nlohmann::json& data) { build(data); }
		bool meshAvailable() override { return true; }

		bool hasExplicitVelocity() const override { return true; }
		Vector3 getVelocity(const Vector3& p) const override;

		Real getFrictionK() override { return m_Kfriction; }

		Real getSDF(const Vector3& p, Vector3* gradient = nullptr) const override;

		void setTransform(const Vector3& t, const Eigen::Quaternion<Real>& r, const Vector3& s, Real dt) override;
		void updateFrame() override;

		void startAsyncLoad(int frameIdx) override { }
		void updateToActive() override { }
		std::string getName() override { return m_name; }
		const Real getThreshold() const override { return m_threshold; }
		std::tuple<Eigen::Matrix<Real, Eigen::Dynamic, 3>, std::vector<std::vector<trimesh::index_t>>> getMesh() const override;

	private:
		void build(nlohmann::json& data);
		void updateAABB();
		bool isInsideAABB(Vector3 pos) const;
	};

	class ObstacleDiscregrid : public Obstacle
	{
		std::string m_name, m_basePath;
		Real m_threshold, m_Kpenalty;
		Real m_Kfriction;
		int m_total_frames;
		bool m_isFrozen = false;
		std::vector<trimesh::triangle_t> m_faces;
		Matrix3N m_vertices;
		Vector3 m_translation = Vector3::Zero();
		Eigen::Quaternion<Real> m_rotation = Eigen::Quaternion<Real>::Identity();
		struct FrameData
		{
			std::shared_ptr<Discregrid::CubicLagrangeDiscreteGrid> sdf;
			Eigen::Matrix<Real, 3, Eigen::Dynamic> vertices;
			std::vector<trimesh::triangle_t> faces;
			bool isValid = false;
			int frameIdx = -1;
		};
		FrameData m_currentData;
		std::future<FrameData> m_loadingFuture;
		bool m_loadMesh = true;

	public:
		ObstacleDiscregrid(nlohmann::json& data) { build(data); }
		bool meshAvailable() override { return true; }

		Real getSDF(const Vector3& p, Vector3* gradient = nullptr) const override;
		bool hasExplicitVelocity() const override { return false; }
		Real getFrictionK() override { return m_Kfriction; }

		void setTransform(const Vector3& t, const Eigen::Quaternion<Real>& r, const Vector3& s, Real dt) override { };

		void updateFrame() override { }
		void startAsyncLoad(int frameIdx) override;
		void updateToActive() override;

		void freezeMesh(bool freeze) override { m_isFrozen = freeze; }
		std::string getName() override { return m_name; }
		const Real getThreshold() const override { return m_threshold; }

		std::tuple<Eigen::Matrix<Real, Eigen::Dynamic, 3>, std::vector<std::vector<trimesh::index_t>>> getMesh() const override;

	private:
		void build(nlohmann::json& data);
		static FrameData loadFrameData(const std::string& basePath, int frameIdx, bool loadObj);
	};

	class ObstacleSSDF : public Obstacle
	{
		struct FrameData
		{
			std::shared_ptr<SDFVolume> volume;
			Eigen::Matrix<Real, 3, Eigen::Dynamic> vertices;
			std::vector<trimesh::triangle_t> faces;
			bool isValid = false;
			int frameIdx = -1;
		};

		FrameData m_currentData;
		std::future<FrameData> m_loadingFuture;

		std::string m_name;
		std::string m_basePath;
		Real m_threshold;
		Real m_Kfriction;
		int m_total_frames;
		bool m_isFrozen = false;
		bool m_loadMesh = true;

	public:
		ObstacleSSDF(nlohmann::json& data) { build(data); }
		bool meshAvailable() override { return m_loadMesh; }

		Real getSDF(const Vector3& p, Vector3* gradient = nullptr) const override;

		bool hasExplicitVelocity() const override { return true; }
		Vector3 getVelocity(const Vector3& p) const override;

		Real getFrictionK() override { return m_Kfriction; }

		void startAsyncLoad(int frameIdx) override;
		void updateToActive() override;

		std::tuple<Eigen::Matrix<Real, Eigen::Dynamic, 3>, std::vector<std::vector<trimesh::index_t>>> getMesh() const override;

		void setTransform(const Vector3& t, const Eigen::Quaternion<Real>& r, const Vector3& s, Real dt) override { }
		void updateFrame() override { }

		void freezeMesh(bool freeze) override { m_isFrozen = freeze; }
		std::string getName() override { return m_name; }
		const Real getThreshold() const override { return m_threshold; }

		void exportToVTK(const std::string& path) const override;

	private:
		void build(nlohmann::json& data);

		static FrameData loadFrameData(const std::string& basePath, int frameIdx, bool loadObj);

		Vector3 computeFiniteDifferenceNormal(const Vector3& p) const;
	};
	/*
	    namespace Discregrid
	    {
	        class CubicLagrangeDiscreteGrid;
	    }

	    class ObstacleDiscregrid : public Obstacle
	    {
	        std::string m_name, m_path;
	        Real m_threshold, m_Kpenalty;

	        Real m_Kfriction;

	        std::shared_ptr<Discregrid::CubicLagrangeDiscreteGrid> sdf;
	        bool m_isFrozen = false;

	    public:
	        ObstacleDiscregrid(nlohmann::json& data);
	        Real getSDF(const Vector3& p, Vector3* gradient = nullptr) const override;

	        bool hasExplicitVelocity() const override { return false; }

	        Real getFrictionK() override { return m_Kfriction; }

	        void setTransform(const Vector3& t, const Eigen::Quaternion<Real>& r, const Vector3& s) override { }
	        void updateFrame() override { }
	        void startAsyncLoad(int frameIdx) override { }
	        void updateToActive() override { }
	        void freezeMesh(bool freeze) override { m_isFrozen = freeze; }
	        std::string getName() override { return m_name; }
	        const Real getThreshold() const override { return m_threshold; }

	    private:
	        void build(nlohmann::json& data);
	    };

	    class ObstacleSSDF : public Obstacle
	    {
	        std::string m_name, m_sdfBaseDir;
	        std::string m_currentRenderPath;
	        Real m_threshold;

	        Real m_frictionK;

	        float m_cellSize;
	        bool m_isFrozen = false;
	        bool m_isLoading = false;

	        std::unique_ptr<SDFVolume> sdfobject;
	        std::unique_ptr<SDFVolume> sdfobjectnext;
	        std::future<std::unique_ptr<SDFVolume>> m_loadingTask;
	        std::string m_nextRenderPath;

	    public:
	        ObstacleSSDF(const nlohmann::json& config);
	        Real getSDF(const Vector3& p, Vector3* gradient = nullptr) const override;

	        bool hasExplicitVelocity() const override { return true; }
	        Vector3 getVelocity(const Vector3& p) const override;

	        Real getFrictionK() override { return m_frictionK; }

	        void setTransform(const Vector3& t, const Eigen::Quaternion<Real>& r, const Vector3& s) override { }
	        void updateFrame() override { }

	        void startAsyncLoad(int frameIdx) override;
	        void updateToActive() override;
	        void exportToVTK(const std::string& path) const override;
	        void freezeMesh(bool freeze) override { m_isFrozen = freeze; }

	        std::string getName() override { return m_name; }
	        const Real getThreshold() const override { return m_threshold; }

	        std::string getRenderMeshPath() const { return m_currentRenderPath; }

	    private:
	        float sampleSDF(const Vector3& p) const;
	        Vector3 computeNormal(const Vector3& p) const;
	        void stopAndFreeze();
	    };
	    */
	class ObstacleSMPLSDF : public Obstacle
	{
	public:
		ObstacleSMPLSDF(nlohmann::json& data) { build(data); }

		void build(nlohmann::json& data);

		Real getSDF(const Vector3& p, Vector3* gradient = nullptr) const override;

		Vector3 getVelocity(const Vector3& p) const override;

		bool inBounds(const Vector3& p, Real threshold) const override;

		// 🔥 IMPORTANT: per-frame update
		void updateTransforms(const std::vector<Eigen::Matrix4d>& T);

		Real getFrictionK() override { return m_frictionK; };
		bool meshAvailable() override;
		void setTransform(const Vector3& t, const Eigen::Quaternion<Real>& r, const Vector3& s, Real dt) override { };
		void updateFrame() override { };
		void startAsyncLoad(int frameIdx) override;
		void updateToActive() override;
		std::string getName() override { return "SMPL SDF"; };
		const Real getThreshold() const override { return m_threshold; };
		struct FrameData
		{
			int frameIdx = -1;
			std::vector<Eigen::Matrix4d> transforms;
			bool isValid = false;
		};

		// 🔥 DECLARE FUNCTION
		FrameData loadFrameData(const std::string& basePath, int frameIdx);
		FrameData m_currentData;
		std::future<FrameData> m_loadingFuture;

		std::string m_basePath;
		int m_total_frames = 1;

		std::tuple<
		    Eigen::Matrix<Real, Eigen::Dynamic, 3>,
		    std::vector<std::vector<trimesh::index_t>>>
		getMesh() const override;

	private:
		using KDTree = nanoflann::KDTreeSingleIndexAdaptor<
		    nanoflann::L2_Simple_Adaptor<double, SMPLPointCloud>,
		    SMPLPointCloud,
		    3>;
		std::vector<int> m_dominantJoint;
		std::unique_ptr<KDTree> m_kdtree;
		std::unique_ptr<SMPLPointCloud> m_cloud;
		Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor> m_faces;
		// --- Rest SDF ---
		std::vector<Vector3> m_restJoints;
		std::shared_ptr<SDFVolume> m_restSDF;
		Real m_frictionK, m_threshold;
		// --- SMPL data ---
		Eigen::MatrixXd m_restVertices; // (N × 3)
		std::vector<float> m_weights; // (N × J)

		std::vector<Eigen::Matrix4d> m_T;
		std::vector<Eigen::Matrix4d> m_T_inv;

		int m_numVertices = 0;
		int m_numJoints = 0;

		// --- Core helpers (we’ll implement later) ---
		int findClosestVertex(const Vector3& p) const;

		Vector3 inverseSkinning(const Vector3& p) const;
	};
	std::optional<std::unique_ptr<Obstacle>> getObstacle(const std::string& obsName, nlohmann::json& data);

}
#endif