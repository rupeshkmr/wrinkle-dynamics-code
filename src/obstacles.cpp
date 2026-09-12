#include "halfedge/trimesh_types.h"
#include <Eigen/Geometry>
#include <argus/obstacles.hpp>
#include <argus/util.hpp>
#include <filesystem>
namespace argus
{
	std::optional<std::unique_ptr<Obstacle>> getObstacle(const std::string& obsName, nlohmann::json& data)
	{
		auto it = obstacleMap.find(obsName);
		if (it == obstacleMap.end())
		{
			std::cerr << "Error: Obstacle type '" << obsName << "' not found in map." << std::endl;
			return std::nullopt;
		}

		ObstacleType obsT = it->second;

		try
		{
			switch (obsT)
			{
			case ObstacleType::PLANE:
				return std::make_unique<Plane>(data);
			case ObstacleType::SPHERE:
				return std::make_unique<Sphere>(data);
			case ObstacleType::CYLINDER:
				return std::make_unique<Cylinder>(data);
			case ObstacleType::SDF:
			{
				std::string sdftype = "";
				if (data.contains("sdftype"))
				{
					sdftype = data["sdftype"];
				}

				if (sdftype == "dg")
				{
					return std::make_unique<ObstacleDiscregrid>(data);
				}
				else if (sdftype == "ssdf")
				{
					return std::make_unique<ObstacleSSDF>(data);
				}
				else if (sdftype == "smpl")
				{
					return std::make_unique<ObstacleSMPLSDF>(data);
				}
				else
				{
					std::cerr << "Exception creating obstacle '" << obsName << "': " << " sdf support for this extension is not present!!" << std::endl;
					return std::nullopt;
				}
			}
			default:
				std::cerr << "Error: Unhandled ObstacleType enum for " << obsName << std::endl;
				return std::nullopt;
			}
		}
		catch (const std::exception& e)
		{
			std::cerr << "Exception creating obstacle '" << obsName << "': " << e.what() << std::endl;
			return std::nullopt;
		}
		catch (const char* exp)
		{
			std::cerr << "Error creating obstacle '" << obsName << "': " << exp << std::endl;
			return std::nullopt;
		}
		catch (...)
		{
			std::cerr << "Unknown error creating obstacle '" << obsName << "'" << std::endl;
			return std::nullopt;
		}
	}

	Real Plane::getSDF(const Vector3& p, Vector3* gradient) const
	{
		Vector3 d = p - OBB.center;

		Vector3 q;
		q.x() = OBB.axes[0].dot(d);
		q.y() = OBB.axes[1].dot(d);
		q.z() = OBB.axes[2].dot(d);

		Real clampedX = std::max(-OBB.extents.x(), std::min(q.x(), OBB.extents.x()));
		Real clampedY = std::max(-OBB.extents.y(), std::min(q.y(), OBB.extents.y()));

		Vector3 closestOnRect(clampedX, clampedY, 0.0);

		Vector3 v = q - closestOnRect;
		Real distToSpine = v.norm();

		Real finalDist = distToSpine - OBB.extents.z();

		if (gradient)
		{
			if (distToSpine > 1e-9)
			{
				Vector3 localNormal = v / distToSpine;

				*gradient = (localNormal.x() * OBB.axes[0] + localNormal.y() * OBB.axes[1] + localNormal.z() * OBB.axes[2]); // Normalized by definition of axes
			}
			else
			{
				*gradient = OBB.axes[2];
			}
		}

		return finalDist;
	}

	Vector3 Plane::getVelocity(const Vector3& p) const
	{
		return m_linearVelocity + m_angularVelocity.cross(p - OBB.center);
	}

	void Plane::setTransform(const Vector3& t, const Eigen::Quaternion<Real>& r, const Vector3& s, Real dt)
	{
		Vector3 old_center = OBB.center;
		OBB.center = t;
		m_rotation = r;

		Eigen::Matrix3d R = r.toRotationMatrix();
		OBB.axes[0] = R.col(0);
		OBB.axes[1] = R.col(1);
		OBB.axes[2] = R.col(2);

		const Vector3& axisU = OBB.axes[0];
		const Vector3& axisV = OBB.axes[1];
		Real extU = OBB.extents(0);
		Real extV = OBB.extents(1);
		m_v0 = OBB.center - axisU * extU - axisV * extV;
		m_v1 = OBB.center + axisU * extU - axisV * extV;
		m_v2 = OBB.center + axisU * extU + axisV * extV;
		m_v3 = OBB.center - axisU * extU + axisV * extV;

		if (dt > 1e-9)
		{
			m_linearVelocity = (OBB.center - old_center) / dt;
			m_hasExplicitVelocity = true;
		}
	}

	void Plane::updateFrame()
	{
		m_center_prev = OBB.center;
	}

	bool Plane::isInsideOBB(Vector3 pos) const
	{
		Vector3 vec = pos - OBB.center;
		Real dist0 = vec.dot(OBB.axes[0]);
		if (std::abs(dist0) >= OBB.extents(0))
			return false;
		Real dist1 = vec.dot(OBB.axes[1]);
		if (std::abs(dist1) >= OBB.extents(1))
			return false;
		Real dist2 = vec.dot(OBB.axes[2]);
		if (std::abs(dist2) >= OBB.extents(2))
			return false;
		return true;
	}

	bool Plane::inBounds(const Vector3& p, Real threshold) const
	{
		Vector3 d = p - OBB.center;
		Real x = std::abs(OBB.axes[0].dot(d));
		Real y = std::abs(OBB.axes[1].dot(d));
		Real z = std::abs(OBB.axes[2].dot(d));

		if (x > OBB.extents[0] + threshold)
			return false;
		if (y > OBB.extents[1] + threshold)
			return false;
		if (z > OBB.extents[2] + threshold)
			return false;

		return true;
	}

	void Plane::build(nlohmann::json& data)
	{
		m_threshold = data.value("threshold", 0.01);
		m_rotation = Eigen::Quaternion<Real>::Identity();
		m_Kpenalty = data.value("Kpenalty", 1000.0);
		m_thickness = data.value("thickness", 0.1);
		m_name = data.value("name", "Unnamed_Plane");

		if (data.contains("Kfriction"))
		{
			m_Kfriction = data["Kfriction"].get<Real>();
		}
		else
		{
			m_Kfriction = 0;
		}

		m_linearVelocity = Vector3::Zero();
		m_angularVelocity = Vector3::Zero();
		m_hasExplicitVelocity = false;

		if (data.contains("linear_velocity"))
		{
			m_linearVelocity = util::getVector3(data["linear_velocity"]).value_or(Vector3::Zero());
		}

		if (data.contains("angular_velocity"))
		{
			m_angularVelocity = util::getVector3(data["angular_velocity"]).value_or(Vector3::Zero());
		}

		if (m_linearVelocity.squaredNorm() > 1e-9 || m_angularVelocity.squaredNorm() > 1e-9)
		{
			m_hasExplicitVelocity = true;
		}

		auto v0 = util::getVector3(data["v0"]);
		if (!v0)
			throw std::runtime_error("Plane Error (" + m_name + "): Missing 'v0' in configuration.");
		m_v0 = v0.value();

		auto v1 = util::getVector3(data["v1"]);
		if (!v1)
			throw std::runtime_error("Plane Error (" + m_name + "): Missing 'v1' in configuration.");
		m_v1 = v1.value();

		auto v2 = util::getVector3(data["v2"]);
		if (!v2)
			throw std::runtime_error("Plane Error (" + m_name + "): Missing 'v2' in configuration.");
		m_v2 = v2.value();

		auto v3 = util::getVector3(data["v3"]);
		if (!v3)
			throw std::runtime_error("Plane Error (" + m_name + "): Missing 'v3' in configuration.");
		m_v3 = v3.value();

		updateOBB();
		m_baseExtents = OBB.extents;

		m_center_prev = OBB.center;
		m_center_prev = OBB.center;
	}

	void Plane::updateOBB()
	{
		OBB.center = (m_v0 + m_v1 + m_v2 + m_v3) / 4.0;
		Vector3 eij = m_v1 - m_v0;
		Vector3 eik = m_v3 - m_v0;
		Vector3 crossP = eij.cross(eik);
		Vector3 n = (crossP.norm() > 1e-6) ? crossP.normalized() : Vector3::UnitZ();

		OBB.axes[0] = eij.normalized();
		OBB.axes[1] = eik.normalized();
		OBB.axes[2] = n;
		OBB.extents(0) = 0.5 * eij.norm();
		OBB.extents(1) = 0.5 * eik.norm();
		OBB.extents(2) = m_thickness;
	}

	std::tuple<Eigen::Matrix<Real, Eigen::Dynamic, 3>, std::vector<std::vector<trimesh::index_t>>> Plane::getMesh() const
	{

		Vector3 edge1 = m_v1 - m_v0;
		Vector3 edge2 = m_v2 - m_v0;
		Vector3 normal = edge1.cross(edge2).normalized();

		Vector3 offset = -normal * m_thickness;

		Eigen::Matrix<Real, Eigen::Dynamic, 3> V(8, 3);

		V.row(0) = m_v0;
		V.row(1) = m_v1;
		V.row(2) = m_v2;
		V.row(3) = m_v3;

		V.row(4) = m_v0 + offset;
		V.row(5) = m_v1 + offset;
		V.row(6) = m_v2 + offset;
		V.row(7) = m_v3 + offset;

		std::vector<std::vector<trimesh::index_t>> F;

		auto addQuad = [&](int a, int b, int c, int d)
		{
			F.push_back({ (trimesh::index_t)a, (trimesh::index_t)b, (trimesh::index_t)c });
			F.push_back({ (trimesh::index_t)a, (trimesh::index_t)c, (trimesh::index_t)d });
		};

		addQuad(0, 1, 2, 3);

		addQuad(4, 7, 6, 5);

		addQuad(0, 4, 5, 1);

		addQuad(1, 5, 6, 2);

		addQuad(2, 6, 7, 3);

		addQuad(3, 7, 4, 0);

		return { V, F };
	}

	void Sphere::build(nlohmann::json& data)
	{
		m_threshold = data.value("threshold", 0.01);
		m_Kpenalty = data.value("Kpenalty", 1000.0);
		m_radius = data["radius"];
		m_name = data.value("name", "Unnamed_Sphere");
		m_rotation = Eigen::Quaternion<Real>::Identity();
		if (data.contains("Kfriction"))
		{
			m_Kfriction = data["Kfriction"].get<Real>();
		}
		else
		{
			m_Kfriction = 0;
		}

		auto c = util::getVector3(data["center"]);
		if (c.has_value())
		{
			m_center = c.value();
		}
		else
		{
			throw std::runtime_error("Sphere Error: Missing 'center' in configuration.");
		}

		m_linearVelocity = Vector3::Zero();
		m_angularVelocity = Vector3::Zero();
		m_hasExplicitVelocity = false;

		if (data.contains("linear_velocity"))
		{
			m_linearVelocity = util::getVector3(data["linear_velocity"]).value_or(Vector3::Zero());
		}
		if (data.contains("angular_velocity"))
		{
			m_angularVelocity = util::getVector3(data["angular_velocity"]).value_or(Vector3::Zero());
		}

		if (m_linearVelocity.squaredNorm() > 1e-9 || m_angularVelocity.squaredNorm() > 1e-9)
		{
			m_hasExplicitVelocity = true;
		}

		updateAABB();
		m_center_prev = m_center;
	}

	Real Sphere::getSDF(const Vector3& p, Vector3* gradient) const
	{
		Vector3 d = p - m_center;
		Real distToCenter = d.norm();
		Real finalDist = distToCenter - m_radius;

		if (gradient)
		{
			if (distToCenter > 1e-9)
			{
				*gradient = d / distToCenter;
			}
			else
			{
				*gradient = Vector3(0, 0, 1);
			}
		}
		return finalDist;
	}

	Vector3 Sphere::getVelocity(const Vector3& p) const
	{
		return m_linearVelocity + m_angularVelocity.cross(p - m_center);
	}

	void Sphere::setTransform(const Vector3& t, const Eigen::Quaternion<Real>& r, const Vector3& s, Real dt)
	{
		Vector3 old_center = m_center;
		Eigen::Quaternion<Real> old_rotation = m_rotation;

		m_center = t;
		m_rotation = r;

		if (dt > 1e-9)
		{
			m_linearVelocity = (m_center - old_center) / dt;

			Eigen::Quaternion<Real> q_diff = m_rotation * old_rotation.inverse();
			if (q_diff.w() < 0)
			{
				q_diff.w() = -q_diff.w();
				q_diff.vec() = -q_diff.vec();
			}
			m_angularVelocity = (q_diff.vec() * 2.0) / dt;

			m_hasExplicitVelocity = true;
		}
		updateAABB();
	}

	void Sphere::updateFrame()
	{
		m_center_prev = m_center;
	}

	void Sphere::updateAABB()
	{
		Real buffer = m_radius + m_threshold;
		x_min = m_center(0) - buffer;
		x_max = m_center(0) + buffer;
		y_min = m_center(1) - buffer;
		y_max = m_center(1) + buffer;
		z_min = m_center(2) - buffer;
		z_max = m_center(2) + buffer;
	}

	bool Sphere::inBounds(const Vector3& p, Real threshold) const
	{

		if (p(0) < x_min - threshold || p(0) > x_max + threshold)
			return false;
		if (p(1) < y_min - threshold || p(1) > y_max + threshold)
			return false;
		if (p(2) < z_min - threshold || p(2) > z_max + threshold)
			return false;

		return true;
	}
	std::tuple<Eigen::Matrix<Real, Eigen::Dynamic, 3>, std::vector<std::vector<trimesh::index_t>>> Sphere::getMesh() const
	{
		int stacks = 20;
		int slices = 20;

		int numVertices = 2 + (stacks - 1) * slices;
		Eigen::Matrix<Real, Eigen::Dynamic, 3> V(numVertices, 3);

		auto transform = [&](const Vector3& local) -> Vector3
		{
			return m_center + m_rotation * local;
		};

		V.row(0) = transform(Vector3(0, 0, m_radius));

		int vIdx = 1;
		for (int i = 1; i < stacks; ++i)
		{
			Real phi = M_PI * (Real)i / (Real)stacks; // 0 to Pi
			for (int j = 0; j < slices; ++j)
			{
				Real theta = 2.0 * M_PI * (Real)j / (Real)slices; // 0 to 2Pi

				Real x = m_radius * std::sin(phi) * std::cos(theta);
				Real y = m_radius * std::sin(phi) * std::sin(theta);
				Real z = m_radius * std::cos(phi);

				V.row(vIdx++) = transform(Vector3(x, y, z));
			}
		}

		int southPoleIdx = vIdx;
		V.row(southPoleIdx) = transform(Vector3(0, 0, -m_radius));

		std::vector<std::vector<trimesh::index_t>> F;

		for (int j = 0; j < slices; ++j)
		{
			trimesh::index_t next = (j + 1) % slices;
			F.push_back({ 0, (trimesh::index_t)(1 + j), (trimesh::index_t)(1 + next) });
		}

		for (int i = 0; i < stacks - 2; ++i)
		{
			trimesh::index_t rowStart = 1 + i * slices;
			trimesh::index_t nextRowStart = rowStart + slices;

			for (int j = 0; j < slices; ++j)
			{
				trimesh::index_t next = (j + 1) % slices;

				trimesh::index_t v0 = rowStart + j;
				trimesh::index_t v1 = rowStart + next;
				trimesh::index_t v2 = nextRowStart + next;
				trimesh::index_t v3 = nextRowStart + j;

				F.push_back({ v0, v1, v2 });
				F.push_back({ v0, v2, v3 });
			}
		}

		trimesh::index_t lastRowStart = 1 + (stacks - 2) * slices;
		for (int j = 0; j < slices; ++j)
		{
			trimesh::index_t next = (j + 1) % slices;
			F.push_back({ (trimesh::index_t)southPoleIdx, (trimesh::index_t)(lastRowStart + next), (trimesh::index_t)(lastRowStart + j) });
		}

		return { V, F };
	}

	void Cylinder::build(nlohmann::json& data)
	{
		m_threshold = data.value("threshold", 0.01);
		m_Kpenalty = data.value("Kpenalty", 1000.0);
		m_radius = data["radius"];
		m_height = data["height"];
		m_name = data.value("name", "Unnamed_Cylinder");
		m_rotation = Eigen::Quaternion<Real>::Identity();
		if (data.contains("Kfriction"))
		{
			m_Kfriction = data["Kfriction"].get<Real>();
		}
		else
		{
			m_Kfriction = 0;
		}

		auto c = util::getVector3(data["center"]);
		if (c.has_value())
			m_center = c.value();
		else
			throw std::runtime_error("Cylinder Error: Missing 'center'.");
		m_local_offset = m_center;
		auto a = util::getVector3(data["axis"]);
		if (a.has_value())
			m_axis = a.value().normalized();
		else
			throw std::runtime_error("Cylinder Error: Missing 'axis'.");

		m_center_init = m_center;
		m_axis_init = m_axis;

		m_linearVelocity = Vector3::Zero();
		m_angularVelocity = Vector3::Zero();
		m_hasExplicitVelocity = false;

		if (data.contains("linear_velocity"))
		{
			m_linearVelocity = util::getVector3(data["linear_velocity"]).value_or(Vector3::Zero());
		}
		if (data.contains("angular_velocity"))
		{
			m_angularVelocity = util::getVector3(data["angular_velocity"]).value_or(Vector3::Zero());
		}

		if (m_linearVelocity.squaredNorm() > 1e-9 || m_angularVelocity.squaredNorm() > 1e-9)
		{
			m_hasExplicitVelocity = true;
		}

		updateAABB();
		m_center_prev = m_center;
		m_axis_prev = m_axis;
	}

	Real Cylinder::getSDF(const Vector3& p, Vector3* gradient) const
	{
		Vector3 d = p - m_center;

		Real y = d.dot(m_axis);

		Vector3 radialVec = d - y * m_axis;
		Real x = radialVec.norm();

		Real d_side = x - m_radius;
		Real d_cap = std::abs(y) - (m_height * 0.5);

		Real outsideDist = std::sqrt(std::pow(std::max(d_side, 0.0), 2) + std::pow(std::max(d_cap, 0.0), 2));

		Real insideDist = std::min(std::max(d_side, d_cap), 0.0);

		Real finalDist = outsideDist + insideDist;

		if (gradient)
		{
			if (finalDist > 0 || insideDist < 0)
			{
				Vector3 n_radial;
				if (x > 1e-9)
					n_radial = (radialVec / x);
				else
					n_radial = Vector3::Zero();

				Vector3 n_axial = (y > 0) ? m_axis : -m_axis;

				if (finalDist > 0)
				{
					Real w_side = std::max(d_side, 0.0);
					Real w_cap = std::max(d_cap, 0.0);
					*gradient = (n_radial * w_side + n_axial * w_cap).normalized();
				}
				else
				{
					if (d_side > d_cap)
						*gradient = n_radial;
					else
						*gradient = n_axial;
				}
			}
			else
			{
				if (d_side > d_cap)
					*gradient = (radialVec.norm() > 1e-9) ? radialVec.normalized() : Vector3(1, 0, 0);
				else
					*gradient = m_axis;
			}
		}

		return finalDist;
	}

	Vector3 Cylinder::getVelocity(const Vector3& p) const
	{
		return m_linearVelocity + m_angularVelocity.cross(p - m_center);
	}

	void Cylinder::setTransform(const Vector3& t, const Eigen::Quaternion<Real>& r, const Vector3& s, Real dt)
	{
		Vector3 old_center = m_center;
		Eigen::Quaternion<Real> old_rotation = m_rotation;
		m_center = t + (r * m_local_offset);
		m_rotation = r;

		m_axis = r * m_axis_init;
		m_axis.normalize();

		if (dt > 1e-9)
		{
			m_linearVelocity = (m_center - old_center) / dt;

			Eigen::Quaternion<Real> q_diff = m_rotation * old_rotation.inverse();

			if (q_diff.w() < 0)
			{
				q_diff.w() = -q_diff.w();
				q_diff.vec() = -q_diff.vec();
			}

			m_angularVelocity = (q_diff.vec() * 2.0) / dt;

			m_hasExplicitVelocity = true;
		}
		else
		{
			m_linearVelocity.setZero();
			m_angularVelocity.setZero();
		}

		updateAABB();
	}

	void Cylinder::updateFrame()
	{
		m_center_prev = m_center;
		m_axis_prev = m_axis;
	}

	void Cylinder::updateAABB()
	{

		Real halfH = m_height * 0.5;
		Vector3 absAxis = m_axis.cwiseAbs();

		Real ex = halfH * absAxis.x() + m_radius * std::sqrt(1.0 - m_axis.x() * m_axis.x());
		Real ey = halfH * absAxis.y() + m_radius * std::sqrt(1.0 - m_axis.y() * m_axis.y());
		Real ez = halfH * absAxis.z() + m_radius * std::sqrt(1.0 - m_axis.z() * m_axis.z());

		Real pad = m_threshold;

		x_min = m_center.x() - ex - pad;
		x_max = m_center.x() + ex + pad;
		y_min = m_center.y() - ey - pad;
		y_max = m_center.y() + ey + pad;
		z_min = m_center.z() - ez - pad;
		z_max = m_center.z() + ez + pad;
	}

	bool Cylinder::isInsideAABB(Vector3 pos) const
	{
		if (pos.x() < x_min || pos.x() > x_max)
			return false;
		if (pos.y() < y_min || pos.y() > y_max)
			return false;
		if (pos.z() < z_min || pos.z() > z_max)
			return false;
		return true;
	}

	std::tuple<Eigen::Matrix<Real, Eigen::Dynamic, 3>, std::vector<std::vector<trimesh::index_t>>> Cylinder::getMesh() const
	{
		int segments = 32;
		Real hh = m_height * 0.5;

		int numVertices = 2 + (2 * segments);
		Eigen::Matrix<Real, Eigen::Dynamic, 3> V(numVertices, 3);

		auto transform = [&](const Vector3& local) -> Vector3
		{
			return m_center + m_rotation * local;
		};

		V.row(0) = transform(Vector3(0, 0, hh)); // Top Center (Index 0)
		V.row(1) = transform(Vector3(0, 0, -hh)); // Bot Center (Index 1)

		// Rings
		for (int i = 0; i < segments; ++i)
		{
			Real theta = 2.0 * M_PI * (Real)i / (Real)segments;
			Real c = std::cos(theta);
			Real s = std::sin(theta);

			V.row(2 + i) = transform(Vector3(m_radius * c, m_radius * s, hh));

			V.row(2 + segments + i) = transform(Vector3(m_radius * c, m_radius * s, -hh));
		}

		std::vector<std::vector<trimesh::index_t>> F;
		F.reserve(segments * 4); // 2 caps + 2 tris per side quad

		for (int i = 0; i < segments; ++i)
		{
			trimesh::index_t next = (i + 1) % segments;

			trimesh::index_t topCenter = 0;
			trimesh::index_t botCenter = 1;
			trimesh::index_t t0 = 2 + i; // Top Current
			trimesh::index_t t1 = 2 + next; // Top Next
			trimesh::index_t b0 = 2 + segments + i; // Bot Current
			trimesh::index_t b1 = 2 + segments + next; // Bot Next

			F.push_back({ topCenter, t0, t1 });

			F.push_back({ botCenter, b1, b0 });

			F.push_back({ t0, b0, b1 });
			F.push_back({ t0, b1, t1 });
		}

		return { V, F };
	}

	ObstacleDiscregrid::FrameData ObstacleDiscregrid::loadFrameData(const std::string& basePath, int frameIdx, bool loadObj)
	{
		FrameData data;
		data.frameIdx = frameIdx;

		std::string frameStr = std::to_string(frameIdx);
		std::string sdfPath = basePath + frameStr + ".dg";

		std::ifstream f(sdfPath);
		if (f.good())
		{
			data.sdf = std::make_shared<Discregrid::CubicLagrangeDiscreteGrid>(sdfPath);
		}
		else
		{
			std::cerr << "[SDF] Error: File not found " << sdfPath << std::endl;
			return data;
		}

		if (loadObj)
		{
			std::string objPath = basePath + frameStr + ".obj";
			try
			{
				Matrix2N uvs;
				std::tie(data.faces, data.vertices, uvs) = util::loadMeshFromObj(objPath);
				std::cout << "[SMPL] Loaded frame: " << frameIdx << std::endl;
			}
			catch (...)
			{
				std::cout << "[SDF] Warning: OBJ load failed or file missing: " << objPath << std::endl;
				data.faces.clear();
				data.vertices.resize(0, 3);
			}
		}
		else
		{
			data.faces.clear();
			data.vertices.resize(0, 3);
		}

		data.isValid = true;
		return data;
	}

	void ObstacleDiscregrid::build(nlohmann::json& data)
	{
		m_name = data.value("name", "Unnamed_SDF");
		m_threshold = data.value("threshold", 0.01);
		m_Kpenalty = data.value("Kpenalty", 1000.0);
		m_Kfriction = data.value("Kfriction", 0.1);
		m_total_frames = data.value("total_frames", 1);

		m_loadMesh = data.value("load_mesh", false);

		std::string relPath = data["path"].get<std::string>();
		m_basePath = argus::project_root + "/" + relPath;

		std::cout << "[SDF] Initializing " << m_name << " (Frame 0)..." << std::endl;

		m_currentData = loadFrameData(m_basePath, 0, m_loadMesh);

		if (!m_currentData.isValid || !m_currentData.sdf)
		{
			throw std::runtime_error("CRITICAL: Failed to load initial SDF frame 0 for " + m_name);
		}
	}

	void ObstacleDiscregrid::startAsyncLoad(int frameIdx)
	{
		if (m_isFrozen)
			return;
		if (frameIdx >= m_total_frames)
			return;
		if (m_currentData.frameIdx == frameIdx)
			return;

		m_loadingFuture = std::async(std::launch::async,
		    &ObstacleDiscregrid::loadFrameData,
		    m_basePath,
		    frameIdx,
		    m_loadMesh);
	}

	void ObstacleDiscregrid::updateToActive()
	{
		if (m_isFrozen)
			return;

		if (m_loadingFuture.valid())
		{
			FrameData nextData = m_loadingFuture.get();

			if (nextData.isValid && nextData.sdf)
			{
				m_currentData = std::move(nextData);
			}
			else
			{
				std::cerr << "[SDF] Warning: Failed to swap to next frame. Keeping Frame " << m_currentData.frameIdx << std::endl;
			}
		}
	}

	std::tuple<Eigen::Matrix<Real, Eigen::Dynamic, 3>, std::vector<std::vector<trimesh::index_t>>> ObstacleDiscregrid::getMesh() const
	{
		if (!m_currentData.isValid)
			return {};

		std::vector<std::vector<trimesh::index_t>> polyFaces;
		polyFaces.reserve(m_currentData.faces.size());

		for (const auto& t : m_currentData.faces)
		{
			polyFaces.push_back({ t.v[0], t.v[1], t.v[2] });
		}

		return { m_currentData.vertices.transpose(), polyFaces };
	}
	Real ObstacleDiscregrid::getSDF(const Vector3& p, Vector3* gradient) const
	{
		if (!m_currentData.sdf)
		{
			if (gradient)
				gradient->setZero();
			return 0;
		}

		if (m_currentData.sdf->domain().contains(p))
		{
			Real dist = (Real)m_currentData.sdf->interpolate(0, p, gradient ? gradient : NULL);
			if (gradient)
			{

				if (gradient->squaredNorm() < 1e-12)
					*gradient = Vector3::Zero();
				else
					gradient->normalize();
			}
			return dist;
		}
		else
		{
			if (gradient)
				gradient->setZero();
			return 0;
		}
	}

	ObstacleSSDF::FrameData ObstacleSSDF::loadFrameData(const std::string& basePath, int frameIdx, bool loadObj)
	{
		FrameData data;
		// frameIdx += 200;
		data.frameIdx = frameIdx;

		std::string frameStr = std::to_string(frameIdx);
		std::string ssdfPath = basePath + frameStr + ".ssdf";

		std::ifstream f(ssdfPath);
		if (f.good())
		{
			f.close();
			data.volume = std::make_shared<SDFVolume>(ssdfPath);
		}
		else
		{
			std::cerr << "[SSDF] Error: File not found " << ssdfPath << std::endl;
			return data;
		}

		if (loadObj)
		{
			std::string objPath = basePath + frameStr + ".obj";
			try
			{
				Matrix2N uvs;
				std::tie(data.faces, data.vertices, uvs) = util::loadMeshFromObj(objPath);
			}
			catch (...)
			{
				data.faces.clear();
				data.vertices.resize(0, 3);
			}
		}
		else
		{
			data.faces.clear();
			data.vertices.resize(0, 3);
		}

		data.isValid = true;
		return data;
	}

	void ObstacleSSDF::build(nlohmann::json& data)
	{
		m_name = data.value("name", "Unnamed_SSDF");
		m_threshold = data.value("threshold", 0.01);
		m_Kfriction = data.value("Kfriction", 0.1);
		m_total_frames = data.value("total_frames", 1);
		m_loadMesh = data.value("load_mesh", true);

		std::string relPath = data["path"].get<std::string>();
		m_basePath = /* argus::project_root + "/" + */ relPath;

		std::cout << "[SSDF] Initializing " << m_name << " (Frame 0)..." << std::endl;

		m_currentData = loadFrameData(m_basePath, 0, m_loadMesh);

		if (!m_currentData.isValid || !m_currentData.volume)
		{
			throw std::runtime_error("CRITICAL: Failed to load initial SSDF frame 0 for " + m_name);
		}
	}

	void ObstacleSSDF::startAsyncLoad(int frameIdx)
	{
		if (m_isFrozen)
			return;
		if (frameIdx >= m_total_frames)
			return;
		if (m_currentData.frameIdx == frameIdx)
			return;

		m_loadingFuture = std::async(std::launch::async,
		    &ObstacleSSDF::loadFrameData,
		    m_basePath,
		    frameIdx,
		    m_loadMesh);
	}

	void ObstacleSSDF::updateToActive()
	{
		if (m_isFrozen)
			return;

		if (m_loadingFuture.valid())
		{
			FrameData nextData = m_loadingFuture.get();
			if (nextData.isValid && nextData.volume)
			{
				m_currentData = std::move(nextData);
			}
			else
			{
				std::cerr << "[SSDF] Warning: Failed to swap to next frame. Keeping Frame " << m_currentData.frameIdx << std::endl;
			}
		}
	}

	Real ObstacleSSDF::getSDF(const Vector3& p, Vector3* gradient) const
	{
		if (!m_currentData.volume)
		{
			if (gradient)
				gradient->setZero();
			return 0;
		}
		Eigen::Vector3d p_d = p.cast<double>();

		double dist = m_currentData.volume->getDistance(p_d);

		if (gradient)
		{
			*gradient = computeFiniteDifferenceNormal(p);
		}

		return (Real)dist;
	}

	Vector3 ObstacleSSDF::getVelocity(const Vector3& p) const
	{
		if (!m_currentData.volume)
			return Vector3::Zero();

		Eigen::Vector3f v = m_currentData.volume->getVelocity(p.cast<double>());

		return v.cast<Real>();
	}

	Vector3 ObstacleSSDF::computeFiniteDifferenceNormal(const Vector3& p) const
	{
		const double h = 1e-4;
		Eigen::Vector3d p_d = p.cast<double>();

		double dx = m_currentData.volume->getDistance(p_d + Eigen::Vector3d(h, 0, 0)) - m_currentData.volume->getDistance(p_d - Eigen::Vector3d(h, 0, 0));

		double dy = m_currentData.volume->getDistance(p_d + Eigen::Vector3d(0, h, 0)) - m_currentData.volume->getDistance(p_d - Eigen::Vector3d(0, h, 0));

		double dz = m_currentData.volume->getDistance(p_d + Eigen::Vector3d(0, 0, h)) - m_currentData.volume->getDistance(p_d - Eigen::Vector3d(0, 0, h));

		Vector3 normal((Real)dx, (Real)dy, (Real)dz);

		if (normal.squaredNorm() > 1e-12)
			normal.normalize();
		else
			normal = Vector3::UnitX();

		return normal;
	}

	std::tuple<Eigen::Matrix<Real, Eigen::Dynamic, 3>, std::vector<std::vector<trimesh::index_t>>> ObstacleSSDF::getMesh() const
	{
		if (!m_currentData.isValid)
			return {};

		std::vector<std::vector<trimesh::index_t>> polyFaces;
		polyFaces.reserve(m_currentData.faces.size());

		for (const auto& t : m_currentData.faces)
		{
			polyFaces.push_back({ t.v[0], t.v[1], t.v[2] });
		}

		return { m_currentData.vertices.transpose(), polyFaces };
	}

	void ObstacleSSDF::exportToVTK(const std::string& path) const
	{
		if (m_currentData.volume)
		{
			m_currentData.volume->exportToVTK(path);
		}
	}
	void ObstacleSMPLSDF::build(nlohmann::json& data)
	{
		// --- 1. Load rest SDF ---
		std::string sdfPath = data["sdf_path"];
		m_restSDF = std::make_shared<SDFVolume>(sdfPath);

		if (!m_restSDF)
		{
			std::cerr << "Failed to load rest SDF: " << sdfPath << std::endl;
		}

		// --- 2. Load rest vertices ---
		std::string restPath = data["rest_vertices_path"];
		std::ifstream vin(restPath, std::ios::binary);

		if (!vin.is_open())
		{
			std::cerr << "Failed to open rest vertices: " << restPath << std::endl;
			return;
		}

		vin.read((char*)&m_numVertices, sizeof(int));

		// m_restVertices.resize(m_numVertices, 3);
		// vin.read((char*)m_restVertices.data(), m_numVertices * 3 * sizeof(float));

		m_restVertices.resize(m_numVertices, 3);

		std::vector<float> buffer(m_numVertices * 3);
		vin.read((char*)buffer.data(), buffer.size() * sizeof(float));

		for (int i = 0; i < m_numVertices; i++)
		{
			m_restVertices(i, 0) = (Real)buffer[3 * i + 0];
			m_restVertices(i, 1) = (Real)buffer[3 * i + 1];
			m_restVertices(i, 2) = (Real)buffer[3 * i + 2];
		}
		vin.close();

		std::cout
		    << "[SMPL] Loaded rest vertices: " << m_numVertices << std::endl;

		// --- 3. Load weights ---
		std::string weightPath = data["weights_path"];
		std::ifstream win(weightPath, std::ios::binary);

		if (!win.is_open())
		{
			std::cerr << "Failed to open weights: " << weightPath << std::endl;
			return;
		}

		win.read((char*)&m_numJoints, sizeof(int));

		m_weights.resize(m_numVertices * m_numJoints);
		win.read((char*)m_weights.data(), m_weights.size() * sizeof(float));

		win.close();

		std::cout << "[SMPL] Loaded weights: " << m_numJoints << " joints" << std::endl;
		m_basePath = data["transforms_path"];
		m_total_frames = data.value("total_frames", 1);

		// load first frame immediately
		m_currentData = loadFrameData(m_basePath, 0);

		if (m_currentData.isValid)
		{
			updateTransforms(m_currentData.transforms);
		}
		std::string facesPath = data["faces_path"];
		std::ifstream fin(facesPath, std::ios::binary);

		if (!fin.is_open())
		{
			std::cerr << "Failed to open faces: " << facesPath << std::endl;
			return;
		}

		int numFaces;
		fin.read((char*)&numFaces, sizeof(int));

		m_faces.resize(numFaces, 3);
		fin.read((char*)m_faces.data(), numFaces * 3 * sizeof(int));

		fin.close();

		std::cout << "[SMPL] Loaded faces: " << numFaces << std::endl;
		if (m_currentData.isValid)
		{
			updateTransforms(m_currentData.transforms);

			// store rest joints
			m_restJoints.resize(m_numJoints);
			for (int j = 0; j < m_numJoints; j++)
			{
				m_restJoints[j] = m_T[j].block<3, 1>(0, 3);
			}
		}
		m_dominantJoint.resize(m_numVertices);

		for (int i = 0; i < m_numVertices; i++)
		{
			const float* w = &m_weights[i * m_numJoints];

			int bestJoint = 0;
			float maxW = 0.0f;

			for (int j = 0; j < m_numJoints; j++)
			{
				if (w[j] > maxW)
				{
					maxW = w[j];
					bestJoint = j;
				}
			}

			m_dominantJoint[i] = bestJoint;
		}

		// m_cloud = std::make_unique<SMPLPointCloud>(m_restVertices);
		m_cloud = std::make_unique<SMPLPointCloud>(&m_restVertices);
		m_kdtree = std::make_unique<KDTree>(
		    3, *m_cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));

		m_kdtree->buildIndex();

		std::cout << "[SMPL] KD-tree built\n";

		if (!m_kdtree)
		{
			std::cerr << "KD-tree is NULL!" << std::endl;
			exit(0);
		}
	}
	bool ObstacleSMPLSDF::meshAvailable()
	{
		return (m_numVertices > 0 && m_faces.rows() > 0);
	}

	std::tuple<
	    Eigen::Matrix<Real, Eigen::Dynamic, 3>,
	    std::vector<std::vector<trimesh::index_t>>>
	ObstacleSMPLSDF::getMesh() const
	{
		int N = m_numVertices;

		Eigen::Matrix<Real, Eigen::Dynamic, 3> V(N, 3);

		for (int i = 0; i < N; i++)
		{
			Vector3 v_rest = m_restVertices.row(i).cast<Real>();

			Eigen::Matrix<Real, 4, 1> vh;
			vh.head<3>() = v_rest;
			vh[3] = (Real)1.0;

			Vector3 v_world = Vector3::Zero();

			const float* w = &m_weights[i * m_numJoints];

			for (int j = 0; j < m_numJoints; j++)
			{
				if (w[j] < 1e-5)
					continue; // skip tiny weights

				Vector3 vj = (m_T[j] * vh).head<3>();

				v_world += w[j] * vj;
			}

			V.row(i) = v_world;

			V.row(i) = v_world;
		}

		// --- faces ---
		std::vector<std::vector<trimesh::index_t>> F;
		F.reserve(m_faces.rows());

		for (int i = 0; i < m_faces.rows(); i++)
		{
			std::vector<trimesh::index_t> face(3);
			face[0] = m_faces(i, 0);
			face[1] = m_faces(i, 1);
			face[2] = m_faces(i, 2);
			F.push_back(face);
		}

		return { V, F };
	}
	void ObstacleSMPLSDF::updateTransforms(const std::vector<Eigen::Matrix4d>& T)
	{
		m_T = T;

		int J = (int)T.size();

		if (J != m_numJoints)
		{
			std::cerr << "[SMPL] Transform size mismatch! Expected "
			          << m_numJoints << ", got " << J << std::endl;
			return;
		}

		m_T_inv.resize(J);

		for (int j = 0; j < J; j++)
		{
			m_T_inv[j] = m_T[j].inverse();
		}
		// std::cout << "[SMPL] root joint: "
		//           << m_T[0].block<3, 1>(0, 3).transpose()
		//           << std::endl;
	}
	int ObstacleSMPLSDF::findClosestVertex(const Vector3& p) const
	{
		if (!m_kdtree)
		{
			std::cerr << "[SMPL] KD-tree null!" << std::endl;
			return 0;
		}

		if (m_restVertices.rows() == 0)
		{
			std::cerr << "[SMPL] KD-tree empty!" << std::endl;
			return 0;
		}

		size_t ret_index = 0;
		double out_dist_sqr = 0.0;

		double query_pt[3] = { p[0], p[1], p[2] };

		nanoflann::KNNResultSet<double> resultSet(1);
		resultSet.init(&ret_index, &out_dist_sqr);

		nanoflann::SearchParameters params;
		params.sorted = true;

		m_kdtree->findNeighbors(resultSet, query_pt, params);

		if (ret_index >= (size_t)m_numVertices)
		{
			std::cerr << "[SMPL] KD returned invalid index: "
			          << ret_index << std::endl;
			return 0;
		}

		return (int)ret_index;
	}
	Vector3 ObstacleSMPLSDF::inverseSkinning(const Vector3& p_world) const
	{
		// --- STEP 1: initial guess using root joint ---
		// (simple, stable, works well in practice)
		if (m_T_inv.empty())
		{
			std::cerr << "[SMPL] ERROR: m_T_inv is empty!" << std::endl;
			return p_world;
		}
		int root = 0;

		Eigen::Matrix<Real, 4, 1> ph;
		ph.head<3>() = p_world;
		ph[3] = (Real)1.0;
		Vector3 p_rest = (m_T_inv[root] * ph).head<3>();
		if (!p_rest.allFinite())
		{
			std::cerr << "[SMPL] p_rest is invalid: " << p_rest.transpose() << std::endl;
			exit(0);
			return p_world;
		}
		// --- STEP 2: KD-tree lookup in correct space (rest space) ---
		int vid = findClosestVertex(p_rest);
		if (vid < 0 || vid >= m_numVertices)
		{
			std::cerr << "[SMPL] Invalid vid: " << vid << std::endl;
			return p_world;
		}
		// --- STEP 3: refine using dominant joint of closest vertex ---
		int bestJoint = m_dominantJoint[vid];

		p_rest = (m_T_inv[bestJoint] * ph).head<3>();

		return p_rest;
	}
	Real ObstacleSMPLSDF::getSDF(const Vector3& p, Vector3* gradient) const
	{
		if (!m_restSDF)
		{
			if (gradient)
				gradient->setZero();
			return 0;
		}

		// --- 1. map to rest pose ---
		Vector3 p_rest = inverseSkinning(p);

		// --- 2. query SDF ---
		Eigen::Vector3d pd = p_rest.cast<double>();

		float dist_f;
		Eigen::Vector3f vel;

		m_restSDF->sample(pd, dist_f, vel);

		Real dist = (Real)dist_f;

		// --- 3. gradient (simple + robust) ---
		if (gradient)
		{
			const Real h = 1e-4;

			Vector3 gx = inverseSkinning(p + Vector3(h, 0, 0));
			Vector3 gx2 = inverseSkinning(p - Vector3(h, 0, 0));

			Vector3 gy = inverseSkinning(p + Vector3(0, h, 0));
			Vector3 gy2 = inverseSkinning(p - Vector3(0, h, 0));

			Vector3 gz = inverseSkinning(p + Vector3(0, 0, h));
			Vector3 gz2 = inverseSkinning(p - Vector3(0, 0, h));

			Real dx = m_restSDF->getDistance(gx.cast<double>()) - m_restSDF->getDistance(gx2.cast<double>());
			Real dy = m_restSDF->getDistance(gy.cast<double>()) - m_restSDF->getDistance(gy2.cast<double>());
			Real dz = m_restSDF->getDistance(gz.cast<double>()) - m_restSDF->getDistance(gz2.cast<double>());

			Vector3 g(dx, dy, dz);

			if (g.squaredNorm() > 1e-12)
				g.normalize();
			else
				g = Vector3::UnitX();

			*gradient = g;
		}

		return dist;
	}
	bool ObstacleSMPLSDF::inBounds(const Vector3& p, Real threshold) const
	{
		return true;
	}
	Vector3 ObstacleSMPLSDF::getVelocity(const Vector3& p) const
	{
		return Vector3::Zero();
	}
	ObstacleSMPLSDF::FrameData
	ObstacleSMPLSDF::loadFrameData(const std::string& basePath, int frameIdx)
	{
		FrameData data;
		data.frameIdx = frameIdx;
		// std::stringstream ss;
		// ss << std::setw(5) << std::setfill('0') << frameIdx;

		// std::string file = basePath + "/transforms_" + ss.str() + ".bin";

		std::string file = basePath + std::to_string(frameIdx) + ".bin";
		// std::cout << "[SMPL] Loading frame " << frameIdx << std::endl;
		std::ifstream in(file, std::ios::binary);
		if (!in.is_open())
		{
			std::cerr << "[SMPL] Missing transforms: " << file << std::endl;
			return data;
		}

		int J = m_numJoints;
		data.transforms.resize(J);

		for (int j = 0; j < J; j++)
		{
			Eigen::Matrix4f Tf;
			in.read((char*)Tf.data(), 16 * sizeof(float));

			// data.transforms[j] = Tf.cast<double>();
			data.transforms[j] = Tf.transpose().cast<double>();
		}
		// std::cout << "[SMPL] root joint: "
		//           << data.transforms[0].block<3, 1>(0, 3).transpose()
		//           << std::endl;
		// std::cout << "[SMPL] Loaded frame " << frameIdx << std::endl;

		data.isValid = true;
		return data;
	}
	void ObstacleSMPLSDF::startAsyncLoad(int frameIdx)
	{
		if (frameIdx >= m_total_frames)
			return;

		// avoid reloading same frame
		if (m_currentData.frameIdx == frameIdx)
			return;

		// launch async load
		m_loadingFuture = std::async(
		    std::launch::async,
		    &ObstacleSMPLSDF::loadFrameData,
		    this,
		    m_basePath,
		    frameIdx);
	}
	void ObstacleSMPLSDF::updateToActive()
	{
		if (m_loadingFuture.valid())
		{
			FrameData next = m_loadingFuture.get();

			if (next.isValid)
			{
				m_currentData = std::move(next);

				// 🔥 THIS LINE IS THE KEY
				updateTransforms(m_currentData.transforms);

				// std::cout << "[SMPL] Applied frame "
				//           << m_currentData.frameIdx << std::endl;
			}
		}
	}
}