#ifndef COMMON_HPP
#define COMMON_HPP
#define EIGEN_RUNTIME_NO_MALLOC
#include <optional>
#include <Eigen/Eigen>
#include <halfedge/trimesh.h>
#ifdef BOGUS
#include <bogus/Core/Block.impl.hpp>
#endif
#include <TinyAD/Scalar.hh>
#include <nlohmann/json.hpp>
#define _USE_MATH_DEFINES
#include <math.h>
#include "../../vendor/SecondFundamentalForm/SecondFundamentalFormDiscretization.h"
#include "../../vendor/SecondFundamentalForm/MidedgeAverageFormulation.h"

static const double PI = 3.1415926535898;
// useful
namespace nlohmann
{
	template <typename Scalar, int Rows, int Options>
	struct adl_serializer<Eigen::Matrix<Scalar, Rows, 1, Options>>
	{
		// Function to convert an Eigen Vector -> JSON
		static void to_json(json& j, const Eigen::Matrix<Scalar, Rows, 1, Options>& vector)
		{
			// This is the same logic as before
			for (int i = 0; i < vector.size(); ++i)
			{
				j.push_back(vector(i));
			}
		}

		// Function to convert JSON -> Eigen Vector
		static void from_json(const json& j, Eigen::Matrix<Scalar, Rows, 1, Options>& vector)
		{
			// Same logic as before
			const auto values = j.get<std::vector<Scalar>>();
			vector = Eigen::Map<const Eigen::Matrix<Scalar, Rows, 1>>(values.data(), values.size());
		}
	};
}
namespace argus
{
#ifdef ARGUS_USE_FLOATS
	using Real = float;
	using AReal15 = TinyAD::Float<15>;
	using AReal1 = TinyAD::Float<1>;
	using AReal2 = TinyAD::Float<2>;
	using AReal3 = TinyAD::Float<3>;
	using AReal4 = TinyAD::Float<4>;
	using AReal8 = TinyAD::Float<8>;
#else
	using Real = double;
	using AReal15 = TinyAD::Double<15>;
	using AReal1 = TinyAD::Double<1>;
	using AReal2 = TinyAD::Double<2>;
	using AReal3 = TinyAD::Double<3>;
	using AReal4 = TinyAD::Double<4>;
	using AReal8 = TinyAD::Double<7>;

#endif
	// project root
	const std::string project_root = PROJECT_ROOT_DIR;

	extern nlohmann::json simConf;
	extern std::shared_ptr<SecondFundamentalFormDiscretization> sffobject;

	template <typename T>
	inline auto toPassive(T x) { return TinyAD::to_passive(x); }

	// TODO: Check for Row major vs Col major performance
	// using MatrixN3 = Eigen::Matrix<Real, Eigen::Dynamic, 3>;
	using Matrix3N = Eigen::Matrix<Real, 3, Eigen::Dynamic>;
	using Matrix2N = Eigen::Matrix<Real, 2, Eigen::Dynamic>;
	using MatrixNN = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;
	using Matrix22 = Eigen::Matrix<Real, 2, 2>;
	using Matrix23 = Eigen::Matrix<Real, 2, 3>;
	using Matrix33 = Eigen::Matrix<Real, 3, 3>;
	using Matrix44 = Eigen::Matrix<Real, 4, 4>;
	using Matrix32 = Eigen::Matrix<Real, 3, 2>;
	using DiagonalMatrix33 = Eigen::DiagonalMatrix<Real, 3>;
	using VectorN = Eigen::Vector<Real, Eigen::Dynamic>;
	using RowVector3 = Eigen::Matrix<Real, 1, 3>;
	using Vector3 = Eigen::Vector<Real, 3>;
	using Vector4 = Eigen::Vector<Real, 4>;
	using Vector2 = Eigen::Vector<Real, 2>;
	using SparseMatrix = Eigen::SparseMatrix<Real>;
	using SparseTriplet = Eigen::Triplet<Real>;
	using VectorNi = Eigen::Vector<int, Eigen::Dynamic>;
	// Using row major as these occur on LHS
	// profiling results also show that row major is faster
	// We should also inspect other matrix types for performance
	// Refer http://gdaviet.fr/doc/bogus/master/doxygen/block.html
	// this library has many caveats that need to be kept in mind
	// while using which are hard to debug
	// using BlockSparseMatrix = bogus::SparseBlockMatrix<Matrix33>;
	// TODO: decide if SYMMETRIC flag should be used, we can reduce
	// memory usage by half but lose performance due to worse cache locality
	// using SymmetricBlockSparseMatrix = bogus::SparseBlockMatrix<Matrix33, bogus::flags::SYMMETRIC>;
	// for computing second fundamental form
	struct ClusterInfo
	{
		std::vector<std::vector<std::pair<int, int>>> clusters;
		std::vector<int> jump;
	};

	template <class T>
	using Map = Eigen::Map<T>;

	template <class T>
	inline Map<VectorN> FlattenMatrixN3(T& vector)
	{
		return Map<VectorN>(vector.data(), 3 * vector.cols());
	}
	template <class T>
	inline Map<const VectorN> FlattenMatrixN3(const T& vector)
	{
		return Map<const VectorN>(vector.data(), 3 * vector.cols());
	}
	template <class T>
	inline Map<const VectorN> FlattenMatrixN2(const T& vector)
	{
		return Map<const VectorN>(vector.data(), 2 * vector.cols());
	}

	template <class T>
	inline Map<Matrix3N> GroupMatrixN3(T& vector)
	{
		return Map<Matrix3N>(vector.data(), 3, vector.size() / 3);
	}
	template <class T>
	inline Map<const Matrix3N> GroupMatrixN3(const T& vector)
	{
		return Map<const Matrix3N>(vector.data(), 3, vector.size() / 3);
	}

	inline Matrix33 ConvertToSkew(const Vector3& vec)
	{
		Matrix33 skew;
		skew << 0, -vec(2), vec(1),
		    vec(2), 0, -vec(0),
		    -vec(1), vec(0), 0;
		return skew;
	}

	inline void set3x3(SparseMatrix& sparseMat, Matrix33& other, size_t i, size_t j)
	{
		for (int x = 0; x < 3; x++)
		{
			for (int y = 0; y < 3; y++)
			{
				sparseMat.coeffRef(x + i, y + j) = other(x, y);
			}
		}
	}

	inline void add3x3(SparseMatrix& sparseMat, Matrix33& other, const trimesh::index_t i, const trimesh::index_t j)
	{
		for (int x = 0; x < 3; x++)
		{
			for (int y = 0; y < 3; y++)
			{
				sparseMat.coeffRef(x + i, y + j) += other(x, y);
			}
		}
	}

	inline void sub3x3(SparseMatrix& sparseMat, Matrix33& other, const trimesh::index_t i, const trimesh::index_t j)
	{
		for (int x = 0; x < 3; x++)
		{
			for (int y = 0; y < 3; y++)
			{
				sparseMat.coeffRef(x + i, y + j) -= other(x, y);
			}
		}
	}

	inline Real delta(trimesh::index_t i, trimesh::index_t j)
	{
		return i == j ? (Real)1 : (Real)0;
	}

	// https://stackoverflow.com/a/46111789
	template <typename XprType>
	auto DuplicateVectorElems(const XprType& xpr, Eigen::Index K)
	{
		return xpr(Eigen::ArrayXi::LinSpaced(xpr.size() * K, 0, xpr.size() - 1));
	}

	// TODO: https://stackoverflow.com/a/12787599
	// No need to "copy" data as the buffer acquires
	// meaning after call to getVectorN
	// Edit: I realised that what I want was a raw T array...
	class SharedVector
	{
		std::vector<Real> m_Buffer;

	public:
		Map<VectorN> getVectorN(size_t n)
		{
			if (n > m_Buffer.size())
			{
				m_Buffer.resize(n);
			}
			return Map<VectorN>(m_Buffer.data(), n);
		};
		void setCapacity(size_t n) { m_Buffer.resize(n); }
	};

	typedef struct
	{
		Real m_SpringConstant;
		Real m_DampingConstant;
	} SpringProperties;
	struct QuadraturePoints
	{
		Real u;
		Real v;
		Real weight;
	};
	typedef struct
	{
		Real weights[25][3];
	} probeWeights;
	// 1. Define the point structure ONCE.
	// Aligning to 128 bytes ensures it occupies exactly two cache lines.
	struct alignas(128) OptimizedPoint
	{
		Vector3 pos; // 24 bytes (Current position)
		Real radius; // 8 bytes  (Precomputed rest-state radius)
		Real b[3]; // 24 bytes (Barycentric weights)
		int vIdx[3]; // 12 bytes (Vertex indices)
		int faceId; // 4 bytes  (Parent triangle ID)
		Real wInv;
		// Total used: 72 bytes.
		// Padding ensures the next point in the vector starts at a 128-byte boundary.
		Real padding[7]; // 56 bytes
	};
	struct PrecomputedSelfCollisionData
	{
		static constexpr int K = 15; // Number of active probes

		// Precomputed weights (constant for all faces)
		alignas(64) Real flattenedWeights[K * 3];

		// Precomputed radii (unique per face/probe)
		std::vector<Real> flattenedRadii;

		// The runtime point cloud scratchpad
		// This will be resized to (numFaces * K)
		std::vector<OptimizedPoint> pointCloudBuffer;

		// --- Added parameters ---
		Real gridCellSize; // This is your 'h'
		Real maxRadius; // Useful for debugging or relative scaling
		Real velocity_damp_factor;
	};
	struct UnionFind
	{
		std::vector<int> parent;
		std::vector<int> sign;
		UnionFind(int items)
		{
			parent.resize(items);
			sign.resize(items);
			for (int i = 0; i < items; i++)
			{
				parent[i] = i;
				sign[i] = 1;
			}
		}

		std::pair<int, int> find(int i)
		{
			if (parent[i] != i)
			{
				std::pair<int, int> newparent = find(parent[i]);
				sign[i] *= newparent.second;
				parent[i] = newparent.first;
			}

			return { parent[i], sign[i] };
		}

		void dounion(int i, int j, int usign)
		{
			std::pair<int, int> xroot = find(i);
			std::pair<int, int> yroot = find(j);
			if (xroot.first != yroot.first)
			{
				parent[xroot.first] = yroot.first;
				sign[xroot.first] = usign * xroot.second * yroot.second;
			}
		}
	};
	// info for vbd
	struct SheetInfo
	{
		std::string sheetPath = "";
		std::string seam_path = "";
		bool wrinkle_mesh = false;
		std::vector<int> clampedVertices = {};
		std::map<int, std::unordered_map<int, Vector3>> clampedVerticesMotion;
		Eigen::MatrixXi vertexFacesSet = Eigen::MatrixXi::Zero(0, 0);
		Eigen::MatrixXi coloredSets = Eigen::MatrixXi::Zero(0, 0);
		Eigen::VectorXi n_adjacentFaces = Eigen::VectorXi::Zero(0);
		VectorN axt = VectorN::Zero(0);
		VectorN aat = VectorN::Zero(0);
		int ndofs = 0;
	};

	struct TimestepperInfo
	{
		std::string ampsposupdate = "separate"; // default
		int obstacle_collision_maxiter = 1;
		int self_collisions_maxsubsteps = 1;
		bool self_collisions = false;
		bool self_collisions_ccd = false;
		Real self_collision_stiffness = 0.001;
		Real self_collision_damping = 0.0;
		Real self_collision_threshold = 0.001;
		std::string sdf_collider = "basic";
		Real kdamp_amps = 0;
		Real kdamp_pos = 0;
	};
}

#endif /* COMMON_HPP */
