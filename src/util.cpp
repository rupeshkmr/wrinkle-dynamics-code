#include <argus/util.hpp>

#include <fstream>

#include <boost/filesystem.hpp>
#include <stdexcept>

namespace argus
{
	namespace util
	{
		/* 	VectorN convertSTDVecVecN(const std::vector<argus::Real>& ip)
		    {
		        VectorN op;
		        op.resize(ip.size());
		        for (int i=0; i<ip.size(); i++)
		            op(i) = ip.at(i);
		        return op;
		    } */
		Eigen::SparseMatrix<int> buildAdjacencyMatrix(const Eigen::MatrixXi& F, int numVertices)
		{
			int numFaces = F.rows();

			// 1. Build Vertex-to-Face Adjacency List
			// logical map: vertex_index -> [face_index_1, face_index_2, ...]
			std::vector<std::vector<int>> vertToFaces(numVertices);
			for (int i = 0; i < numFaces; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					int vIdx = F(i, j);
					vertToFaces[vIdx].push_back(i);
				}
			}

			// 2. Build Face-to-Face Adjacency Triplets
			std::vector<Eigen::Triplet<int>> triplets;
			// Reserve memory estimation (avg valence ~6, so 3 verts * 6 neighbors = 18)
			triplets.reserve(numFaces * 18);

			for (int f = 0; f < numFaces; ++f)
			{
				// A face is always adjacent to itself
				triplets.emplace_back(f, f, 1);

				// Check all 3 vertices of this face
				for (int j = 0; j < 3; ++j)
				{
					int vIdx = F(f, j);
					// Look at all OTHER faces connected to this vertex
					const auto& neighbors = vertToFaces[vIdx];
					for (int neighborFace : neighbors)
					{
						if (neighborFace != f)
						{
							triplets.emplace_back(f, neighborFace, 1);
						}
					}
				}
			}

			// 3. Construct Matrix
			Eigen::SparseMatrix<int> adjMat(numFaces, numFaces);
			adjMat.setFromTriplets(triplets.begin(), triplets.end());
			return adjMat;
		}
		std::optional<Vector3> getVector3(const nlohmann::json& v)
		{
			/**
			 * Convert list data present in json to eigen vector
			 */
			Vector3 op;
			if (v.is_array() && v.size() == 3)
			{
				op(0) = v[0].get<double>();
				op(1) = v[1].get<double>();
				op(2) = v[2].get<double>();
				return op;
			}
			return std::nullopt;
		}

		void testScalarGradients(std::function<Real(VectorN&, VectorN&)> func)
		{
			VectorN X, dx, g, Xinit, g1, g2;
			Xinit = VectorN::Random() * 9.9;
			dx = VectorN::Random().normalized();
			Real h = 0.1;
			Real e = func(Xinit, g);
			for (int i = 0; i < 10; i++)
			{
				X = Xinit - dx * h / 2.0;
				Real s1 = func(X, g1);
				X = Xinit + dx * h / 2.0;
				Real s2 = func(X, g2);
				Real lhs = (s2 - s1) / h;
				Real rhs = dx.dot(g);
				std::cout << "Iter " << i << " Error " << fabs(lhs - rhs) << std::endl;
				h = h / 10.0;
			}
		}

		Matrix2N convertSTDVecMatrix2N(const std::vector<Vector2>& ip)
		{
			Matrix2N op;
			int n = ip.size();
			op.resize(2, n);
			for (int i = 0; i < n; i++)
			{
				op.col(i) = ip.at(i);
			}
			return op;
		}

		Matrix3N convertSTDVecMatrix3N(const std::vector<Vector3>& ip)
		{
			Matrix3N op;
			int n = ip.size();
			op.resize(3, n);
			for (int i = 0; i < n; i++)
			{
				op.col(i) = ip.at(i);
			}
			return op;
		}

		std::vector<argus::Vector2> convertMatrix2NToSTDVec(const Matrix2N& ip)
		{
			int n = ip.cols();
			std::vector<argus::Vector2> op;
			for (int i = 0; i < n; i++)
			{
				op.push_back(ip.col(i));
			}
			return op;
		}

		std::vector<argus::Vector3> convertMatrix3NToSTDVec(const Matrix3N& ip)
		{
			int n = ip.cols();
			std::vector<argus::Vector3> op;
			for (int i = 0; i < n; i++)
			{
				op.push_back(ip.col(i));
			}
			return op;
		}

		std::chrono::time_point<std::chrono::steady_clock> captureCurrentTime()
		{
			std::chrono::time_point<std::chrono::steady_clock> t = std::chrono::steady_clock::now();
			return t;
		}

		std::tuple<std::vector<trimesh::triangle_t>, Matrix3N, Matrix2N> generateTesselatedSquare(Vector3 center, size_t n, Real sideLength, Alignment alignment)
		{
			std::vector<trimesh::triangle_t> triangles;
			triangles.reserve(n * n * 2);
			trimesh::index_t verticesPerRow = n + 1;
			Matrix3N positions = center.replicate(1, verticesPerRow * verticesPerRow);
			Matrix2N uvCoords = Matrix2N::Zero(2, verticesPerRow * verticesPerRow);
			Real tileLength = sideLength / n;
			Real halfLength = sideLength / 2;
			for (trimesh::index_t i = 0; i < n; i++)
			{
				for (trimesh::index_t j = 0; j < n; j++)
				{
					trimesh::index_t prev = i * verticesPerRow;
					triangles.push_back({
					    j + prev,
					    j + prev + verticesPerRow,
					    j + prev + verticesPerRow + 1,
					});
					triangles.push_back({
					    j + prev + 1,
					    j + prev,
					    j + prev + verticesPerRow + 1,
					});
				}
			}
			for (size_t i = 0; i <= n; i++)
			{
				for (size_t j = 0; j <= n; j++)
				{
					uvCoords.col(i * verticesPerRow + j) = Vector2(i * sideLength, j * sideLength) / n;
					switch (alignment)
					{
					case argus::util::Alignment::Z:
						positions.col(i * verticesPerRow + j) += Vector3({ i * tileLength - halfLength,
						    j * tileLength - halfLength,
						    0 });
						break;
					case argus::util::Alignment::X:
						positions.col(i * verticesPerRow + j) += Vector3({ 0, j * tileLength - halfLength, i * tileLength - halfLength });
						break;
					case argus::util::Alignment::Y:
						positions.col(i * verticesPerRow + j) += Vector3({ j * tileLength - halfLength, 0, i * tileLength - halfLength });
						break;
					default:
						break;
					}
				}
			}
			return { triangles, positions, uvCoords };
		}

		std::tuple<std::vector<trimesh::triangle_t>, Matrix3N, Matrix2N> generateCylindricalSheet(Vector3 center, size_t n, Real radius, Real theta, Real height)
		{
			std::vector<trimesh::triangle_t> triangles;
			triangles.reserve(n * n * 2);
			trimesh::index_t verticesPerRow = n + 1;
			Matrix3N positions = center.replicate(1, verticesPerRow * verticesPerRow);
			Matrix2N uvCoords = Matrix2N::Zero(2, verticesPerRow * verticesPerRow);
			for (trimesh::index_t i = 0; i < n; i++)
			{
				for (trimesh::index_t j = 0; j < n; j++)
				{
					trimesh::index_t prev = i * verticesPerRow;
					triangles.push_back({
					    j + prev,
					    j + prev + verticesPerRow,
					    j + prev + verticesPerRow + 1,
					});
					triangles.push_back({
					    j + prev + 1,
					    j + prev,
					    j + prev + verticesPerRow + 1,
					});
				}
			}

			Real heightTileLength = height / n;
			Real heightHalfLength = height / 2;
			Real halfAngle = theta / 2;
			Real angleStep = theta / n;
			Real clothWidth = radius * theta;
			for (size_t i = 0; i <= n; i++)
			{
				Real angle = -halfAngle + i * angleStep;
				for (size_t j = 0; j <= n; j++)
				{
					uvCoords.col(i * verticesPerRow + j) = Vector2(j * clothWidth, i * height) / n;
					// uvCoords.col(i * verticesPerRow + j) = Vector2(i * clothWidth, j * height) / n;
					positions.col(i * verticesPerRow + j) += Vector3({ -radius * cos(angle),
					    j * heightTileLength - heightHalfLength,
					    radius * sin(angle) });
				}
			}
			return { triangles, positions, uvCoords };
		}

		// TODO returns incorrect uvs for meshes having non-constant y values.
		std::tuple<std::vector<trimesh::triangle_t>, Matrix3N, Matrix2N> loadMeshFromObj(std::string fileName)
		{
			std::vector<trimesh::triangle_t> triangles;
			Matrix3N positions;
			Matrix2N textureCoords;
			argus::MatrixNN V, TC, N;
			Eigen::MatrixXi F, FTC, FN;
			V.resize(0, 0);
			F.resize(0, 0);

			// Use: V (Vertices), F (Face Indices), TC (UV coords), FTC (UV Indices)
			if (igl::readOBJ(fileName, V, TC, N, F, FTC, FN))
			{
#ifdef ARGUS_DEBUG
				std::cout << "Read Success!\n";
#endif
			}
			else
			{
				std::cout << "Unable to read obj data !! Exiting\n";
				throw std::runtime_error("Error!");
			}

			positions = V.transpose();
			textureCoords.resize(2, positions.cols());
			textureCoords = TC.transpose();

			for (int i = 0; i < F.rows(); i++)
			{
				trimesh::index_t v0 = F(i, 0);
				trimesh::index_t v1 = F(i, 1);
				trimesh::index_t v2 = F(i, 2);
				triangles.push_back({ v0, v1, v2 });
			}

			return std::make_tuple(triangles, positions, textureCoords);
		}
		std::tuple<std::vector<trimesh::triangle_t>, Matrix3N, std::vector<trimesh::triangle_t>, Matrix2N> loadMeshWithUVFromObj(std::string fileName)
		{
			std::vector<trimesh::triangle_t> triangles, UVTris;
			Matrix3N positions;
			Matrix2N textureCoords;
			argus::MatrixNN V, TC, N;
			Eigen::MatrixXi F, FTC, FN;

			// Initialize
			V.resize(0, 0);
			F.resize(0, 0);
			TC.resize(0, 0);
			FTC.resize(0, 0);

			// Read OBJ
			if (igl::readOBJ(fileName, V, TC, N, F, FTC, FN))
			{
#ifdef ARGUS_DEBUG
				std::cout << "Read Success!\n";
#endif
			}
			else
			{
				std::cout << "Unable to read obj data !! Exiting\n";
				throw std::runtime_error("Error: Failed to read OBJ file.");
			}

			// --- FALLBACK LOGIC: Handle Missing UVs via PCA Projection ---
			if (TC.rows() == 0 || FTC.rows() == 0)
			{
				std::cout << "Warning: No UVs found. Generating UVs via PCA projection.\n";
				// 1. Center the vertices
				Eigen::RowVector3d centroid = V.colwise().mean();
				Eigen::MatrixXd V_centered = V.rowwise() - centroid;

				// 2. Compute Covariance / Scatter Matrix
				Eigen::Matrix3d covariance = V_centered.transpose() * V_centered;

				// 3. Eigen Decomposition (PCA)
				// SelfAdjointEigenSolver sorts eigenvalues in INCREASING order.
				// Index 0: Smallest variance (Normal direction of the sheet)
				// Index 1: Medium variance (Tangent / Bitangent)
				// Index 2: Largest variance (Tangent / Bitangent)
				Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(covariance);
				Eigen::Vector3d eigenvals = eigensolver.eigenvalues(); // Sorted: [small, medium, large]

				// Check "Flatness": Ratio of thickness (smallest var) to span (largest var)
				double thickness_ratio = eigenvals[0] / eigenvals[2];

				// If thickness is > 1% of the width, it's likely a 3D volume (cylinder/sphere), not a sheet.
				if (thickness_ratio > 0.01)
				{
					// throw std::runtime_error("Error: Mesh has no UVs and is not planar (significant curvature detected). Automatic flattening impossible.");
					std::cerr << "Error: Mesh has no UVs and is not planar (significant curvature detected). Automatic flattening impossible.\n";
				}
				else
				{
					Eigen::Matrix3d eigenvectors = eigensolver.eigenvectors();

					// The two dominant axes define the "Tangent Plane" of the mesh
					Eigen::Vector3d axis_u = eigenvectors.col(2); // Axis of most variance
					Eigen::Vector3d axis_v = eigenvectors.col(1); // Axis of second most variance

					// 4. Project vertices onto these axes to get UVs
					TC.resize(V.rows(), 2);
					TC.col(0) = V_centered * axis_u;
					TC.col(1) = V_centered * axis_v;

					// 5. Duplicate Topology
					// Since we generated 1 UV per Vertex, UV indices match Vertex indices
					FTC = F;
				}
			}
			// -------------------------------------------------------------

			// Transpose and Store Positions
			positions = V.transpose();

			// Resize strictly based on the loaded (or generated) TC count
			textureCoords.resize(2, TC.rows());
			textureCoords = TC.transpose();

			// Store Geometry Triangles
			triangles.reserve(F.rows());
			for (int i = 0; i < F.rows(); i++)
			{
				triangles.push_back({ (trimesh::index_t)F(i, 0), (trimesh::index_t)F(i, 1), (trimesh::index_t)F(i, 2) });
			}

			// Store UV Triangles
			UVTris.reserve(FTC.rows());
			for (int i = 0; i < FTC.rows(); ++i)
			{
				trimesh::triangle_t tri;
				tri.v[0] = FTC(i, 0);
				tri.v[1] = FTC(i, 1);
				tri.v[2] = FTC(i, 2);
				UVTris.push_back(tri);
			}

			return std::make_tuple(triangles, positions, UVTris, textureCoords);
		}

		std::tuple<std::vector<trimesh::triangle_t>, Matrix3N, Matrix2N> loadSheetFromObj(std::string fileName)
		{
			std::vector<trimesh::triangle_t> triangles;
			Matrix3N positions(3, 0);
			Matrix2N textureCoords(2, 0);
			std::ifstream objFile(fileName);
			std::string line;
			while (std::getline(objFile, line))
			{
				std::istringstream iss(line);
				std::string type;
				iss >> type;
				if (type == "v")
				{
					Real x, y, z;
					iss >> x >> y >> z;
					positions.conservativeResize(3, positions.cols() + 1);
					positions.col(positions.cols() - 1) << x, y, z;
				}
				else if (type == "vt")
				{
					Real u, v;
					iss >> u >> v;
					textureCoords.conservativeResize(2, textureCoords.cols() + 1);
					textureCoords.col(textureCoords.cols() - 1) << u, v;
				}
				else if (type == "f")
				{
					std::vector<std::string> faceVertexStrings;
					std::string faceVertexString;
					while (iss >> faceVertexString)
					{
						faceVertexStrings.push_back(faceVertexString);
					}
					trimesh::triangle_t triangle;
					for (int i = 0; i < 3; i++)
					{
						std::istringstream fvss(faceVertexStrings[i]);
						int positionIndex, textureCoordIndex, normalIndex;
						char delimiter;
						fvss >> positionIndex >> delimiter >> textureCoordIndex >> delimiter >> normalIndex;
						if (positionIndex != textureCoordIndex)
						{
							throw "argus needs positionIndex and textureCoordIndex to be same";
						}
						triangle.v[i] = positionIndex - 1; // Subtract 1 to convert from 1-based to 0-based indexing
					}
					triangles.push_back(triangle);
				}
			}
			return std::make_tuple(triangles, positions, textureCoords);
		}

		std::tuple<std::vector<trimesh::triangle_t>, Matrix3N> generateThinBeam(Vector3 center, Real length, Real width, Alignment alignment)
		{
			std::vector<trimesh::triangle_t> triangles = {
				{ 0, 4, 1 },
				{ 1, 4, 5 },
				{ 1, 5, 2 },
				{ 2, 5, 6 },
				{ 2, 6, 3 },
				{ 3, 6, 7 },
				{ 3, 7, 0 },
				{ 0, 7, 4 },
			};
			Matrix3N positions = center.replicate(1, 8);
			Real halfLength = length / 2;
			switch (alignment)
			{
			case Alignment::Z:
				positions.col(0) += Vector3({ width, width, halfLength });
				positions.col(1) += Vector3({ width, -width, halfLength });
				positions.col(2) += Vector3({ -width, -width, halfLength });
				positions.col(3) += Vector3({ -width, width, halfLength });
				positions.col(4) += Vector3({ width, width, -halfLength });
				positions.col(5) += Vector3({ width, -width, -halfLength });
				positions.col(6) += Vector3({ -width, -width, -halfLength });
				positions.col(7) += Vector3({ -width, width, -halfLength });
				break;
			case Alignment::X:
				positions.col(0) += Vector3({
				    halfLength,
				    width,
				    width,
				});
				positions.col(1) += Vector3({
				    halfLength,
				    width,
				    -width,
				});
				positions.col(2) += Vector3({
				    halfLength,
				    -width,
				    -width,
				});
				positions.col(3) += Vector3({
				    halfLength,
				    -width,
				    width,
				});
				positions.col(4) += Vector3({ -halfLength, width, width });
				positions.col(5) += Vector3({ -halfLength, width, -width });
				positions.col(6) += Vector3({ -halfLength, -width, -width });
				positions.col(7) += Vector3({ -halfLength, -width, width });
				break;
			case Alignment::Y:
				positions.col(0) += Vector3({
				    width,
				    halfLength,
				    width,
				});
				positions.col(1) += Vector3({
				    width,
				    halfLength,
				    -width,
				});
				positions.col(2) += Vector3({
				    -width,
				    halfLength,
				    -width,
				});
				positions.col(3) += Vector3({
				    -width,
				    halfLength,
				    width,
				});
				positions.col(4) += Vector3({ width, -halfLength, width });
				positions.col(5) += Vector3({ width, -halfLength, -width });
				positions.col(6) += Vector3({ -width, -halfLength, -width });
				positions.col(7) += Vector3({ -width, -halfLength, width });
				break;
			}
			return { triangles, positions };
		}
		// void saveClothData(std::string saveDir, int frame,const Matrix3N& positions, const VectorN& amplitudes, const VectorN& dphis) {
		// 	// taken from https://github.com/AleksandarHaber/Save-and-Load-Eigen-Cpp-Matrices-Arrays-to-and-from-CSV-files/blob/master/source_file.cpp
		// 	// make sure it is a directory
		// 	if(saveDir[saveDir.length()-1] != '/')
		// 		saveDir = saveDir + "/";

		// 	boost::filesystem::create_directories(saveDir);
		// 	const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");

		// 	std::string ampSavepath= saveDir + "amplitudes.csv." + std::to_string(frame);
		// 	std::string dphiSavepath= saveDir + "dphis.csv." + std::to_string(frame);
		// 	std::string positionsSavepath= saveDir + "positions.csv." + std::to_string(frame);

		// 	std::ofstream file;
		// 	// saving amplitudes
		// 	file.open(ampSavepath);
		// 	file << amplitudes.format(CSVFormat);
		// 	file.close();

		// 	// saving dphis
		// 	file.open(dphiSavepath);
		// 	file << dphis.format(CSVFormat);
		// 	file.close();

		// 	// saving positions
		// 	file.open(positionsSavepath);
		// 	file << positions.format(CSVFormat);
		// 	file.close();

		// }

		// void loadClothData(std::string loadDir, int frame, Matrix3N& positions, VectorN& amplitudes, VectorN& dphis) {
		// 	std::string ampSavepath= loadDir+ "/amplitudes.csv." + std::to_string(frame);
		// 	std::string dphiSavepath= loadDir+ "/dphis.csv." + std::to_string(frame);
		// 	std::string positionsSavepath= loadDir+ "/positions.csv." + std::to_string(frame);
		// 	// try
		// 	// {
		// 		positions = openData(positionsSavepath).reshaped(3, sheet.getPositionsN3().cols());
		// 		amplitudes = openData(ampSavepath);
		// 		dphis = openData(dphiSavepath);
		// 		std::cout << "Loaded frame " << frame << std::endl;
		// 	// }
		// 	// catch(const std::runtime_error& error)
		// 	// {
		// 	// 	std::cout << "Failed to load frame " << frame << std::endl;
		// 	// }

		// }

		Eigen::MatrixXi convertTrianglesToMatrix(const std::vector<trimesh::triangle_t>& triangles)
		{
			Eigen::MatrixXi mat(triangles.size(), 3);
			for (int i = 0; i < triangles.size(); i++)
			{
				mat(i, 0) = triangles[i].v[0];
				mat(i, 1) = triangles[i].v[1];
				mat(i, 2) = triangles[i].v[2];
			}
			return mat;
		}

		/* void saveCheckpoint(std::vector<std::string> data_to_save, std::string parent_directory, int frame, argus::Sheet& sheet)
		{
		    // Only save checkpoints if required
		    if (sheet.saveFrames())
		    {
		        std::map<std::string, int>& m = sheet.getParamsToSave();
		        m.clear();
		        for (auto it = data_to_save.begin(); it != data_to_save.end(); it++)
		            sheet.addParamToSave(*it);
		        argus::util::saveSheetState(parent_directory+ sheet.getCheckpointDirName(), frame, sheet, &m);
		    }
		} */

		bool getConfigData(const nlohmann::json& data, std::string key, std::string& out, bool path)
		{
			if (data.contains(key))
			{
				out = data[key].get<std::string>();
				if (path)
					out = project_root + out;
				return true;
			}
			return false;
		}
		void loadConfig(std::string path, nlohmann::json& data)
		{
			std::ifstream f(path);
			data = nlohmann::json::parse(f);
		}

		void loadClampedVertices(std::string path, std::vector<int>& clampedVertices)
		{
			std::ifstream file(path);
			std::string str;
			while (std::getline(file, str))
			{
				clampedVertices.push_back(std::atoi(str.c_str()));
			}
		}
		void loadClampedVertices(std::string path, std::vector<std::vector<int>>& clampedVertices)
		{
			std::ifstream file(path);
			std::string str;
			while (std::getline(file, str))
			{
				std::stringstream ss(str);
				std::vector<int> d;
				std::string ch;
				while (std::getline(ss, ch, ' '))
				{
					d.push_back(std::atoi(ch.c_str()));
				}
				clampedVertices.push_back(d);
			}
		}
		// Get pd hessians
		Matrix44 hessianProjection(const Matrix44& mat, bool vbd)
		{
			Matrix44 posHess = mat;
			Vector4 evals;
			Matrix44 D, V;
			// check if the matrix is exactly symmetric
			if (posHess == posHess.transpose())
			{
				Eigen::SelfAdjointEigenSolver<Matrix44> es;
				es.computeDirect(posHess);
				evals = es.eigenvalues();
				V = es.eigenvectors();
			}
			else
			{
#ifdef ARGUS_DEBUG
				std::cerr << "Warning! The matrix passed to hessianProjection is not exactly symmetric!" << std::endl;
				std::cerr << "Using general EigenSolver, performance might be affected." << std::endl;
#endif
				Eigen::EigenSolver<Matrix44> es;
				es.compute(posHess);
				evals = es.eigenvalues().real();
				V = es.eigenvectors().real();
			}

			for (int i = 0; i < evals.size(); i++)
			{
				if (evals(i) < 0)
					evals(i) *= -1.0;
				/* if (vbd)
				 *useless
				{
				    if (evals(i) <= 1e-4)
				        evals(i) = 0;
				} */
			}
			D = evals.asDiagonal();
			posHess = V * D * V.inverse();
			return posHess;
		}

		// Get pd hessians
		Matrix33 hessianProjection(const Matrix33& mat, bool vbd)
		{
			Matrix33 posHess = mat;
			Vector3 evals;
			Matrix33 D, V;
			// check if the matrix is exactly symmetric
			if (posHess == posHess.transpose())
			{
				Eigen::SelfAdjointEigenSolver<Matrix33> es;
				es.computeDirect(posHess);
				evals = es.eigenvalues();
				V = es.eigenvectors();
			}
			else
			{
#ifdef ARGUS_DEBUG
				std::cerr << "Warning! The matrix passed to hessianProjection is not exactly symmetric!" << std::endl;
				std::cerr << "Using general EigenSolver, performance might be affected." << std::endl;
#endif
				Eigen::EigenSolver<Matrix33> es;
				es.compute(posHess);
				evals = es.eigenvalues().real();
				V = es.eigenvectors().real();
			}
			for (int i = 0; i < evals.size(); i++)
			{
				if (evals(i) < 0)
					evals(i) *= -1.0;
				/* if (vbd)
				 *useless
				{
				    if (evals(i) <= 1e-4)
				        evals(i) = 0;
				} */
			}
			D = evals.asDiagonal();
			posHess = V * D * V.inverse();
			return posHess;
		}

		// Get pd hessians
		MatrixNN lowRankApprox(MatrixNN A, bool vbd)
		{
			if ((A - A.transpose()).norm() > A.norm() * 1e-4)
			{
#ifdef ARGUS_DEBUG
				throw "Error! The matrix passed to lowRankApprox is not exactly symmetric!";
#endif
#ifndef ARGUS_QUIET
				std::cerr << "Warning! The matrix passed to lowRankApprox is not exactly symmteric!" << std::endl;
				std::cerr << "Warning! Returning the matrix as it is! Optimization might be affected." << std::endl;
#endif
				return A;
			}
			MatrixNN posHess = A;
			Eigen::SelfAdjointEigenSolver<MatrixNN> es;
			es.compute(posHess);
			VectorN evals = es.eigenvalues();

			for (int i = 0; i < evals.size(); i++)
			{
				if (evals(i) < 0)
					evals(i) *= -1.0;
				if (vbd)
				{
					if (evals(i) <= 1e-4)
						evals(i) = 0;
				}
			}
			MatrixNN D = evals.asDiagonal();
			MatrixNN V = es.eigenvectors();
			posHess = V * D * V.inverse();
			return posHess;
		}

		nlohmann::json loadSimulationConfig(std::string configpath)
		{
			nlohmann::json data;
			try
			{
				argus::util::loadConfig(configpath, data);
			}
			catch (const char* e)
			{
				std::cout << "Error " << e << std::endl;
			}
			catch (...)
			{
				std::cout << "Error while parsing json" << std::endl;
			}
			try
			{
				data["obstacles_present"] = false;
				if (data.contains("obstacles") == true)
				{
					data["obstacles_present"] = true;
				}
				if (data.contains("cloth") == false)
					data["cloth"] = NULL;
				if (data["cloth"].contains("wrinkle_self_coll"))
					simConf["wrinkle_self_coll"] = data["wrinkle_self_coll"].get<nlohmann::json>();

				if (data["cloth"].contains("clamped_vertices_path") == false)
					data["cloth"]["clamped_vertices_path"] = "";
				if (data["cloth"].contains("clamped_vertices_motion_path") == false)
					data["cloth"]["clamped_vertices_motion_path"] = "";

				if (data["cloth"].contains("seam_vertices_path") == false)
					data["cloth"]["seam_vertices_path"] = "";
				if (data.contains("name") == false)
					data["name"] = "general simulation";
				if (data.contains("max_iter") == false)
					data["max_iter"] = 10;

				// Default values
				if (data.contains("rho") == false)
					data["rho"] = 0.002;
				if (data.contains("Y") == false)
					data["Y"] = 1e3;
				if (data.contains("mu") == false)
					data["mu"] = 0.1;
				if (data.contains("h") == false)
					data["h"] = 0.005;
				if (data.contains("pos_dampingK") == false)
					data["pos_dampingK"] = 1e-2;
				if (data.contains("amp_dampingK") == false)
					data["amp_dampingK"] = 1e-4;
				if (data.contains("total_frames") == false)
					data["total_frames"] = 0;
				if (data.contains("frictionK") == false)
					data["frictionK"] = 0;
				if (data.contains("timestepper_config"))
				{
					data["timestepper_config"]["kdamp_amps"] = data["amp_dampingK"].get<Real>();
					data["timestepper_config"]["kdamp_pos"] = data["pos_dampingK"].get<Real>();
				}
			}
			catch (const char* err)
			{
				std::cout << "Error in parsing config data:\n"
				          << err << "!\nExiting.\n";
				throw "Error while loading simulation parameters.\n";
			}
#ifdef ARGUS_DEBUG
			std::cout << "Mesh path " << data["mesh_path"] << std::endl;
#endif
			return data;
		}

		std::pair<MatrixNN, Eigen::MatrixXi> generateTesselatedSphere(const Vector3& center, Real radius, unsigned int subdivisions)
		{
			// Golden ratio
			float t = (1.0f + std::sqrt(5.0f)) / 2.0f;

			// Initialize vertices as std::vector for dynamic growth
			std::vector<Vector3> verts;

			// Create 12 vertices of icosahedron
			verts.push_back(Vector3(-1, t, 0).normalized());
			verts.push_back(Vector3(1, t, 0).normalized());
			verts.push_back(Vector3(-1, -t, 0).normalized());
			verts.push_back(Vector3(1, -t, 0).normalized());

			verts.push_back(Vector3(0, -1, t).normalized());
			verts.push_back(Vector3(0, 1, t).normalized());
			verts.push_back(Vector3(0, -1, -t).normalized());
			verts.push_back(Vector3(0, 1, -t).normalized());

			verts.push_back(Vector3(t, 0, -1).normalized());
			verts.push_back(Vector3(t, 0, 1).normalized());
			verts.push_back(Vector3(-t, 0, -1).normalized());
			verts.push_back(Vector3(-t, 0, 1).normalized());

			// Create 20 triangular faces
			std::vector<Eigen::Vector3i> faces;
			faces.push_back(Eigen::Vector3i(0, 11, 5));
			faces.push_back(Eigen::Vector3i(0, 5, 1));
			faces.push_back(Eigen::Vector3i(0, 1, 7));
			faces.push_back(Eigen::Vector3i(0, 7, 10));
			faces.push_back(Eigen::Vector3i(0, 10, 11));

			faces.push_back(Eigen::Vector3i(1, 5, 9));
			faces.push_back(Eigen::Vector3i(5, 11, 4));
			faces.push_back(Eigen::Vector3i(11, 10, 2));
			faces.push_back(Eigen::Vector3i(10, 7, 6));
			faces.push_back(Eigen::Vector3i(7, 1, 8));

			faces.push_back(Eigen::Vector3i(3, 9, 4));
			faces.push_back(Eigen::Vector3i(3, 4, 2));
			faces.push_back(Eigen::Vector3i(3, 2, 6));
			faces.push_back(Eigen::Vector3i(3, 6, 8));
			faces.push_back(Eigen::Vector3i(3, 8, 9));

			faces.push_back(Eigen::Vector3i(4, 9, 5));
			faces.push_back(Eigen::Vector3i(2, 4, 11));
			faces.push_back(Eigen::Vector3i(6, 2, 10));
			faces.push_back(Eigen::Vector3i(8, 6, 7));
			faces.push_back(Eigen::Vector3i(9, 8, 1));

			// Subdivide
			for (int iter = 0; iter < subdivisions; iter++)
			{
				std::vector<Eigen::Vector3i> newFaces;
				std::map<std::pair<int, int>, int> midpointCache;

				auto getMidpoint = [&](int v1, int v2) -> int
				{
					if (v1 > v2)
						std::swap(v1, v2);
					auto key = std::make_pair(v1, v2);
					auto it = midpointCache.find(key);

					if (it != midpointCache.end())
					{
						return it->second;
					}

					Vector3 mid = ((verts[v1] + verts[v2]) * 0.5f).normalized();
					int index = verts.size();
					verts.push_back(mid);
					midpointCache[key] = index;

					return index;
				};

				for (const auto& face : faces)
				{
					int v1 = face[0], v2 = face[1], v3 = face[2];
					int m1 = getMidpoint(v1, v2);
					int m2 = getMidpoint(v2, v3);
					int m3 = getMidpoint(v3, v1);

					newFaces.push_back(Eigen::Vector3i(v1, m1, m3));
					newFaces.push_back(Eigen::Vector3i(m1, v2, m2));
					newFaces.push_back(Eigen::Vector3i(m3, m2, v3));
					newFaces.push_back(Eigen::Vector3i(m1, m2, m3));
				}

				faces = newFaces;
			}

			// Convert to Eigen matrices: Nverts x 3 and Nfaces x 3
			int numVerts = verts.size();
			int numFaces = faces.size();

			MatrixNN vertices(numVerts, 3);
			Eigen::MatrixXi triangles(numFaces, 3);

			for (int i = 0; i < numVerts; i++)
			{
				vertices.row(i) = (verts[i] * radius + center).transpose();
			}

			for (int i = 0; i < numFaces; i++)
			{
				triangles.row(i) = faces[i].transpose();
			}

			return { vertices, triangles };
		}
	}

}
