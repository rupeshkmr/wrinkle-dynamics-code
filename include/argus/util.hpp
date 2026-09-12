#ifndef UTIL_HPP
#define UTIL_HPP

#include "common.hpp"
#include <time.h>
#include <string>
#include <sstream>
#include <fstream>
#include <igl/readOBJ.h>
#include <nlohmann/json.hpp>

namespace argus
{
	namespace util
	{
		enum class Alignment
		{
			Z = 0,
			X,
			Y,
		};

		// get eigen vector3 from json data
		std::optional<Vector3> getVector3(const nlohmann::json& v);
		// Gradient tests
		void testScalarGradients(std::function<Real(VectorN&, VectorN&)> func);
		Eigen::SparseMatrix<int> buildAdjacencyMatrix(const Eigen::MatrixXi& F, int numVerts);
		// TODO: Cleanup use templates here
		std::vector<argus::Vector3> convertMatrix3NToSTDVec(const Matrix3N& ip);
		std::vector<argus::Vector2> convertMatrix2NToSTDVec(const Matrix2N& ip);
		Eigen::MatrixXi convertTrianglesToMatrix(const std::vector<trimesh::triangle_t>& triangles);
		Matrix2N convertSTDVecMatrix2N(const std::vector<Vector2>& ip);
		Matrix3N convertSTDVecMatrix3N(const std::vector<Vector3>& ip);
		template <typename T>
		Eigen::Matrix<T, 1, Eigen::Dynamic> convertSTDVec2VecN(std::vector<T>& ip)
		{
			Eigen::Matrix<T, 1, Eigen::Dynamic> op;
			op.resize(1, ip.size());
			for (int i = 0; i < ip.size(); i++)
			{
				op(i) = ip.at(i);
			}
			return op;
		}

		template <typename T>
		bool getConfigData(const nlohmann::json& data, std::string key, T& out)
		{
			if (data.contains(key))
			{
				out = data[key].get<T>();
				return true;
			}
			return false;
		}
		bool getConfigData(const nlohmann::json& data, std::string key, std::string& out, bool path = true);

		template <typename T>
		void writeEigenData(std::string filename, T data)
		{
			const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
			std::ofstream file;
			file.open(filename);
			file << data.format(CSVFormat);
		}
		template <typename T>
		bool openEigenData(std::string fileToOpen, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>* data)
		{

			// the inspiration for creating this function was drawn from here (I did NOT copy and paste the code)
			// https://stackoverflow.com/questions/34247057/how-to-read-csv-file-and-assign-to-eigen-matrix
			std::vector<T> matrixEntries;
			try
			{
				// in this object we store the data from the matrix
				std::ifstream matrixDataFile(fileToOpen);
				// this variable is used to store the row of the matrix that contains commas
				std::string matrixRowString;

				// this variable is used to store the matrix entry;
				std::string matrixEntry;

				// this variable is used to track the number of rows
				int matrixRowNumber = 0;

				while (getline(matrixDataFile, matrixRowString)) // here we read a row by row of matrixDataFile and store every line into the string variable matrixRowString
				{
					std::stringstream matrixRowStringStream(matrixRowString); // convert matrixRowString that is a string to a stream variable.

					while (getline(matrixRowStringStream, matrixEntry, ',')) // here we read pieces of the stream matrixRowStringStream until every comma, and store the resulting character into the matrixEntry
					{
						matrixEntries.push_back(stod(matrixEntry)); // here we convert the string to double and fill in the row vector storing all the matrix entries
					}
					matrixRowNumber++; // update the column numbers
				}

				// here we convet the vector variable into the matrix and return the resulting object,
				// note that matrixEntries.data() is the pointer to the first memory location at which the entries of the vector matrixEntries are stored;
				if (matrixRowNumber != 0)
					(*data) = Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(matrixEntries.data(), matrixRowNumber, matrixEntries.size() / matrixRowNumber);
				else
					return false;
				return true;
			}
			catch (const std::exception& e)
			{
				std::cerr << "Failed to read data!\n";
				return false;
			}
		}
		template <typename T>
		bool openEigenData(std::string fileToOpen, Eigen::Vector<T, Eigen::Dynamic>* data)
		{
			MatrixNN readData;
			openEigenData(fileToOpen, &readData);
			try
			{
				readData.reshaped();
				(*data) = Map<Eigen::Vector<T, Eigen::Dynamic>>(readData.data(), readData.size());
				return true;
			}
			catch (const std::exception& e)
			{
				std::cerr << "Failed to read data from " << fileToOpen << std::endl;
				return false;
			}
			return false;
		}
		std::chrono::time_point<std::chrono::steady_clock> captureCurrentTime();
		std::tuple<std::vector<trimesh::triangle_t>, Matrix3N, Matrix2N> generateTesselatedSquare(Vector3 center, size_t n, Real sideLength, Alignment alignment);
		std::pair<MatrixNN, Eigen::MatrixXi> generateTesselatedSphere(const Vector3& center, Real radius, unsigned int subdivisions = 2);
		// Generates a sheet with a cylindrical shape
		// center: center of the cylinder
		// n: number of nodes along each dimension
		// radius: radius of the cylinder
		// theta: the angle denoting what portion of cylinder is generated, 2pi = full cylinder
		// height: height of the cylinder
		std::tuple<std::vector<trimesh::triangle_t>, Matrix3N, Matrix2N> generateCylindricalSheet(Vector3 center, size_t n, Real radius, Real theta, Real height);
		std::tuple<std::vector<trimesh::triangle_t>, Matrix3N, std::vector<trimesh::triangle_t>, Matrix2N> loadMeshWithUVFromObj(std::string fileName);
		std::tuple<std::vector<trimesh::triangle_t>, Matrix3N, Matrix2N> loadMeshFromObj(std::string fileName);

		std::tuple<std::vector<trimesh::triangle_t>, Matrix3N, Matrix2N> loadSheetFromObj(std::string fileName);
		std::tuple<std::vector<trimesh::triangle_t>, Matrix3N> generateThinBeam(Vector3 center, Real length, Real width, Alignment alignment);

		void saveClothData(std::string saveDir, int frame, const Matrix3N& positions, const VectorN& amplitudes, const VectorN& dphis);
		// void loadClothData(std::string loadDir, int frame, Matrix3N& positions, VectorN& amplitudes, VectorN& dphis);
		void loadConfig(std::string path, nlohmann::json& data);
		void loadClampedVertices(std::string path, std::vector<int>& clampedVertices);
		void loadClampedVertices(std::string path, std::vector<std::vector<int>>& clampedVertices);
		MatrixNN lowRankApprox(MatrixNN A, bool vbd = false);
		// TODO: add analytical equation for 3x3 matrices
		Matrix33 hessianProjection(const Matrix33& mat, bool vbd = false);
		Matrix44 hessianProjection(const Matrix44& mat, bool vbd = false);
		nlohmann::json loadSimulationConfig(std::string configpath);

	}
}

#endif /* UTIL_HPP */
