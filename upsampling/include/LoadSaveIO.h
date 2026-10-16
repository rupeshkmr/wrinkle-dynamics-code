#pragma once
#include "CommonTools.h"
#include <fstream>
#include <iostream>

bool loadEdgeOmega(const std::string &filename, const int &nlines,
                   Eigen::VectorXd &edgeOmega);
bool loadVertexZvals(const std::string &filePath, const int &nlines,
                     std::vector<std::complex<double>> &zvals);
bool loadVertexAmp(const std::string &filePath, const int &nlines,
                   Eigen::VectorXd &amp);

bool saveEdgeOmega(const std::string &filename,
                   const Eigen::VectorXd &edgeOmega);
bool saveVertexZvals(const std::string &filePath,
                     const std::vector<std::complex<double>> &zvals);
bool saveVertexAmp(const std::string &filePath, const Eigen::VectorXd &amp);

/**
 * @brief Reads a column of numbers from a text file into an Eigen Vector.
 * @tparam T The scalar type (int, float, double, etc.)
 * @param filePath Path to the .txt file
 * @param vec The Eigen vector to be populated
 * @return true if successful, false otherwise
 */
template <typename T>
bool readEigenVectorFromTxt(const std::string &filePath,
                            Eigen::Matrix<T, Eigen::Dynamic, 1> &vec) {
  std::ifstream file(filePath);

  if (!file.is_open()) {
    std::cerr << "Error: Could not open file: " << filePath << std::endl;
    return false;
  }

  std::vector<T> buffer;
  T value;

  // Efficiently read values until EOF or invalid input
  while (file >> value) {
    buffer.push_back(value);
  }

  if (buffer.empty()) {
    std::cerr << "Warning: File was empty or contained no valid data: "
              << filePath << std::endl;
    return false;
  }

  // Resize the Eigen vector to match the collected data
  vec.resize(static_cast<Eigen::Index>(buffer.size()));

  // Transfer data using Eigen::Map for performance
  // This avoids a manual loop and is highly optimized by Eigen
  vec = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1>>(buffer.data(),
                                                        buffer.size());

  return true;
}