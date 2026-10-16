#pragma once

#include <Discregrid/cubic_lagrange_discrete_grid.hpp>
#include <igl/boundary_loop.h>
#include <igl/cotmatrix_entries.h>
#include <igl/doublearea.h>
#include <igl/file_dialog_open.h>
#include <igl/file_dialog_save.h>
#include <igl/invert_diag.h>
#include <igl/loop.h>
#include <igl/per_vertex_normals.h>
#include <igl/readOBJ.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/writeOBJ.h>

#include <CLI/CLI.hpp>
#include <chrono>
#include <filesystem>
#include <functional>
#include <future>
#include <iostream>
#include <utility>

#include "../../include/CommonTools.h"
#include "../../include/IntrinsicFormula/InterpolateZvals.h"
#include "../../include/IntrinsicFormula/WrinkleEditingCWF.h"
#include "../../include/IntrinsicFormula/WrinkleEditingModel.h"
#include "../../include/MeshLib/MeshUpsampling.h"
#include "../../include/Optimization/NewtonDescent.h"

#include "../../dep/SecStencils/Subd.h"
#include "../../dep/SecStencils/types.h"
#include "../../dep/SecStencils/utils.h"
#include "../../include/IntrinsicFormula/KnoppelStripePatternEdgeOmega.h"
#include "../../include/WrinkleFieldsEditor.h"

#include "../../include/ComplexLoop/ComplexLoop.h"
#include "../../include/ComplexLoop/ComplexLoopZuenko.h"
#include "../../include/LoadSaveIO.h"
#include "../../include/MeshLib/RegionEdition.h"
#include "../../include/SecMeshParsing.h"
#include "../../include/json.hpp"

#include <random>

using Complex = std::complex<double>;

struct RenderConfig {
  int uvPresent = 0;
  int numFrames = 0;
  int upsampleTimes = 0;
  double ampScale = 1.0;
  double learningRate = 1e-4;
  int initialIters = 2000;
  int subsequentIters = 100;
  int intermediateIters = 1;
  double sdfThreshold = 0.0001;
  double timeStep = 0.005;
  double beta = 6666.666666666667;
  double meshDensityRho = 0.1;
  double meshThickness = 0.001;

  std::string optimizationMethod = "gd";
  std::string initializationMethod = "fixed";

  std::string renderFolder;
  std::string workingFolder;

  std::string seamFile = "";
  std::string referenceFile = "";
  std::string baseMeshFile = "";
  std::string facesFile = "";
  std::string verticesFile = "";
  std::string faceOmegasFile = "";
  std::string amplitudesFile = "";
  std::string sdfValueFile = "";
  std::string obstaclePath = "";
  bool isValid() const {
    return numFrames > 0 && upsampleTimes >= 0 && !facesFile.empty() &&
           !verticesFile.empty();
  }
} config;

void copyComplexEigenToStdVector(const Eigen::VectorXcd &psi_eigen,
                                 std::vector<Complex> &zvals_std) {
  int nverts = psi_eigen.size();
  zvals_std.resize(nverts);
  for (int i = 0; i < nverts; ++i) {
    zvals_std[i] = psi_eigen(i);
  }
}

void copyComplexStdVectorToEigen(const std::vector<Complex> &zvals_std,
                                 Eigen::VectorXcd &psi_eigen) {
  int nverts = zvals_std.size();
  psi_eigen.resize(nverts);
  for (int i = 0; i < nverts; ++i) {
    psi_eigen(i) = zvals_std[i];
  }
}

void initializeRandomUnitNormZvals(
    const Eigen::MatrixXd &triV, std::vector<std::complex<double>> &initZvals) {
  int nverts = triV.rows();

  initZvals.resize(nverts);

  std::random_device rd;
  std::mt19937 gen(rd());

  std::uniform_real_distribution<> distrib(0.0, 2.0 * M_PI);

  for (int i = 0; i < nverts; ++i) {
    double angle = distrib(gen);

    double magnitude = 1.0;

    double real_part = magnitude * std::cos(angle);
    double imag_part = magnitude * std::sin(angle);

    initZvals[i] = std::complex<double>(real_part, imag_part);
  }
}

void initializeFixedUnitNormZvals(
    const Eigen::MatrixXd &triV, std::vector<std::complex<double>> &initZvals) {
  int nverts = triV.rows();

  initZvals.resize(nverts);

  for (int i = 0; i < nverts; ++i) {
    double angle = 0;

    double magnitude = 1.0;

    double real_part = magnitude * std::cos(angle);
    double imag_part = magnitude * std::sin(angle);

    initZvals[i] = std::complex<double>(real_part, imag_part);
  }
}

void convertAmplitudeToComplexZ(const Eigen::VectorXd &curAmp,
                                std::vector<Complex> &augmentedAmp) {
  size_t nverts = curAmp.size();

  augmentedAmp.resize(nverts);

  for (size_t j = 0; j < nverts; ++j) {
    // Complex Z = (Amplitude, 0.0)
    augmentedAmp[j] = Complex(curAmp(j), 0.0);
  }
}

void augmentUpsampledZList(std::vector<std::complex<double>> &curUpZ,
                           std::vector<std::complex<double>> &curUpAmp) {
  if (curUpZ.size() != curUpAmp.size()) {
    std::cerr << "ERROR: Size mismatch in augmentUpsampledZList! curUpZ size: "
              << curUpZ.size() << ", curUpAmp size: " << curUpAmp.size()
              << std::endl;
    return;
  }

  for (size_t j = 0; j < curUpZ.size(); ++j) {
    double scalarAmplitude = curUpAmp[j].real();
    curUpZ[j] *= scalarAmplitude;
  }
}