#include "optimizationFunctions.hpp"
#include "readingFunctions.hpp"
#include "vbdOptimization.hpp"
#include "viz.hpp"
// std::vector<Eigen::MatrixXd> triV, upsampledTriV;
Eigen::MatrixXd curTriV, curUpsampledTriV, refTriV;
Eigen::MatrixXd baseTriV;

Eigen::MatrixXi triF, upsampledTriF, refTriF;

MeshConnectivity triMesh, upsampledTriMesh, refTriMesh;
// MeshConnectivity baseMesh;

Mesh secMesh, subSecMesh;

std::map<int, int> seamMap;
Eigen::VectorXi seamVerts;
Discregrid::CubicLagrangeDiscreteGrid *curSdf;

// need to add edge lengths and edge weights

// initial information
std::vector<std::complex<double>> prevZvals;
std::vector<std::complex<double>> prev2Zvals;

std::vector<std::vector<int>> vert_to_faces;
std::vector<std::vector<int>> V_by_color;

Eigen::VectorXd curAmp, refAmp;
// Eigen::VectorXd uvCurAmp;
Eigen::VectorXd prevAmp;

Eigen::VectorXd curOmega, refOmega;
Eigen::VectorXd prevOmega;

Eigen::VectorXd vertexMasses;
Eigen::VectorXd edgeAreas;
Eigen::VectorXd edgeLengths;

Eigen::MatrixXd curFaceOmegaList2Dinput;

Eigen::MatrixXd curFaceOmegaList3Dinput;
Eigen::MatrixXd edgeGradPhis;
Eigen::MatrixXd vertexGradPhis;
// Eigen::VectorXd curOmega3D;

Eigen::MatrixXd curVertexOmega;

Eigen::MatrixXd curFaceOmega;

std::vector<std::complex<double>> curZ;

Eigen::VectorXd curSubOmega;

Eigen::MatrixXd curSubFaceOmega;

Eigen::VectorXd curPhaseFields;

Eigen::VectorXd curAmpFields;

std::vector<std::complex<double>> curUpZ, curUpAmp;

// std::vector<Eigen::MatrixXd> wrinkledVList;
Eigen::MatrixXd curWrinkledV;

std::vector<SparseMatrixX> S0List, S1List;

std::shared_ptr<ComplexLoop> complexLoopGlobal =
    std::make_shared<ComplexLoopZuenko>();

std::vector<std::vector<int>> vertNeiFaces;

std::vector<Mesh> meshList;

double globalAmpMax = 1;
double globalAmpMin = 0;

bool isFixedBnd = true;

bool isUseV2 = false;

// Default arguments
struct {
  std::string paramsFile = "";
  bool saveCheckpoint = false;
  bool headless = false; // no polyscope window, frames are not kept in memory
  int saveEvery = 1;     // save every Nth frame (the last frame is always saved)
} args;

void loadRenderConfig(const std::string &jsonFilePath) {
  using json = nlohmann::json;
  std::ifstream file(jsonFilePath);
  if (!file.is_open()) {
    std::cerr << "Error: Could not open JSON config file: " << jsonFilePath
              << std::endl;
    return;
  }

  json j;
  try {
    file >> j;
  } catch (const json::parse_error &e) {
    std::cerr << "Error parsing JSON file: " << e.what() << std::endl;
    return;
  }

  std::string filePath = jsonFilePath;

  std::replace(filePath.begin(), filePath.end(), '\\', '/');
  size_t last_slash = filePath.rfind("/");

  if (last_slash != std::string::npos) {
    config.workingFolder = filePath.substr(0, last_slash + 1);
  } else {
    config.workingFolder = "./";
  }

  // 1 means using uv coordinates, 0 means using the basemesh
  config.uvPresent = j.value("uvPresent", 0);
  config.numFrames = j.value("numFrames", 0);
  config.upsampleTimes = j.value("upsampledTimes", 0);
  config.ampScale = j.value("ampScale", 1.0);
  config.initialIters = j.value("initialIters", 2000);
  config.subsequentIters = j.value("subsequentIters", 100);
  config.intermediateIters = j.value("intermediateIters", 1);
  config.learningRate = j.value("learningRate", 1e-4);
  config.sdfThreshold = j.value("sdfThreshold", 0.0001);
  config.timeStep = j.value("timeStep", 0.005);
  config.beta = j.value("beta", 6666.666666666667);
  config.meshDensityRho = j.value("meshDensityRho", 0.1);
  config.meshThickness = j.value("meshThickness", 0.001);

  config.optimizationMethod = j.value("optimizationMethod", "gd");
  config.initializationMethod = j.value("initializationMethod", "fixed");

  config.renderFolder = j.value("renderFolder", "");

  const auto &files = j.value("files", json::object());

  config.seamFile = files.value("seamFile", "");
  config.baseMeshFile = files.value("baseMeshFile", "");
  config.facesFile = files.value("facesFile", "");
  config.verticesFile = files.value("verticesFile", "");
  config.faceOmegasFile = files.value("faceOmegasFile", "");
  config.amplitudesFile = files.value("amplitudesFile", "");
  config.sdfValueFile = files.value("sdfValueFile", "");
  config.referenceFile = files.value("referenceMeshFile", "");
  config.obstaclePath = files.value("obstaclePath", "");

  if (!config.isValid()) {
    std::cerr << "Error: JSON configuration is missing required parameters "
                 "(numFrames, upsampleTimes, vertices_prefix, or faces_file)."
              << std::endl;
  } else {
    std::cout << "Configuration loaded successfully from "
              << config.workingFolder << std::endl;
  }
}

Eigen::Matrix3d computePlaneToXZRotation(const Eigen::MatrixXd &V,
                                         const Eigen::MatrixXi &F) {
  Eigen::Vector3d n(0, 0, 0);
  Eigen::Vector3d t(0, 0, 0);
  Eigen::Vector3d p0(0, 0, 0);

  bool found_face = false;

  // --- find a non-degenerate face ---
  for (int f = 0; f < F.rows(); ++f) {
    Eigen::Vector3d v0 = V.row(F(f, 0));
    Eigen::Vector3d v1 = V.row(F(f, 1));
    Eigen::Vector3d v2 = V.row(F(f, 2));

    Eigen::Vector3d e1 = v1 - v0;
    Eigen::Vector3d e2 = v2 - v0;

    Eigen::Vector3d fn = e1.cross(e2);
    if (fn.norm() > 1e-8) {
      n = fn.normalized();
      p0 = v0;

      // pick a stable in-plane direction
      t = e1 - e1.dot(n) * n;
      t.normalize();

      found_face = true;
      break;
    }
  }

  if (!found_face) {
    throw std::runtime_error(
        "computePlaneToXZRotation: Could not find a non-degenerate face.");
  }

  // --- Planarity check ---
  const double plane_tol = 1e-6;
  // for (int i = 0; i < V.rows(); ++i) {
  //     Eigen::Vector3d vi = V.row(i);
  //     double dist = std::abs(n.dot(vi - p0));
  //     if (dist > plane_tol) {
  //         throw std::runtime_error(
  //             "computePlaneToXZRotation: Mesh is not planar "
  //             "(vertex deviates from plane by more than tolerance).");
  //     }
  // }

  // complete orthonormal basis
  Eigen::Vector3d b = n.cross(t);

  // source frame (mesh plane)
  Eigen::Matrix3d from;
  from.col(0) = t; // → X
  from.col(1) = n; // → Y (plane normal)
  from.col(2) = b; // → Z

  // target frame (XZ plane)
  Eigen::Matrix3d to = Eigen::Matrix3d::Identity();

  // rotation that maps plane → XZ
  return to * from.transpose();
}

Eigen::Vector3d
lineFieldAddition(const std::vector<Eigen::Vector3d> &grad_phi_list) {
  if (grad_phi_list.empty()) {
    return Eigen::Vector3d::Zero();
  }

  Eigen::Vector3d total_vec = Eigen::Vector3d::Zero();

  for (const auto &vec : grad_phi_list) {
    // Assume the vector is approximately in the X-Z plane (common in wrinkle
    // fields)
    double norm = vec.norm();
    // atan2(y, x) -> using Z as the vertical component (y) and X as the
    // horizontal (x)
    double angle = std::atan2(vec[2], vec[0]);

    // Map to double-angle space: [norm*cos(2*angle), 0, norm*sin(2*angle)]
    total_vec[0] +=
        norm * std::cos(2.0 * angle); // X component of doubled vector
    // total_vec[1] += 0.0;                       // Y component remains 0
    total_vec[2] +=
        norm * std::sin(2.0 * angle); // Z component of doubled vector
  }

  // Average the doubled vector
  total_vec /= (double)grad_phi_list.size();

  // Map back to single-angle space
  double final_norm = total_vec.norm();
  // atan2(y, x) -> using Z as y, X as x
  double final_angle = std::atan2(total_vec[2], total_vec[0]);

  // Halve the angle: final_angle / 2
  Eigen::Vector3d final_vec;
  final_vec[0] = final_norm * std::cos(final_angle / 2.0); // X component
  final_vec[1] = 0.0;                                      // Y component
  final_vec[2] = final_norm * std::sin(final_angle / 2.0); // Z component

  // if(final_vec[2]<0) final_vec = -final_vec;

  return final_vec;
}

void getGradPhisPerVertexLineField(const Eigen::MatrixXd &face_grad_phis,
                                   const Eigen::MatrixXd &V,
                                   const Eigen::MatrixXi &F,
                                   Eigen::MatrixXd &vertex_grad_phis) {
  int num_vertices = V.rows();
  int num_faces = F.rows();

  // --- Compute plane ↔ XZ rotation ---
  Eigen::Matrix3d R = computePlaneToXZRotation(V, F);
  Eigen::Matrix3d R_inv = R.transpose();

  // --- Rotate face gradients into XZ ---
  std::vector<Eigen::Vector3d> face_grad_xz(num_faces);
  for (int f = 0; f < num_faces; ++f) {
    face_grad_xz[f] = R * face_grad_phis.row(f).transpose();
  }

  // --- Accumulate per-vertex ---
  std::vector<std::vector<Eigen::Vector3d>> vertex_lists(num_vertices);
  for (int f = 0; f < num_faces; ++f) {
    for (int k = 0; k < 3; ++k) {
      int v = F(f, k);
      vertex_lists[v].push_back(face_grad_xz[f]);
    }
  }

  // --- Line-field averaging (in XZ) ---
  vertex_grad_phis.resize(num_vertices, 3);
  for (int v = 0; v < num_vertices; ++v) {
    Eigen::Vector3d avg_xz = lineFieldAddition(vertex_lists[v]);

    // Rotate result back to original plane
    Eigen::Vector3d avg_world = R_inv * avg_xz;
    vertex_grad_phis.row(v) = avg_world.transpose();
  }

  // if (config.seamFile != "") {
  //   for (auto &pair : seamMap) {
  //     int v1 = pair.first;
  //     int v2 = pair.second;

  //     Eigen::Vector3d grad_phi_v1 = vertex_grad_phis.row(v1).transpose();
  //     Eigen::Vector3d grad_phi_v2 = vertex_grad_phis.row(v2).transpose();

  //     Eigen::Vector3d avg_grad_phi =
  //         lineFieldAddition({grad_phi_v1, grad_phi_v2});
  //     // Eigen::Vector3d avg_grad_phi = (grad_phi_v1 + grad_phi_v2)
  //     // / 2.0;

  //     vertex_grad_phis.row(v1) = avg_grad_phi.transpose();
  //     vertex_grad_phis.row(v2) = avg_grad_phi.transpose();
  //   }
  // }
}

Eigen::Matrix3d computeEdgeLocalPlaneToXZRotation(
    const Eigen::Vector3d &v0, const Eigen::Vector3d &v1,
    const Eigen::Vector3d &v2_0, const Eigen::Vector3d &v2_1_reflected) {
  // Compute normal for face0 (v2_1_reflected already lies in face0's plane)
  Eigen::Vector3d e1 = v1 - v0;
  Eigen::Vector3d e2_0 = v2_0 - v0;
  Eigen::Vector3d n0 = e1.cross(e2_0);

  if (n0.norm() < 1e-12) {
    // Degenerate triangle, return identity
    return Eigen::Matrix3d::Identity();
  }
  n0.normalize();

  // Create orthonormal basis with edge as X-axis and n0 as Z-axis
  Eigen::Vector3d edge_dir = e1.normalized();
  Eigen::Vector3d y_axis = n0;
  Eigen::Vector3d z_axis = edge_dir.cross(y_axis);

  if (z_axis.norm() < 1e-12) {
    // Edge and normal are parallel, pick arbitrary perpendicular to edge as Z
    Eigen::Vector3d perp;
    if (std::abs(edge_dir[0]) < 0.9) {
      perp = Eigen::Vector3d(1, 0, 0);
    } else {
      perp = Eigen::Vector3d(0, 1, 0);
    }
    z_axis = perp - (perp.dot(edge_dir)) * edge_dir;
    z_axis.normalize();
  } else {
    z_axis.normalize();
  }

  // Build rotation matrix: from world to XZ (columns are basis vectors)
  Eigen::Matrix3d Q;
  Q.col(0) = edge_dir;
  Q.col(1) = y_axis;
  Q.col(2) = z_axis;

  // Rotation from world to local XZ plane
  return Q.transpose();
}
void fixOmegas(Eigen::MatrixXd baseTriV, Eigen::VectorXd &curOmega,
               std::vector<std::complex<double>> &zvals) {
  int num_edges = triMesh.nEdges();
  for (int e = 0; e < num_edges; e++) {
    int v0_idx = triMesh.edgeVertex(e, 0);
    int v1_idx = triMesh.edgeVertex(e, 1);
    if (v1_idx < v0_idx)
      std::swap(v0_idx, v1_idx);
    std::complex<double> z0 = zvals[v0_idx];
    std::complex<double> z1 = zvals[v1_idx];
    std::complex<double> expomega = std::polar(1.0, curOmega(e));
    std::complex<double> negexpomega = std::polar(1.0, -curOmega(e));
    double err1 = std::abs(z0 * expomega - z1);
    double err2 = std::abs(z0 * negexpomega - z1);
    if (err1 > err2)
      curOmega(e) *= -1.0;
  }
}

void faceToEdgeOmegas(Eigen::MatrixXd curFaceOmegaList2Dinput,
                      Eigen::MatrixXd baseTriV, Eigen::VectorXd &curOmega) {

  int num_edges = triMesh.nEdges();
  // std::cout << "[faceToEdgeOmegas] START: num_edges = " << num_edges <<
  // std::endl; std::cout << "[faceToEdgeOmegas] baseTriV size: " <<
  // baseTriV.rows() << " x " << baseTriV.cols() << std::endl; std::cout <<
  // "[faceToEdgeOmegas] curFaceOmegaList2Dinput size: " <<
  // curFaceOmegaList2Dinput.rows() << " x " << curFaceOmegaList2Dinput.cols()
  // << std::endl;
  curOmega.setZero(num_edges);

  for (int e = 0; e < num_edges; e++) {
    int face0 = triMesh.edgeFace(e, 0);
    int face1 = triMesh.edgeFace(e, 1);

    // std::cout << "\n[Edge " << e << "] face0=" << face0 << ", face1=" <<
    // face1 << std::endl;

    // If both faces are invalid, skip this edge
    if (face0 < 0 && face1 < 0) {
      // std::cout << "  -> Both faces invalid, skipping" << std::endl;
      continue;
    }

    // If only face0 is valid, use it as the primary face
    // If only face1 is valid, use it as the primary face
    // If both are valid, use face0 as primary
    int primaryFace = (face0 >= 0) ? face0 : face1;
    int secondaryFace = (face0 >= 0) ? face1 : face0;

    // std::cout << "  primaryFace=" << primaryFace << ", secondaryFace=" <<
    // secondaryFace << std::endl;

    // Get edge vertices
    int v0_idx = triMesh.edgeVertex(e, 0);
    int v1_idx = triMesh.edgeVertex(e, 1);

    // std::cout << "  v0_idx=" << v0_idx << ", v1_idx=" << v1_idx << std::endl;

    Eigen::Vector3d v0 = baseTriV.row(v0_idx).transpose();
    Eigen::Vector3d v1 = baseTriV.row(v1_idx).transpose();
    // std::cout << "  v0=" << v0.transpose() << std::endl;
    // std::cout << "  v1=" << v1.transpose() << std::endl;

    Eigen::Vector3d edge_vec = v1 - v0;
    double edge_len = edge_vec.norm();

    // std::cout << "  edge_vec=" << edge_vec.transpose() << ", edge_len=" <<
    // edge_len << std::endl;

    if (edge_len < 1e-12) {
      // std::cout << "  -> edge_len too small, skipping" << std::endl;
      continue;
    }
    Eigen::Vector3d edge_dir = edge_vec / edge_len;
    // std::cout << "  edge_dir=" << edge_dir.transpose() << std::endl;

    // Get third vertex of primary face, preserving triF vertex order
    // Find the position of v0_idx and v1_idx in triF, then get the third vertex
    // accordingly
    int pos_v0 = -1, pos_v1 = -1;
    for (int k = 0; k < 3; k++) {
      if (triF(primaryFace, k) == v0_idx)
        pos_v0 = k;
      if (triF(primaryFace, k) == v1_idx)
        pos_v1 = k;
    }

    // Find the third position (the one that's not pos_v0 or pos_v1)
    int pos_v2_0 = -1;
    for (int k = 0; k < 3; k++) {
      if (k != pos_v0 && k != pos_v1) {
        pos_v2_0 = k;
        break;
      }
    }

    if (pos_v0 < 0 || pos_v1 < 0 || pos_v2_0 < 0) {
      std::cout
          << "  -> Could not find consistent vertex positions in triF, skipping"
          << std::endl;
      continue;
    }

    int v2_0_idx = triF(primaryFace, pos_v2_0);
    // std::cout << "  v2_0_idx (from primaryFace)=" << v2_0_idx << " at triF
    // position " << pos_v2_0 << std::endl;

    Eigen::Vector3d v2_0 = baseTriV.row(v2_0_idx).transpose();
    // std::cout << "  v2_0=" << v2_0.transpose() << std::endl;

    // Compute plane of primary face
    Eigen::Vector3d e1_0 = v1 - v0;
    Eigen::Vector3d e2_0 = v2_0 - v0;
    Eigen::Vector3d n0 = e1_0.cross(e2_0).normalized();
    // std::cout << "  n0 (normalized)=" << n0.transpose() << std::endl;

    // If secondary face exists, reflect v2_secondary into the plane of primary
    // face; otherwise use v2_0
    Eigen::Vector3d v2_1_reflected = v2_0; // Default to v2_0 for boundary edges
    // std::cout << "  v2_1_reflected initial (default)=" <<
    // v2_1_reflected.transpose() << std::endl;

    if (secondaryFace >= 0) {
      // std::cout << "  Processing secondaryFace..." << std::endl;
      // Find the third vertex in secondary face, preserving triF vertex order
      int pos_v0_sec = -1, pos_v1_sec = -1;
      for (int k = 0; k < 3; k++) {
        if (triF(secondaryFace, k) == v0_idx)
          pos_v0_sec = k;
        if (triF(secondaryFace, k) == v1_idx)
          pos_v1_sec = k;
      }

      int pos_v2_1 = -1;
      for (int k = 0; k < 3; k++) {
        if (k != pos_v0_sec && k != pos_v1_sec) {
          pos_v2_1 = k;
          break;
        }
      }

      int v2_1_idx = (pos_v2_1 >= 0) ? triF(secondaryFace, pos_v2_1) : -1;
      // std::cout << "    v2_1_idx=" << v2_1_idx << (pos_v2_1 >= 0 ? (" at triF
      // position " + std::to_string(pos_v2_1)) : " (not found)") << std::endl;

      if (v2_1_idx >= 0) {
        Eigen::Vector3d v2_1 = baseTriV.row(v2_1_idx).transpose();
        // std::cout << "    v2_1=" << v2_1.transpose() << std::endl;

        // Project v2_1 onto edge
        double proj_along_edge =
            (v2_1 - v0).dot(edge_vec) / (edge_len * edge_len);
        Eigen::Vector3d v2_1_proj_on_edge = v0 + proj_along_edge * edge_vec;
        // std::cout << "    proj_along_edge=" << proj_along_edge << std::endl;
        // std::cout << "    v2_1_proj_on_edge=" <<
        // v2_1_proj_on_edge.transpose() << std::endl;

        // Perpendicular from v2_1 to edge
        Eigen::Vector3d perp_v2_1 = v2_1 - v2_1_proj_on_edge;
        double perp_dist = perp_v2_1.norm();
        // std::cout << "    perp_v2_1=" << perp_v2_1.transpose() << ",
        // perp_dist=" << perp_dist << std::endl;

        if (perp_dist > 1e-12) {
          Eigen::Vector3d perp_dir_v2_1 = perp_v2_1 / perp_dist;
          Eigen::Vector3d u = edge_dir.cross(perp_dir_v2_1).normalized();
          // std::cout << "    perp_dir_v2_1=" << perp_dir_v2_1.transpose() <<
          // std::endl; std::cout << "    u=" << u.transpose() << std::endl;

          // Find where the circle (v2_1 rotated about edge) intersects the
          // plane n0
          double a = n0.dot(perp_dir_v2_1);
          double b = n0.dot(u);
          double c = n0.dot(v2_1_proj_on_edge - v0);
          // std::cout << "    a=" << a << ", b=" << b << ", c=" << c <<
          // std::endl;

          double amplitude = std::sqrt(a * a + b * b);
          // std::cout << "    amplitude=" << amplitude << std::endl;

          if (amplitude > 1e-12) {
            double rhs_norm = -c / perp_dist;
            // std::cout << "    rhs_norm=" << rhs_norm << std::endl;

            if (std::abs(rhs_norm) <= amplitude) {
              double cos_theta_minus_phi = rhs_norm / amplitude;
              double phi = std::atan2(b, a);
              double theta_minus_phi =
                  std::acos(std::max(-1.0, std::min(1.0, cos_theta_minus_phi)));
              // std::cout << "    cos_theta_minus_phi=" << cos_theta_minus_phi
              // << std::endl; std::cout << "    phi=" << phi << ",
              // theta_minus_phi=" << theta_minus_phi << std::endl;

              double theta_1 = theta_minus_phi + phi;
              double theta_2 = -theta_minus_phi + phi;
              // std::cout << "    theta_1=" << theta_1 << ", theta_2=" <<
              // theta_2 << std::endl;

              auto computePoint = [&](double theta) {
                return v2_1_proj_on_edge +
                       perp_dist * (std::cos(theta) * perp_dir_v2_1 +
                                    std::sin(theta) * u);
              };

              Eigen::Vector3d p1_cand = computePoint(theta_1);
              Eigen::Vector3d p2_cand = computePoint(theta_2);
              // std::cout << "    p1_cand=" << p1_cand.transpose() <<
              // std::endl; std::cout << "    p2_cand=" << p2_cand.transpose()
              // << std::endl;

              double dist1 = (p1_cand - v2_0).norm();
              double dist2 = (p2_cand - v2_0).norm();
              // std::cout << "    dist1 (from v2_0)=" << dist1 << ", dist2
              // (from v2_0)=" << dist2 << std::endl;

              v2_1_reflected = (dist1 > dist2) ? p1_cand : p2_cand;
              // std::cout << "    v2_1_reflected (chosen)=" <<
              // v2_1_reflected.transpose() << std::endl;
            } else {
              // std::cout << "    -> |rhs_norm| > amplitude, fallback" <<
              // std::endl;
              v2_1_reflected = v2_1;
            }
          } else {
            // std::cout << "    -> amplitude too small, fallback" << std::endl;
            v2_1_reflected = v2_1;
          }
        } else {
          // std::cout << "    -> perp_dist too small, fallback" << std::endl;
          v2_1_reflected = v2_1;
        }
      }
    }
    // std::cout << "  v2_1_reflected final=" << v2_1_reflected.transpose() <<
    // std::endl;

    // Compute rotation to XZ plane based only on these two faces
    Eigen::Matrix3d R =
        computeEdgeLocalPlaneToXZRotation(v0, v1, v2_0, v2_1_reflected);
    // std::cout << "  Rotation matrix R:\n" << R << std::endl;

    // Get vertices in XZ space (using reflected v2_1 if face1 exists)
    Eigen::Vector3d v0_xz = R * v0;
    Eigen::Vector3d v1_xz = R * v1;
    Eigen::Vector3d v2_0_xz = R * v2_0;
    // std::cout << "  v0_xz=" << v0_xz.transpose() << std::endl;
    // std::cout << "  v1_xz=" << v1_xz.transpose() << std::endl;
    // std::cout << "  v2_0_xz=" << v2_0_xz.transpose() << std::endl;

    // Reconstruct omegas in XZ space using triF vertex ordering
    std::vector<Eigen::Vector3d> omega_list_xz;

    // Primary face - use triF ordering
    // Get all three vertices in triF order
    Eigen::Vector3d v_triF0_xz =
        R * baseTriV.row(triF(primaryFace, 0)).transpose();
    Eigen::Vector3d v_triF1_xz =
        R * baseTriV.row(triF(primaryFace, 1)).transpose();
    Eigen::Vector3d v_triF2_xz =
        R * baseTriV.row(triF(primaryFace, 2)).transpose();

    // Build basis from triF ordering
    Eigen::Vector3d e1_0_xz = v_triF1_xz - v_triF0_xz;
    Eigen::Vector3d e2_0_xz = v_triF2_xz - v_triF0_xz;
    // std::cout << "  [Primary Face] e1_0_xz (triF order)=" <<
    // e1_0_xz.transpose() << std::endl; std::cout << "  [Primary Face] e2_0_xz
    // (triF order)=" << e2_0_xz.transpose() << std::endl;

    Eigen::Vector2d c0;
    c0(0) = curFaceOmegaList2Dinput(primaryFace, 0);
    c0(1) = curFaceOmegaList2Dinput(primaryFace, 1);
    // std::cout << "  [Primary Face] c0=" << c0.transpose() << std::endl;

    // Solve Gram matrix system: G * x = c
    Eigen::Matrix2d G0;
    G0(0, 0) = e1_0_xz.dot(e1_0_xz);
    G0(0, 1) = e1_0_xz.dot(e2_0_xz);
    G0(1, 0) = e2_0_xz.dot(e1_0_xz);
    G0(1, 1) = e2_0_xz.dot(e2_0_xz);
    // std::cout << "  [Primary Face] G0=\n" << G0 << std::endl;
    // std::cout << "  [Primary Face] det(G0)=" << G0.determinant() <<
    // std::endl;

    Eigen::Vector3d omega0_xz = Eigen::Vector3d::Zero();
    if (std::abs(G0.determinant()) > 1e-15) {
      Eigen::Vector2d x0 = G0.inverse() * c0;
      // std::cout << "  [Primary Face] x0=" << x0.transpose() << std::endl;
      omega0_xz = x0(0) * e1_0_xz + x0(1) * e2_0_xz;
      // std::cout << "  [Primary Face] omega0_xz=" << omega0_xz.transpose() <<
      // std::endl;
    } else {
      // std::cout << "  [Primary Face] det(G0) too small, omega0_xz=zero" <<
      // std::endl;
    }
    omega_list_xz.push_back(omega0_xz);

    // Secondary face (if exists)
    if (secondaryFace >= 0) {
      // std::cout << "  [Secondary Face] Processing..." << std::endl;
      // Note: v2_1_idx was already found and used during reflection above
      // We can reuse that same logic but we already have v2_1_reflected
      // Just continue with omega reconstruction using v2_1_reflected

      // The v2_1_idx should be the same as found earlier in the reflection
      // block For safety, let's look it up again with the same logic
      int pos_v0_sec = -1, pos_v1_sec = -1;
      for (int k = 0; k < 3; k++) {
        if (triF(secondaryFace, k) == v0_idx)
          pos_v0_sec = k;
        if (triF(secondaryFace, k) == v1_idx)
          pos_v1_sec = k;
      }

      int pos_v2_1_omega = -1;
      for (int k = 0; k < 3; k++) {
        if (k != pos_v0_sec && k != pos_v1_sec) {
          pos_v2_1_omega = k;
          break;
        }
      }

      int v2_1_idx =
          (pos_v2_1_omega >= 0) ? triF(secondaryFace, pos_v2_1_omega) : -1;
      // std::cout << "    v2_1_idx=" << v2_1_idx << std::endl;

      if (v2_1_idx >= 0) {
        // Secondary face - use triF ordering
        Eigen::Vector3d v_triF0_sec_xz =
            R * baseTriV.row(triF(secondaryFace, 0)).transpose();
        Eigen::Vector3d v_triF1_sec_xz =
            R * baseTriV.row(triF(secondaryFace, 1)).transpose();
        Eigen::Vector3d v_triF2_sec_xz =
            R * baseTriV.row(triF(secondaryFace, 2)).transpose();

        // Build basis from triF ordering
        Eigen::Vector3d e1_1_xz = v_triF1_sec_xz - v_triF0_sec_xz;
        Eigen::Vector3d e2_1_xz = v_triF2_sec_xz - v_triF0_sec_xz;
        // std::cout << "  [Secondary Face] e1_1_xz (triF order)=" <<
        // e1_1_xz.transpose() << std::endl; std::cout << "  [Secondary Face]
        // e2_1_xz (triF order)=" << e2_1_xz.transpose() << std::endl;

        Eigen::Vector2d c1;
        c1(0) = curFaceOmegaList2Dinput(secondaryFace, 0);
        c1(1) = curFaceOmegaList2Dinput(secondaryFace, 1);
        // std::cout << "  [Secondary Face] c1=" << c1.transpose() << std::endl;

        // Solve Gram matrix system: G * x = c
        Eigen::Matrix2d G1;
        G1(0, 0) = e1_1_xz.dot(e1_1_xz);
        G1(0, 1) = e1_1_xz.dot(e2_1_xz);
        G1(1, 0) = e2_1_xz.dot(e1_1_xz);
        G1(1, 1) = e2_1_xz.dot(e2_1_xz);
        // std::cout << "  [Secondary Face] G1=\n" << G1 << std::endl;
        // std::cout << "  [Secondary Face] det(G1)=" << G1.determinant() <<
        // std::endl;

        Eigen::Vector3d omega1_xz = Eigen::Vector3d::Zero();
        if (std::abs(G1.determinant()) > 1e-15) {
          Eigen::Vector2d x1 = G1.inverse() * c1;
          // std::cout << "  [Secondary Face] x1=" << x1.transpose() <<
          // std::endl;
          omega1_xz = x1(0) * e1_1_xz + x1(1) * e2_1_xz;
          // std::cout << "  [Secondary Face] omega1_xz=" <<
          // omega1_xz.transpose() << std::endl;
        } else {
          // std::cout << "  [Secondary Face] det(G1) too small, omega1_xz=zero"
          // << std::endl;
        }
        omega_list_xz.push_back(omega1_xz);
      }
    }

    // std::cout << "  omega_list_xz size: " << omega_list_xz.size() <<
    // std::endl; for(size_t i = 0; i < omega_list_xz.size(); i++) {
    //     std::cout << "    omega_list_xz[" << i << "]=" <<
    //     omega_list_xz[i].transpose() << std::endl;
    // }

    // Line field addition in XZ space
    Eigen::Vector3d omega_avg_xz = lineFieldAddition(omega_list_xz);
    // std::cout << "  omega_avg_xz (after lineFieldAddition)=" <<
    // omega_avg_xz.transpose() << std::endl;

    // Project onto edge (in XZ space)
    Eigen::Vector3d edge_xz = R * edge_vec;
    double edge_xz_dot = omega_avg_xz.dot(edge_xz);
    // std::cout << "  edge_xz=" << edge_xz.transpose() << std::endl;
    // std::cout << "  omega_avg_xz.dot(edge_xz)=" << edge_xz_dot << std::endl;

    curOmega(e) = edge_xz_dot;
    // std::cout << "  curOmega(" << e << ")=" << curOmega(e) << std::endl;
  }

  // std::cout << "[faceToEdgeOmegas] END - curOmega values:" << std::endl;
  // for(int e = 0; e < std::min(num_edges, 10); e++) {
  //     std::cout << "  curOmega(" << e << ")=" << curOmega(e) << std::endl;
  // }
  // if(num_edges > 10) {
  //     std::cout << "  ... (showing first 10 of " << num_edges << " edges)" <<
  //     std::endl;
  // }
}

void upsampleMeshAndFields() {
  secMesh = convert2SecMesh(curTriV, triF);

  std::shared_ptr<ComplexLoopZuenko> complexLoopOpt =
      std::make_shared<ComplexLoopZuenko>();
  complexLoopOpt->setBndFixFlag(isFixedBnd);
  complexLoopOpt->SetMesh(secMesh);

  MatrixX X;
  complexLoopOpt->_mesh.GetPos(X);

  std::vector<Complex> augmentedAmp;
  convertAmplitudeToComplexZ(curAmp, augmentedAmp);

  Eigen::VectorXd localOmegaNewPhase = swapEdgeVec(triF, curOmega, 0);
  // Eigen::VectorXd localOmegaNewPhase = edgeVecPhase;
  curUpZ = curZ;

  Eigen::VectorXd localOmegaNewAmps = Eigen::VectorXd::Zero(triMesh.nEdges());
  curUpAmp = augmentedAmp;

  int nverts = secMesh.GetVertCount();

  Eigen::VectorXd amp_phases(nverts);

  for (int i = 0; i < nverts; i++) {
    amp_phases(i) = std::abs(curZ[i]);
  }

  std::vector<std::complex<double>> upZvalsNewPhase;
  std::vector<std::complex<double>> upZvalsNewAmps;
  std::vector<Vector3> points;

  for (int l = 0; l < config.upsampleTimes; ++l) {
    X = S0List[l] * X;

    // complexLoopOpt->updateLoopedZvals(localOmegaNewPhase, curUpZ,
    // upZvalsNewPhase, meshList[l]);
    // complexLoopOpt->updateLoopedZvals(localOmegaNewAmps, curUpAmp,
    // upZvalsNewAmps, meshList[l]);

    complexLoopOpt->updateLoopedZvals(localOmegaNewPhase, curUpZ,
                                      upZvalsNewPhase, localOmegaNewAmps,
                                      curUpAmp, upZvalsNewAmps, meshList[l]);

    localOmegaNewPhase = S1List[l] * localOmegaNewPhase;
    localOmegaNewAmps = S1List[l] * localOmegaNewAmps;

    ConvertToVector3(X, points);
    complexLoopOpt->_mesh.Populate(points);

    curUpZ.swap(upZvalsNewPhase);
    curUpAmp.swap(upZvalsNewAmps);
  }

  parseSecMesh(complexLoopOpt->_mesh, curUpsampledTriV, upsampledTriF,
               meshList[config.upsampleTimes]);
  // upsampledTriMesh = MeshConnectivity(upsampledTriMesh);

  augmentUpsampledZList(curUpZ, curUpAmp);

  getWrinkledMesh(curUpsampledTriV, upsampledTriF, curUpZ, &vertNeiFaces,
                  curWrinkledV, config.ampScale, isUseV2);
}
bool reconstructFaceOmegas(const Eigen::MatrixXd &faceOmegaList2DinputFrame,
                           const Eigen::MatrixXd &V,
                           Eigen::MatrixXd &faceOmegaList3DinputFrame,
                           double scale = 1.0) {

  int num_faces = triF.rows();
  if (faceOmegaList2DinputFrame.rows() != num_faces) {
    std::cerr << "Error: 2D omega input rows must match face count."
              << std::endl;
    return false;
  }

  // Resize the output matrix to store the 3D vectors (Nf x 3)
  faceOmegaList3DinputFrame.resize(num_faces, 3);

  // ... (Initial setup, resizing, loop) ...
  const double GEOM_TOL_SQ = 1e-15; // Tolerance for area squared check

  for (int f_idx = 0; f_idx < num_faces; ++f_idx) {
    // ... (Code to calculate v0, v1, v2, e1, e2, c) ...

    // 1. Get vertex coordinates for face (i, j, k)
    int i = triF(f_idx, 0);
    int j = triF(f_idx, 1);
    int k = triF(f_idx, 2);

    Eigen::Vector3d v0 = V.row(i).transpose();
    Eigen::Vector3d v1 = V.row(j).transpose();
    Eigen::Vector3d v2 = V.row(k).transpose();

    // 2. Define the basis vectors (edges v0->v1 and v0->v2)
    Eigen::Vector3d e1 = v1 - v0; // Edge v0->v1
    Eigen::Vector3d e2 = v2 - v0; // Edge v0->v2

    // 3. Define the system: U * psi = c
    Eigen::MatrixXd U(2, 3);
    U.row(0) = e1;
    U.row(1) = e2;

    // c is the 2x1 vector of input components [c1; c2]
    Eigen::Vector2d c;
    c(0) = faceOmegaList2DinputFrame(f_idx, 0); // component along e1
    c(1) = faceOmegaList2DinputFrame(f_idx, 1); // component along e2

    Eigen::Vector3d psi_vec;

    // --- 1. Check for Trivial Input (c = 0) ---
    if (c.norm() < 1e-12) {
      psi_vec = Eigen::Vector3d::Zero();
    }

    // --- 2. Geometric Reconstruction (2x2 Solve) ---
    else {
      // Gram Matrix G (2x2): G_ij = e_i . e_j
      Eigen::Matrix2d G;
      G(0, 0) = e1.dot(e1);
      G(0, 1) = e1.dot(e2);
      G(1, 0) = e2.dot(e1);
      G(1, 1) = e2.dot(e2);

      // Determinant (proportional to Area squared)
      double detG = G.determinant();

      if (std::abs(detG) < GEOM_TOL_SQ) {
        // Triangle is degenerate or ill-conditioned (Area is zero/near-zero).
        psi_vec = Eigen::Vector3d::Zero();
      } else {
        // Solve G * x = c for the coefficients x = [x1; x2] in the e1, e2
        // basis.
        Eigen::Vector2d x = G.inverse() * c;

        // Reconstruct the final 3D vector: psi_vec = x1*e1 + x2*e2
        // This vector is guaranteed to lie in the plane of the triangle.
        psi_vec = x(0) * e1 + x(1) * e2;
      }
    }

    // 5. Apply optional scaling and store
    faceOmegaList3DinputFrame.row(f_idx) = (psi_vec * scale).transpose();
  }
  return true;
}

// void getEdgeOmegas(const Eigen::MatrixXd& V, const Eigen::MatrixXd&
// vertex_grad_phis, Eigen::VectorXd& edgeOmegaList)
// {
//     int num_edges = triMesh.nEdges();
//     edgeOmegaList.resize(num_edges);

//     for (int edge_id = 0; edge_id < num_edges; ++edge_id)
//     {
//         // Get the two vertices (i, j) connected by the edge
//         int i = triMesh.edgeVertex(edge_id, 0);
//         int j = triMesh.edgeVertex(edge_id, 1);

//         // 1. Edge vector
//         Eigen::Vector3d V_i = V.row(i).transpose();
//         Eigen::Vector3d V_j = V.row(j).transpose();
//         Eigen::Vector3d edge_vec = V_j - V_i; // Direction i -> j

//         // 2. Average per-vertex omega field
//         Eigen::Vector3d grad_phi_i = vertex_grad_phis.row(i).transpose();
//         Eigen::Vector3d grad_phi_j = vertex_grad_phis.row(j).transpose();

//         //CONSIDER CHANGING THIS TO LINE FIELD ADDITION!!!
//         Eigen::Vector3d avg_grad_phi = lineFieldAddition({grad_phi_i,
//         grad_phi_j});
//         //Eigen::Vector3d avg_grad_phi = (grad_phi_i + grad_phi_j)/2.0;

//         // 3. Projection (Dot Product)
//         // grad_along_edge = omega_ij
//         double grad_along_edge = avg_grad_phi.dot(edge_vec);

//         // Store the result
//         edgeOmegaList(edge_id) = grad_along_edge;
//     }

// }

// void collisionResolution() {
//     if (!curSdf) {
//         return; // No SDF loaded, skip collision resolution
//     }

//     int nverts = curWrinkledV.rows();

//     // Store parallel results: (isCollision, correctedPosition)
//     std::vector<std::pair<bool, Eigen::Vector3d>> collisionResults(nverts);

//     // Collision detection (parallel)
//     #ifdef _OPENMP
//     #pragma omp parallel for
//     #endif
//     for (int i = 0; i < nverts; i++) {
//         auto domain = curSdf->domain();
//         Eigen::Vector3d point = curWrinkledV.row(i).transpose();

//         if (domain.contains(point)) {
//             Eigen::Vector3d gradient;
//             double sdfValue = curSdf->interpolate(0, point, &gradient);

//             // Assuming threshold of 0.0 (can be made configurable)
//             const double threshold = config.sdfThreshold;
//             sdfValue = sdfValue - threshold;

//             if (sdfValue < 0) {
//                 sdfValue = -sdfValue;
//                 gradient.normalize();
//                 // Move point outside by distance (sdfValue + threshold)
//                 along gradient collisionResults[i] = std::make_pair(true,
//                 point + (sdfValue + threshold) * gradient);
//             } else {
//                 collisionResults[i] = std::make_pair(false, point);
//             }
//         } else {
//             // Point outside SDF domain, no collision
//             collisionResults[i] = std::make_pair(false, point);
//         }
//     }

//     // Update positions (sequential to avoid race conditions)
//     bool hadCollision = false;
//     for (int i = 0; i < nverts; i++) {
//         if (collisionResults[i].first) {
//             curWrinkledV.row(i) = collisionResults[i].second.transpose();
//             hadCollision = true;
//         }
//     }

//     if (hadCollision) {
//         std::cout << "  [Collision] Resolved collisions for upsampled mesh"
//         << std::endl;
//     }
// }

void calculateVertexMassesAndEdgeAreas(const Eigen::MatrixXd &V,
                                       const Eigen::MatrixXi &F,
                                       const MeshConnectivity &mesh, double rho,
                                       Eigen::VectorXd &vertexMasses,
                                       Eigen::VectorXd &edgeAreas) {

  vertexMasses.setZero(V.rows());

  int num_edges = mesh.nEdges();
  edgeAreas.setZero(num_edges);

  for (int f_idx = 0; f_idx < F.rows(); ++f_idx) {
    int i = F(f_idx, 0);
    int j = F(f_idx, 1);
    int k = F(f_idx, 2);

    const Eigen::Vector3d &vi = V.row(i);
    const Eigen::Vector3d &vj = V.row(j);
    const Eigen::Vector3d &vk = V.row(k);

    double area = 0.5 * ((vj - vi).cross(vk - vi)).norm();

    edgeAreas[mesh.faceEdge(f_idx, 0)] += area / 3.0;
    edgeAreas[mesh.faceEdge(f_idx, 1)] += area / 3.0;
    edgeAreas[mesh.faceEdge(f_idx, 2)] += area / 3.0;

    vertexMasses[i] += rho * area / 3.0;
    vertexMasses[j] += rho * area / 3.0;
    vertexMasses[k] += rho * area / 3.0;
  }
}

void calculateEdgeLengths(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                          const MeshConnectivity &mesh,
                          Eigen::VectorXd &edgeLengths) {

  edgeLengths.resize(mesh.nEdges());

  for (int e_idx = 0; e_idx < mesh.nEdges(); ++e_idx) {
    int i = mesh.edgeVertex(e_idx, 0);
    int j = mesh.edgeVertex(e_idx, 1);

    const Eigen::Vector3d &vi = V.row(i);
    const Eigen::Vector3d &vj = V.row(j);

    edgeLengths[e_idx] = (vi - vj).norm();
  }
}

bool loadSolvedProblem(const int &i, const int &uvPresent,
                       std::string &seamFile, std::string &baseMeshFile,
                       std::string &verticesFile, std::string &facesFile,
                       std::string &amplitudesFile, std::string &faceOmegasFile,
                       std::string &sdfValueFile, std::string &referenceFile) {
  std::string curVerticesFile =
      config.workingFolder + verticesFile + '.' + std::to_string(i);
  readVertices(curVerticesFile, curTriV);

  if (i == 0) {

    std::string objPathFull = config.workingFolder + baseMeshFile;

    // readObjVerticesAndFaces(objPathFull, uvPresent, baseTriV, triF, uvVerts);
    //  Eigen::MatrixXd dummyNormals;
    //  Eigen::MatrixXi dummyNormalFaces;
    //  igl::readOBJ(objPathFull, baseTriV, uvVerts, dummyNormals, triF,
    //  uvFaces, dummyNormalFaces);

    // buildUVto3DMapping(triF, uvFaces, baseTriV.rows(), uvVerts.rows(),
    // uvTo3D, threeDtoUV);
    std::string seamPath = config.workingFolder + seamFile;
    readEigenVectorFromTxt(seamPath, seamVerts);

    std::string refPath = config.workingFolder + referenceFile;
    igl::readOBJ(refPath, refTriV, refTriF);

    baseTriV = curTriV;

    facesFile = config.workingFolder + facesFile;
    readFaces(facesFile, triF);

    // baseTriV contains the 3d base mesh, 42 vertices
    // triF contains the face mapping in the 3d vertex context
    // uvVerts contain the 2d positions of vertices in uv context, 63 vertices

    triMesh = MeshConnectivity(triF); // need to review if this is to be changed
    // uvMesh = MeshConnectivity(uvFaces);
    refTriMesh = MeshConnectivity(refTriF);
    if (config.initializationMethod == "random")
      initializeRandomUnitNormZvals(baseTriV, prevZvals);
    else
      initializeFixedUnitNormZvals(baseTriV, prevZvals);

    prev2Zvals = prevZvals;

    int nverts = prevZvals.size();
    build_vertex_to_face_adjacency(triF, nverts, vert_to_faces);
    vertex_coloring_mcs(triF, nverts, V_by_color);

    calculateVertexMassesAndEdgeAreas(baseTriV, triF, triMesh,
                                      config.meshDensityRho, vertexMasses,
                                      edgeAreas);

    calculateEdgeLengths(baseTriV, triF, triMesh, edgeLengths);

    secMesh = convert2SecMesh(baseTriV, triF);

    complexLoopGlobal->setBndFixFlag(isFixedBnd);

    complexLoopGlobal->SetMesh(secMesh);

    // initializing s0 and s1
    S0List.resize(config.upsampleTimes);
    S1List.resize(config.upsampleTimes);
    meshList.resize(config.upsampleTimes + 1);

    meshList[0] = secMesh;
    std::vector<Vector3> points;

    for (int j = 0; j < config.upsampleTimes; j++) {

      complexLoopGlobal->BuildS0(S0List[j]);
      complexLoopGlobal->BuildS1(S1List[j]);

      MatrixX X;
      complexLoopGlobal->GetMesh().GetPos(X);

      X = S0List[j] * X;

      ConvertToVector3(X, points);

      std::vector<std::vector<int>> edgeToVert;
      complexLoopGlobal->GetSubdividedEdges(edgeToVert);

      std::vector<std::vector<int>> faceToVert;
      complexLoopGlobal->GetSubdividedFaces(faceToVert);

      complexLoopGlobal->_mesh.Populate(points, faceToVert, edgeToVert);

      meshList[j + 1] = complexLoopGlobal->GetMesh();
    }

    buildVertexNeighboringFaces(upsampledTriMesh, points.size(), vertNeiFaces);
  }

  std::string curAmplitudesFile =
      config.workingFolder + amplitudesFile + '.' + std::to_string(i);
  loadVertexAmp(curAmplitudesFile, curTriV.rows(), curAmp);

  // map3DVertexDataToUV(curAmp, uvTo3D, uvCurAmp);

  // read face omega list
  std::string curFaceOmegasFile =
      config.workingFolder + faceOmegasFile + '.' + std::to_string(i);
  loadFaceOmegas(curFaceOmegasFile, triF.rows(), curFaceOmegaList2Dinput);

  return true;
}

// void fixCurOmega() {
//     for(int e_idx = 0; e_idx < uvMesh.nEdges(); ++e_idx) {

//         double w_ij = curOmega[e_idx];

//         Complex exp_ij = std::exp(Complex(0.0, w_ij));
//         Complex exp_ji = std::conj(exp_ij);

//         int i = uvMesh.edgeVertex(e_idx, 0);
//         int j = uvMesh.edgeVertex(e_idx, 1);

//         double error1 = std::abs(curZ[i] - exp_ij * curZ[j]);
//         double error2 = std::abs(curZ[i] - exp_ji * curZ[j]); //correct
//         forward energy

//         double edash1 = std::abs(curOmega[e_idx] - prevOmega[e_idx]);
//         double edash2 = std::abs(curOmega[e_idx] + prevOmega[e_idx]);

//         if(std::abs(curZ[i]-curZ[j]) < std::abs(curZ[i]*std::exp(Complex(0.0,
//         M_PI*0.005)))) {
//             if(edash2 < edash1) curOmega[e_idx] = -curOmega[e_idx];
//         }
//         else{
//             if(error2 > error1) curOmega[e_idx] = -curOmega[e_idx];
//         }

//     }
// }

// void resolveSeams(const std::vector<std::vector<int>>& threeDtoUV,
// std::vector<std::complex<double>>& psi_current_complex){
//     for(const auto& x: threeDtoUV) {
//         if(x.empty()) continue;

//         std::complex<double> prod_z(1.0, 0.0);
//         int count = 0;

//         for(int y : x) {
//             if(y >= 0 && y < psi_current_complex.size()) {
//                 prod_z *= psi_current_complex[y];
//                 count++;
//             }
//         }

//         if (count == 0) continue;

//         // Geometric mean: (z1 * z2 * ... * zn)^(1/n)
//         // Since we want to average phases, geometric meam is appropriate for
//         unit complex numbers.
//         // However, we should be careful about branch cuts.
//         // A simpler way for phases might be summing normalized vectors if
//         amplitudes are ignored,
//         // but sticking to geometric mean as per original logic.

//         std::complex<double> avg_z = std::pow(prod_z, 1.0 / count);
//         // Ensure unit norm if the original logic intended just phase
//         averaging
//         // avg_z /= std::abs(avg_z);

//         for(int y : x) {
//             if(y >= 0 && y < psi_current_complex.size()) {
//                  psi_current_complex[y] = avg_z;
//             }
//         }
//     }
// }

// void resolveVertexGradPhis(const std::vector<std::vector<int>>& threeDtoUV,
// Eigen::MatrixXd& curVertexOmega) {
//     for(const auto& x: threeDtoUV) {
//         if(x.empty()) continue;
//         std::vector<Eigen::Vector3d> grad_phi_list;
//         for(int y : x) {
//              if(y >= 0 && y < curVertexOmega.rows())
//                 grad_phi_list.emplace_back(curVertexOmega.row(y).transpose());
//         }

//         if(grad_phi_list.empty()) continue;

//         Eigen::Vector3d avg_grad_phi = lineFieldAddition(grad_phi_list);

//         for(int y : x) {
//              if(y >= 0 && y < curVertexOmega.rows())
//                 curVertexOmega.row(y) = avg_grad_phi.transpose();
//         }
//     }
// }

// Smallest rotation taking the (flat) mesh's plane to the XZ plane, which is
// where lineFieldAddition works. The plane normal is oriented towards +Y, so
// the rotation is exactly the identity for meshes already in XZ.
Eigen::Matrix3d computeMinimalPlaneToXZRotation(const Eigen::MatrixXd &V,
                                                const Eigen::MatrixXi &F) {
  Eigen::Vector3d n = Eigen::Vector3d::Zero();
  for (int f = 0; f < F.rows(); ++f) {
    Eigen::Vector3d v0 = V.row(F(f, 0));
    Eigen::Vector3d v1 = V.row(F(f, 1));
    Eigen::Vector3d v2 = V.row(F(f, 2));
    n += (v1 - v0).cross(v2 - v0);
  }
  if (n.norm() < 1e-12) {
    throw std::runtime_error(
        "computeMinimalPlaneToXZRotation: degenerate reference mesh normal.");
  }
  n.normalize();
  if (n[1] < 0)
    n = -n;
  const Eigen::Vector3d up(0, 1, 0);
  if (n == up)
    return Eigen::Matrix3d::Identity();
  return Eigen::Quaterniond::FromTwoVectors(n, up).toRotationMatrix();
}

void getEdgeOmegasInReferenceMesh(const Eigen::MatrixXd &V,
                                  const Eigen::MatrixXi &F,
                                  const Eigen::MatrixXd &vertex_grad_phis,
                                  const MeshConnectivity &tmesh,
                                  Eigen::VectorXd &edgeOmegaList) {
  int num_edges = tmesh.nEdges();
  edgeOmegaList.resize(num_edges);
  edgeGradPhis.resize(num_edges, 3);
  // lineFieldAddition only sees x and z, so average in the XZ plane and
  // rotate back (reference meshes in e.g. the XY plane otherwise get every
  // edge gradient along x)
  const Eigen::Matrix3d R = computeMinimalPlaneToXZRotation(V, F);
  for (int edge_id = 0; edge_id < num_edges; ++edge_id) {
    // Get the two vertices (i, j) connected by the edge
    int i = tmesh.edgeVertex(edge_id, 0);
    int j = tmesh.edgeVertex(edge_id, 1);

    // 1. Edge vector
    Eigen::Vector3d V_i = V.row(i).transpose();
    Eigen::Vector3d V_j = V.row(j).transpose();
    Eigen::Vector3d edge_vec = V_j - V_i; // Direction i -> j

    // 2. Average per-vertex omega field
    Eigen::Vector3d grad_phi_i = vertex_grad_phis.row(i).transpose();
    Eigen::Vector3d grad_phi_j = vertex_grad_phis.row(j).transpose();

    // CONSIDER CHANGING THIS TO LINE FIELD ADDITION!!!
    Eigen::Vector3d avg_grad_phi =
        R.transpose() * lineFieldAddition({R * grad_phi_i, R * grad_phi_j});
    // Eigen::Vector3d avg_grad_phi = (grad_phi_i + grad_phi_j)/2.0;
    edgeGradPhis.row(edge_id) = avg_grad_phi;
    // 3. Projection (Dot Product)
    // grad_along_edge = omega_ij
    double grad_along_edge = avg_grad_phi.dot(edge_vec);

    // Store the result
    edgeOmegaList(edge_id) = grad_along_edge;
  }
}
bool isSeam(int vid) {
  auto it =
      std::find(seamVerts.data(), seamVerts.data() + seamVerts.size(), vid);
  if (it != seamVerts.data() + seamVerts.size()) {
    return true;
  }
  return false;
}
void getEdgeOmegasNew(Eigen::MatrixXd &faceGradPhiMatrix) {
  curOmega.resize(triMesh.nEdges());
  int nfaces = refTriF.rows();
  Eigen::Matrix<double, Eigen::Dynamic, 3> faceGradPhis;
  faceGradPhis.resize(nfaces, 3);
  // construct face grad phis
  for (int faceid = 0; faceid < nfaces; faceid++) {
    Eigen::Vector<double, 2> dphi = curFaceOmegaList2Dinput.row(faceid);
    Eigen::Matrix<double, 3, 2> A;
    Eigen::Vector3d vi = refTriV.row(refTriF(faceid, 0));
    Eigen::Vector3d vj = refTriV.row(refTriF(faceid, 1));
    Eigen::Vector3d vk = refTriV.row(refTriF(faceid, 2));
    A << vj - vi, vk - vi;
    faceGradPhis.row(faceid) = A * (A.transpose() * A).inverse() * dphi;
  }
  faceGradPhiMatrix = faceGradPhis;
  // now compute line field interpolation result

  Eigen::VectorXd edgeOmegas;
  getGradPhisPerVertexLineField(faceGradPhis, refTriV, refTriF, vertexGradPhis);
  getEdgeOmegasInReferenceMesh(refTriV, refTriF, vertexGradPhis, refTriMesh,
                               edgeOmegas);
  // map edge omegas in the reference mesh to base mesh
  // finally map edge omegas from reference mesh to the curOmega vector using
  for (int faceid = 0; faceid < nfaces; faceid++) {
    int vi = refTriF(faceid, 0);
    int vj = refTriF(faceid, 1);
    int vk = refTriF(faceid, 2);
    bool i_seam, j_seam, k_seam;
    i_seam = isSeam(vi);
    j_seam = isSeam(vj);
    k_seam = isSeam(vk);
    // edge opposite to vi vj-vk
    if (!(j_seam || k_seam)) {
      // fill currOmega
      int edgeid = triMesh.faceEdge(faceid, 0);
      int refEdgeId = refTriMesh.faceEdge(faceid, 0);
      curOmega(edgeid) = edgeOmegas(refEdgeId);
    }

    // edge opposite to vj vk-vi
    if (!(i_seam || k_seam)) {
      // fill currOmega
      int edgeid = triMesh.faceEdge(faceid, 1);
      int refEdgeId = refTriMesh.faceEdge(faceid, 1);
      curOmega(edgeid) = edgeOmegas(refEdgeId);
    }

    // edge opposite to vk vi-vj
    if (!(i_seam || j_seam)) {
      // fill currOmega
      int edgeid = triMesh.faceEdge(faceid, 2);
      int refEdgeId = refTriMesh.faceEdge(faceid, 2);
      curOmega(edgeid) = edgeOmegas(refEdgeId);
    }
  }
}

/* // get edge omegas using uv mesh approach
void getEdgeOmegas() {
  // we need to fill curOmega vector
  // first compute dphis per face in the reference mesh
  // refTriV, refTriF, refTriMesh
  // for the base mesh, we have curFaceOmegaList2Dinput

  // then get line field interpolation to get the edge gradphi

  // then compute edge omegas for the reference mesh
  // finally map edge omegas from reference mesh to the curOmega vector using
  // seam information
} */

int main(int argc, char **argv) {
  CLI::App app("Wrinkle Interpolation");
  app.add_option("input,-i,--input", args.paramsFile, "Input model")
      ->required()
      ->check(CLI::ExistingFile);
  app.add_option("checkpoint,-c,--c", args.saveCheckpoint, "save checkpoint?")
      ->required();
  app.add_flag("--headless", args.headless,
               "run without the polyscope viewer (no display needed)");
  app.add_option("--save-every", args.saveEvery,
                 "save every Nth frame; the last frame is always saved")
      ->check(CLI::PositiveNumber);

  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError &e) {
    return app.exit(e);
  }

  loadRenderConfig(args.paramsFile);
  std::string renderFolder = config.workingFolder + config.renderFolder;
  std::string outputFolderWrinkles = renderFolder + "/wrinkledMesh/";
  std::string outputFolderSim = renderFolder + "/simMesh/";
  std::string outputFolderPhases = renderFolder + "/phases/";
  std::string outputFolderUpsampledPhases = renderFolder + "/upsampledPhases/";
  std::string outputFolderVertexGradPhis = renderFolder + "/vertexGradPhis/";
  std::string outputFolderFaceGradPhis = renderFolder + "/faceGradPhis/";
  std::string outputFolderEdgeGradPhis = renderFolder + "/edgeGradPhis/";
  std::string outputFolderEdgeOmegas = renderFolder + "/edgeOmegas/";
  std::string outputFolderAmplitudes = renderFolder + "/amplitudes/";
  std::string outputFolderUpsampledAmplitudes =
      renderFolder + "/upsampledAmplitudes/";

  if (args.saveCheckpoint) {
    mkdir(renderFolder);
    mkdir(outputFolderWrinkles);
    mkdir(outputFolderSim);
    mkdir(outputFolderPhases);
    mkdir(outputFolderUpsampledPhases);
    mkdir(outputFolderVertexGradPhis);
    mkdir(outputFolderFaceGradPhis);
    mkdir(outputFolderEdgeGradPhis);
    mkdir(outputFolderEdgeOmegas);
    mkdir(outputFolderAmplitudes);
    mkdir(outputFolderUpsampledAmplitudes);
  }

  auto start = std::chrono::high_resolution_clock::now();
  auto end = std::chrono::high_resolution_clock::now();
  long long loading_time, optimization_time, updating_time, saving_time;
  VizCallback visualizer(config.numFrames);
  if (!args.headless) {
    polyscope::init();
    polyscope::options::groundPlaneMode =
        polyscope::GroundPlaneMode::ShadowOnly;
    polyscope::state::userCallback = [&]() { visualizer.callback(); };
  }
  for (int i = 0; i < config.numFrames; ++i) {

    std::cout << "--- Processing Frame " << i << " ---" << std::endl;
    auto frame_start = std::chrono::high_resolution_clock::now();

    start = std::chrono::high_resolution_clock::now();

    if (!loadSolvedProblem(
            i, config.uvPresent, config.seamFile, config.baseMeshFile,
            config.verticesFile, config.facesFile, config.amplitudesFile,
            config.faceOmegasFile, config.sdfValueFile, config.referenceFile)) {
      std::cout << "failed to load file." << std::endl;
      return 1;
    }

    end = std::chrono::high_resolution_clock::now();
    loading_time =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();

    start = std::chrono::high_resolution_clock::now();
    // get face grad phis
    // TODO: add get edge omegas here
    // reconstructFaceOmegas(curFaceOmegaList2Dinput, baseTriV,
    //                       curFaceOmegaList3Dinput);
    // getEdgeOmegas();
    getEdgeOmegasNew(curFaceOmegaList3Dinput);
    // edge omegas from line field interpolation in the reference mesh
    // uses edge flap method
    // faceToEdgeOmegas(curFaceOmegaList2Dinput, baseTriV, curOmega);

    if (config.optimizationMethod == "bfgs") {
      LBFGSpp::LBFGSParam<double> param;
      param.max_iterations = config.intermediateIters;
      param.epsilon = 1e-5;
      param.delta = 1e-5;
      optimizeAndSaveZvalsBFGS(curFaceOmegaList2Dinput, curZ, curAmp, prevZvals,
                               triF, config, param, i);
    } else if (config.optimizationMethod == "vbd") {

      optimizeAndSaveZvalsVBDObj(curOmega, curZ, curAmp, vertexMasses,
                                 edgeLengths, edgeAreas, prevZvals, prev2Zvals,
                                 triMesh, config, i, vert_to_faces, V_by_color);

    } else {
      optimizeAndSaveZvalsGDObj(curOmega, curZ, curAmp, prevZvals, triMesh,
                                config, i);
    }
    // fix orientation
    // TODO: NOt needed in the uv mesh based appraoch
    // fixOmegas(baseTriV, curOmega, curZ);

    if (i == 0) {
      prevOmega = curOmega;
    }

    end = std::chrono::high_resolution_clock::now();
    optimization_time =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();

    start = std::chrono::high_resolution_clock::now();

    upsampleMeshAndFields();
    // upsampledTriMesh = MeshConnectivity(upsampledTriMesh);

    end = std::chrono::high_resolution_clock::now();
    updating_time =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    // visualize
    if (!args.headless)
      visualizer.updateMesh(curWrinkledV, upsampledTriF);

    // ---------------- SAVING ----------------
    const bool saveFrame =
        args.saveCheckpoint &&
        (i % args.saveEvery == 0 || i == config.numFrames - 1);
    if (saveFrame) {
      start = std::chrono::high_resolution_clock::now();

      // 1. Save Upsampled Wrinkled Mesh
      igl::writeOBJ(outputFolderWrinkles + "wrinkledMesh_" + std::to_string(i) +
                        ".obj",
                    curWrinkledV, upsampledTriF);

      // 2. Save Coarse Simulation Mesh (Geometry for phases)
      igl::writeOBJ(outputFolderSim + "simMesh_" + std::to_string(i) + ".obj",
                    curTriV, triF);

      // 3. Save Coarse Phases (Data on simMesh)
      std::string phaseFilepath =
          outputFolderPhases + "phases_" + std::to_string(i) + ".txt";
      std::ofstream phaseFile(phaseFilepath);

      if (phaseFile.is_open()) {
        phaseFile << std::scientific << std::setprecision(16);
        // Iterate over curZ (Coarse)
        for (const auto &val : curZ) {
          phaseFile << std::arg(val) << "\n";
        }
        phaseFile.close();
      } else {
        std::cerr << "Error: Unable to save phases to " << phaseFilepath
                  << std::endl;
      }

      // 4. Save Upsampled Phases (Data on upsampled mesh)
      std::string upsampledPhaseFilepath = outputFolderUpsampledPhases +
                                           "upsampledPhases_" +
                                           std::to_string(i) + ".txt";
      std::ofstream upsampledPhaseFile(upsampledPhaseFilepath);

      if (upsampledPhaseFile.is_open()) {
        upsampledPhaseFile << std::scientific << std::setprecision(16);
        for (const auto &val : curUpZ) {
          upsampledPhaseFile << std::arg(val) << "\n";
        }
        upsampledPhaseFile.close();
      } else {
        std::cerr << "Error: Unable to save upsampled phases to "
                  << upsampledPhaseFilepath << std::endl;
      }

      // 6. Save Face Grad Phis
      std::string faceGradPhisFilepath = outputFolderFaceGradPhis +
                                         "faceGradPhis_" + std::to_string(i) +
                                         ".txt";
      std::ofstream faceGradPhisFile(faceGradPhisFilepath);

      if (faceGradPhisFile.is_open()) {
        faceGradPhisFile << std::scientific << std::setprecision(16);
        for (int f = 0; f < curFaceOmegaList3Dinput.rows(); ++f) {
          faceGradPhisFile << curFaceOmegaList3Dinput(f, 0) << " "
                           << curFaceOmegaList3Dinput(f, 1) << " "
                           << curFaceOmegaList3Dinput(f, 2) << "\n";
        }
        faceGradPhisFile.close();
      } else {
        std::cerr << "Error: Unable to save face grad phis to "
                  << faceGradPhisFilepath << std::endl;
      }

      // 6. Save Edge Grad Phis
      std::string edgeGradPhisFilepath = outputFolderEdgeGradPhis +
                                         "edgeGradPhis_" + std::to_string(i) +
                                         ".txt";
      std::ofstream edgeGradPhisFile(edgeGradPhisFilepath);

      if (edgeGradPhisFile.is_open()) {
        edgeGradPhisFile << std::scientific << std::setprecision(16);
        for (int f = 0; f < edgeGradPhis.rows(); ++f) {
          edgeGradPhisFile << edgeGradPhis(f, 0) << " " << edgeGradPhis(f, 1)
                           << " " << edgeGradPhis(f, 2) << "\n";
        }
        edgeGradPhisFile.close();
      } else {
        std::cerr << "Error: Unable to save edge grad phis to "
                  << edgeGradPhisFilepath << std::endl;
      }

      // 6. Save Vetex Grad Phis
      std::string vertexGradPhisFilepath = outputFolderVertexGradPhis +
                                           "vertexGradPhis_" +
                                           std::to_string(i) + ".txt";
      std::ofstream vertexGradPhisFile(vertexGradPhisFilepath);

      if (vertexGradPhisFile.is_open()) {
        vertexGradPhisFile << std::scientific << std::setprecision(16);
        for (int f = 0; f < vertexGradPhis.rows(); ++f) {
          vertexGradPhisFile << vertexGradPhis(f, 0) << " "
                             << vertexGradPhis(f, 1) << " "
                             << vertexGradPhis(f, 2) << "\n";
        }
        vertexGradPhisFile.close();
      } else {
        std::cerr << "Error: Unable to save vertex grad phis to "
                  << vertexGradPhisFilepath << std::endl;
      }

      // 7. Save Edge Omegas
      std::string edgeOmegasFilepath =
          outputFolderEdgeOmegas + "edgeOmegas_" + std::to_string(i) + ".txt";
      std::ofstream edgeOmegasFile(edgeOmegasFilepath);

      if (edgeOmegasFile.is_open()) {
        edgeOmegasFile << std::scientific << std::setprecision(16);
        for (int e = 0; e < curOmega.size(); ++e) {
          edgeOmegasFile << curOmega(e) << "\n";
        }
        edgeOmegasFile.close();
      } else {
        std::cerr << "Error: Unable to save edge omegas to "
                  << edgeOmegasFilepath << std::endl;
      }

      // 8. Save Coarse Amplitudes (Data on simMesh vertices)
      std::string amplitudesFilepath =
          outputFolderAmplitudes + "amplitudes_" + std::to_string(i) + ".txt";
      std::ofstream amplitudesFile(amplitudesFilepath);

      if (amplitudesFile.is_open()) {
        amplitudesFile << std::scientific << std::setprecision(16);
        for (int v = 0; v < curAmp.size(); ++v) {
          amplitudesFile << std::abs(curAmp[v]) << "\n";
        }
        amplitudesFile.close();
      } else {
        std::cerr << "Error: Unable to save amplitudes to "
                  << amplitudesFilepath << std::endl;
      }

      // 9. Save Upsampled Amplitudes (Data on upsampled mesh vertices)
      std::string upsampledAmplitudesFilepath =
          outputFolderUpsampledAmplitudes + "upsampledAmplitudes_" +
          std::to_string(i) + ".txt";
      std::ofstream upsampledAmplitudesFile(upsampledAmplitudesFilepath);

      if (upsampledAmplitudesFile.is_open()) {
        upsampledAmplitudesFile << std::scientific << std::setprecision(16);
        for (int v = 0; v < curUpAmp.size(); ++v) {
          upsampledAmplitudesFile << std::abs(curUpAmp[v]) << "\n";
        }
        upsampledAmplitudesFile.close();
      } else {
        std::cerr << "Error: Unable to save upsampled amplitudes to "
                  << upsampledAmplitudesFilepath << std::endl;
      }

      end = std::chrono::high_resolution_clock::now();
      saving_time =
          std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
              .count();
    }

    prev2Zvals = prevZvals;
    prevZvals = curZ;

    // prevAmp = uvCurAmp;
    prevAmp = curAmp;

    prevOmega = curOmega;

    auto frame_end = std::chrono::high_resolution_clock::now();
    long long total_frame_time =
        std::chrono::duration_cast<std::chrono::milliseconds>(frame_end -
                                                              frame_start)
            .count();

    std::cout << "  Loading time:      " << loading_time << " ms" << std::endl;
    std::cout << "  Optimization time: " << optimization_time << " ms"
              << std::endl;
    std::cout << "  Updating time:     " << updating_time << " ms" << std::endl;
    if (saveFrame) {
      std::cout << "  Saving time:       " << saving_time << " ms" << std::endl;
    }
    std::cout << "  Total Frame time:  " << total_frame_time << " ms"
              << std::endl;
    std::cout << "-------------------------" << std::endl;
  }
  if (!args.headless)
    polyscope::show();

  // Cleanup before exit
  delete curSdf;
  curSdf = nullptr;

  return 0;
}