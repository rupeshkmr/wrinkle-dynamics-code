#pragma once
#include "graphColoringHelpers.hpp"

// --- Helper: Compute Target H using Face Data (Standard Signature) ---
// Complex computeOptimalVertexTargetStandard(
//     int v_idx,
//     const Eigen::VectorXcd& psi,
//     const MeshConnectivity& mesh,
//     const Eigen::MatrixXd& face_grad_phis_2d,
//     const Eigen::VectorXd& wa,
//     const Eigen::VectorXd& edgeLengths,
//     const Eigen::VectorXd& edgeAreas,
//     const std::vector<std::vector<int>>& vert_to_faces
// ) {
//     using std::exp;
//     using std::norm;

//     Complex H(0.0, 0.0); // The accumulator for the target vector
//     int nverts = psi.size();
//     const Eigen::MatrixXi& F = mesh.faces();

//     for (int f_idx : vert_to_faces[v_idx]) {
//         int i = F(f_idx, 0);
//         int j = F(f_idx, 1);
//         int k = F(f_idx, 2);

//         if (i < 0 || i >= nverts || j < 0 || j >= nverts || k < 0 || k >=
//         nverts)
//             continue;

//         // --- Retrieve Edge Info ---
//         int e_ij = mesh.faceEdge(f_idx, 2);
//         int e_jk = mesh.faceEdge(f_idx, 0);
//         int e_ki = mesh.faceEdge(f_idx, 1);

//         double len_ij = edgeLengths[mesh.faceEdge(f_idx, 2)];
//         double len_jk = edgeLengths[mesh.faceEdge(f_idx, 0)];
//         double len_ki = edgeLengths[mesh.faceEdge(f_idx, 1)];

//         double area_ij = edgeAreas[mesh.faceEdge(f_idx, 2)];
//         double area_jk = edgeAreas[mesh.faceEdge(f_idx, 0)];
//         double area_ki = edgeAreas[mesh.faceEdge(f_idx, 1)];

//         double ori_ij = (mesh.edgeVertex(e_ij, 0) == i) ? 1.0 : -1.0;
//         double ori_jk = (mesh.edgeVertex(e_jk, 0) == j) ? 1.0 : -1.0;
//         double ori_ki = (mesh.edgeVertex(e_ki, 0) == k) ? 1.0 : -1.0;

//         // --- Derive Edge Omegas from 2D Face Omegas ---
//         // This matches the logic in vertexComplexGrad
//         double c1 = face_grad_phis_2d(f_idx, 0);
//         double c2 = face_grad_phis_2d(f_idx, 1);

//         double w_ij = c1;
//         double w_ki = -c2;
//         double w_jk = c2 - c1;

//         // Weights
//         double s_ij = 0.5 * (wa[i] + wa[j]);
//         double s_jk = 0.5 * (wa[j] + wa[k]);
//         double s_ki = 0.5 * (wa[k] + wa[i]);

//         // Squared weights for Least Squares target
//         double w_ij_sq = s_ij * s_ij;
//         double w_jk_sq = s_jk * s_jk;
//         double w_ki_sq = s_ki * s_ki;

//         // Rotators
//         Complex rot_ij = std::exp(Complex(0.0, w_ij));
//         Complex rot_jk = std::exp(Complex(0.0, w_jk));
//         Complex rot_ki = std::exp(Complex(0.0, w_ki));

//         Complex rot_ij_inv = std::conj(rot_ij);
//         Complex rot_jk_inv = std::conj(rot_jk);
//         Complex rot_ki_inv = std::conj(rot_ki);

//         if (v_idx == i) {
//             H += area_ij * w_ij_sq * (psi(j) * rot_ij_inv) / (len_ij *
//             len_ij);  // Fwd H += area_ij * w_ij_sq * (psi(j) * rot_ij) /
//             (len_ij * len_ij);  // Fwd

//             H += area_ki * w_ki_sq * (psi(k) * rot_ki_inv) / (len_ki *
//             len_ki);      // Fwd H += area_ki * w_ki_sq * (psi(k) * rot_ki) /
//             (len_ki * len_ki);      // Fwd

//         }
//         else if (v_idx == j) {
//             H += area_jk * w_jk_sq * (psi(k) * rot_jk_inv) / (len_jk *
//             len_jk); H += area_ij * w_ij_sq * (psi(i) * rot_ij) / (len_ij *
//             len_ij); H += area_jk * w_jk_sq * (psi(k) * rot_jk) / (len_jk *
//             len_jk); H += area_ij * w_ij_sq * (psi(i) * rot_ij_inv) / (len_ij
//             * len_ij);
//         }
//         else if (v_idx == k) {
//             H += area_ki * w_ki_sq * (psi(i) * rot_ki_inv) / (len_ki *
//             len_ki); H += area_jk * w_jk_sq * (psi(j) * rot_jk) / (len_jk *
//             len_jk); H += area_ki * w_ki_sq * (psi(i) * rot_ki) / (len_ki *
//             len_ki); H += area_jk * w_jk_sq * (psi(j) * rot_jk_inv) / (len_jk
//             * len_jk);
//         }
//     }

//     return H;
// }

// void optimizeAndSaveZvalsVBD(
//     Eigen::MatrixXd& faceOmega2D,
//     std::vector<Complex>& curZ,
//     Eigen::VectorXd& wrinkleAmps,
//     Eigen::VectorXd& vertexMasses,
//     Eigen::VectorXd& edgeLengths,
//     Eigen::VectorXd& edgeAreas,
//     const std::vector<Complex>& prevZvals,
//     const std::vector<Complex>& prev2Zvals,
//     Eigen::VectorXd& prevAmps,
//     const MeshConnectivity& mesh,
//     const RenderConfig& config,
//     const int& i,
//     std::vector<std::vector<int>>& vert_to_faces,
//     std::vector<std::vector<int>>& V_by_color
// )
// {
//     int max_outer_iters = (i == 0) ? config.initialIters :
//     config.subsequentIters; const double ENERGY_TOL = 1e-6;

//     copyComplexStdVectorToEigen(prevZvals, psi_current_complex);
//     Eigen::VectorXcd dummy_grad;

//     // --- Initial Energy Log ---
//     // Note: Using the standard 'computeEnergyAndGradient' that takes
//     faceOmega2D
//     //double initial_energy = computeEnergyAndGradient(psi_current_complex,
//     triF, faceOmega2D, wrinkleAmps, dummy_grad);
//     //std::cout << "Iter " << -1 << " | Energy: " << initial_energy <<
//     std::endl;

//     for (int outer_iter = 0; outer_iter < max_outer_iters; ++outer_iter)
//     {
//         Eigen::VectorXcd psi_before_pass = psi_current_complex;

//         for (size_t color = 0; color < V_by_color.size(); ++color)
//         {
//             const std::vector<int>& vertex_batch = V_by_color[color];
//             if (vertex_batch.empty()) continue;

//             tbb::parallel_for
//             (
//                 tbb::blocked_range<int>(0, vertex_batch.size()),
//                 [&](const tbb::blocked_range<int>& r)
//                 {
//                     for (int v_local_idx = r.begin(); v_local_idx < r.end();
//                     ++v_local_idx)
//                     {
//                         int v_idx = vertex_batch[v_local_idx];

//                         // 1. Calculate Target H using the Standard helper
//                         Complex H =
//                         computeOptimalVertexTargetStandard(
//                             v_idx,
//                             psi_current_complex,
//                             mesh,
//                             faceOmega2D,
//                             wrinkleAmps,
//                             edgeLengths,
//                             edgeAreas,
//                             vert_to_faces) / 16.0;

//                         H *= config.meshThickness * config.beta;

//                         // DECIDE IF WE HAVE TO KEEP THE FACTOR OF 0.5!!!!!
//                         Complex Q = 0.5 * vertexMasses[v_idx] *
//                         wrinkleAmps[v_idx] * wrinkleAmps[v_idx] *
//                         (prevZvals[v_idx])
//                         / (config.timeStep * config.timeStep);

//                         double kdamp = 0.3;
//                         if(i>1 && prevAmps[v_idx] > config.meshThickness) {
//                             Q *= std::pow(prevZvals[v_idx] /
//                             prev2Zvals[v_idx], 1 - kdamp);
//                         }

//                         H += Q;

//                         // 2. Exact "Newton" Step: Projection to Unit Circle
//                         if (std::abs(H) > 1e-12) {
//                             psi_current_complex(v_idx) = H / std::abs(H);
//                         }
//                     }
//                 }
//             );
//         }

//         double change = (psi_current_complex - psi_before_pass).norm();

//         if (change < ENERGY_TOL && outer_iter > 0) {
//             std::cout << "VBD converged after " << outer_iter + 1 << " outer
//             iterations." << std::endl; break;
//         }
//     }

//     copyComplexEigenToStdVector(psi_current_complex, curZ);
//     //std::cout << "Frame " << i << " Z-values successfully optimized (VBD
//     Analytic)." << std::endl;
// }

Complex computeOptimalVertexTargetObj(
    int v_idx, const Eigen::VectorXcd &psi, const MeshConnectivity &mesh,
    const Eigen::VectorXd &edgeOmegas, const Eigen::VectorXd &wa,
    const Eigen::VectorXd &edgeLengths, const Eigen::VectorXd &edgeAreas,
    const std::vector<std::vector<int>> &vert_to_faces) {
  using std::exp;
  using std::norm;

  Complex H(0.0, 0.0); // The accumulator for the target vector
  int nverts = psi.size();
  const Eigen::MatrixXi &F = mesh.faces();

  for (int f_idx : vert_to_faces[v_idx]) {
    int i = F(f_idx, 0);
    int j = F(f_idx, 1);
    int k = F(f_idx, 2);

    if (i < 0 || i >= nverts || j < 0 || j >= nverts || k < 0 || k >= nverts)
      continue;

    // --- Retrieve Edge Info ---
    int e_ij = mesh.faceEdge(f_idx, 2);
    int e_jk = mesh.faceEdge(f_idx, 0);
    int e_ki = mesh.faceEdge(f_idx, 1);

    double len_ij = edgeLengths[mesh.faceEdge(f_idx, 2)];
    double len_jk = edgeLengths[mesh.faceEdge(f_idx, 0)];
    double len_ki = edgeLengths[mesh.faceEdge(f_idx, 1)];

    double area_ij = edgeAreas[mesh.faceEdge(f_idx, 2)];
    double area_jk = edgeAreas[mesh.faceEdge(f_idx, 0)];
    double area_ki = edgeAreas[mesh.faceEdge(f_idx, 1)];

    double ori_ij = (mesh.edgeVertex(e_ij, 0) == i) ? 1.0 : -1.0;
    double ori_jk = (mesh.edgeVertex(e_jk, 0) == j) ? 1.0 : -1.0;
    double ori_ki = (mesh.edgeVertex(e_ki, 0) == k) ? 1.0 : -1.0;

    double w_ij = edgeOmegas(e_ij) * ori_ij;
    double w_jk = edgeOmegas(e_jk) * ori_jk;
    double w_ki = edgeOmegas(e_ki) * ori_ki;

    // Weights (squared because energy is |s * residual|^2)
    double s_ij = 0.5 * (wa[i] + wa[j]);
    double s_jk = 0.5 * (wa[j] + wa[k]);
    double s_ki = 0.5 * (wa[k] + wa[i]);

    double w_ij_sq = s_ij * s_ij;
    double w_jk_sq = s_jk * s_jk;
    double w_ki_sq = s_ki * s_ki;

    // Precompute complex rotators
    // For standard direction: exp(i * w)
    Complex rot_ij = std::exp(Complex(0.0, w_ij));
    Complex rot_jk = std::exp(Complex(0.0, w_jk));
    Complex rot_ki = std::exp(Complex(0.0, w_ki));

    // For reverse direction: exp(-i * w) (conjugate)
    Complex rot_ij_inv = std::conj(rot_ij);
    Complex rot_jk_inv = std::conj(rot_jk);
    Complex rot_ki_inv = std::conj(rot_ki);

    if (v_idx == i) {
      H += area_ij * w_ij_sq * (psi(j) * rot_ij_inv) / (len_ij * len_ij); // Fwd
      // H += area_ij * w_ij_sq * (psi(j) * rot_ij) / (len_ij * len_ij);     //
      // Fwd

      H += area_ki * w_ki_sq * (psi(k) * rot_ki) / (len_ki * len_ki); // back
      // H +=
      //     area_ki * w_ki_sq * (psi(k) * rot_ki_inv) / (len_ki * len_ki); //
      //     back

    } else if (v_idx == j) {
      H += area_jk * w_jk_sq * (psi(k) * rot_jk_inv) / (len_jk * len_jk);
      // H += area_jk * w_jk_sq * (psi(k) * rot_jk) / (len_jk * len_jk);

      H += area_ij * w_ij_sq * (psi(i) * rot_ij) / (len_ij * len_ij);
      // H += area_ij * w_ij_sq * (psi(i) * rot_ij_inv) / (len_ij * len_ij);

    } else if (v_idx == k) {

      H += area_ki * w_ki_sq * (psi(i) * rot_ki_inv) / (len_ki * len_ki);
      // H += area_ki * w_ki_sq * (psi(i) * rot_ki) / (len_ki * len_ki);

      H += area_jk * w_jk_sq * (psi(j) * rot_jk) / (len_jk * len_jk);
      // H += area_jk * w_jk_sq * (psi(j) * rot_jk_inv) / (len_jk * len_jk);
    }

    // if (v_idx == i) {
    //     H += w_ij_sq * (psi(j) * rot_ij_inv);  // Fwd
    //     H += w_ki_sq * (psi(k) * rot_ki);      // Fwd
    // }
    // else if (v_idx == j) {
    //     H += w_jk_sq * (psi(k) * rot_jk_inv);
    //     H += w_ij_sq * (psi(i) * rot_ij);
    // }
    // else if (v_idx == k) {
    //     H += w_ki_sq * (psi(i) * rot_ki_inv);
    //     H += w_jk_sq * (psi(j) * rot_jk);
    // }
  }

  return H;
}

void optimizeAndSaveZvalsVBDObj(
    const Eigen::VectorXd &edgeOmegas, std::vector<Complex> &curZ,
    Eigen::VectorXd &wrinkleAmps, Eigen::VectorXd &vertexMasses,
    Eigen::VectorXd &edgeLengths, Eigen::VectorXd &edgeAreas,
    const std::vector<Complex> &prevZvals,
    const std::vector<Complex> &prev2Zvals, const MeshConnectivity &mesh,
    const RenderConfig &config, const int &i,
    std::vector<std::vector<int>> &vert_to_faces,
    std::vector<std::vector<int>> &V_by_color) {
  int max_outer_iters = (i == 0) ? config.initialIters : config.subsequentIters;
  const double ENERGY_TOL = 1e-6;
  // double lambda = config.temporalPenalty; //temporal penalty coefficient

  copyComplexStdVectorToEigen(prevZvals, psi_current_complex);
  Eigen::VectorXcd dummy_grad;

  // Initial Energy Log
  double initial_energy = computeEnergyAndGradientObj(
      psi_current_complex, mesh, edgeOmegas, wrinkleAmps, dummy_grad);
  // std::cout << "Iter " << -1 << " | Energy: " << initial_energy << std::endl;

  for (int outer_iter = 0; outer_iter < max_outer_iters; ++outer_iter) {
    Eigen::VectorXcd psi_before_pass = psi_current_complex;

    for (int color = 0; color < V_by_color.size(); ++color) {
      const std::vector<int> &vertex_batch = V_by_color[color];
      if (vertex_batch.empty())
        continue;

      tbb::parallel_for(
          tbb::blocked_range<int>(0, vertex_batch.size()),
          [&](const tbb::blocked_range<int> &r) {
            for (int v_local_idx = r.begin(); v_local_idx < r.end();
                 ++v_local_idx) {
              int v_idx = vertex_batch[v_local_idx];

              // 1. Calculate Target H
              // std::cout << v_idx << std::endl;
              Complex H = computeOptimalVertexTargetObj(
                  v_idx, psi_current_complex, mesh, edgeOmegas, wrinkleAmps,
                  edgeLengths, edgeAreas, vert_to_faces);

              H *= config.meshThickness * config.beta / 4.0;

              // double dampingK = 15; 
              double dampingK = 15; // works for cloth drape under gravity

              // double dampingK = 1 / config.timeStep;
              Complex zNext = prevZvals[v_idx] *
                              std::pow(prevZvals[v_idx] / prev2Zvals[v_idx],
                                       1 - dampingK * config.timeStep);

              // inertial term

              // with velocity
              H += 0.5 * vertexMasses[v_idx] * wrinkleAmps[v_idx] *
                   wrinkleAmps[v_idx] * zNext /
                   (config.timeStep * config.timeStep);

              // with no velocity
              /*               H += 0.5 * vertexMasses[v_idx] *
                 wrinkleAmps[v_idx] * wrinkleAmps[v_idx] * (prevZvals[v_idx]) /
                                 (config.timeStep * config.timeStep); */
              if (std::abs(H) > 1e-16) {
                Complex z_target = H / std::abs(H);
                psi_current_complex(v_idx) = z_target; // no damping
                // add damping here
                // double dampingK = 10;
                /* Complex dzdt = (z_target / psi_current_complex(v_idx));

                psi_current_complex(v_idx) *=
                    std::pow(dzdt, (1 - dampingK * 0.01));
                psi_current_complex(v_idx) /=
                    std::abs(psi_current_complex(v_idx)); */
              }
            }
          });
    }

    double change = (psi_current_complex - psi_before_pass).norm();

    // std::cout << "Iter " << outer_iter << " | Energy: " << current_energy <<
    // " | Change: " << change << std::endl;

    if (change < ENERGY_TOL && outer_iter > 0)
      break;
  }
  copyComplexEigenToStdVector(psi_current_complex, curZ);
}