#pragma once
#include "commonHeader.hpp"

#include <LBFGS.h>

double computeEnergyAndGradient(
    const Eigen::VectorXcd& psi,
    const Eigen::MatrixXi& F,
    const Eigen::MatrixXd& face_grad_phis_2d,
    const Eigen::VectorXd& wa,
    Eigen::VectorXcd& grad) 
{
    int num_vertices = psi.size();
    int num_faces = F.rows();
    
    grad.setZero(num_vertices);
    double total_energy = 0.0;

    for (int f_idx = 0; f_idx < num_faces; ++f_idx) {
        int i = F(f_idx, 0);
        int j = F(f_idx, 1);
        int k = F(f_idx, 2);

        // Get the face-local 2D omega components
        double c1 = face_grad_phis_2d(f_idx, 0); // component 1
        double c2 = face_grad_phis_2d(f_idx, 1); // component 2
        
        // Edge omegas (w) derived from face components
        double w_ij = c1;
        double w_ki = -c2;
        double w_jk = c2 - c1;

        double s_ij = 0.5 * (wa[i] + wa[j]);
        double s_jk = 0.5 * (wa[j] + wa[k]);
        double s_ki = 0.5 * (wa[k] + wa[i]);


        Complex exp_ij_fwd = std::exp(Complex(0.0, w_ij));
        Complex exp_jk_fwd = std::exp(Complex(0.0, w_jk));
        Complex exp_ki_fwd = std::exp(Complex(0.0, w_ki));
        

        Complex r_ij = s_ij * (psi(j) - psi(i) * exp_ij_fwd);
        Complex r_jk = s_jk * (psi(k) - psi(j) * exp_jk_fwd);
        Complex r_ki = s_ki * (psi(i) - psi(k) * exp_ki_fwd);
        
        double E_fwd = std::norm(r_ij) + std::norm(r_jk) + std::norm(r_ki);
        

        Complex exp_ij_rev = std::conj(exp_ij_fwd);
        Complex exp_jk_rev = std::conj(exp_jk_fwd);
        Complex exp_ki_rev = std::conj(exp_ki_fwd);

        Complex nr_ij = s_ij * (psi(j) - psi(i) * exp_ij_rev);
        Complex nr_jk = s_jk * (psi(k) - psi(j) * exp_jk_rev);
        Complex nr_ki = s_ki * (psi(i) - psi(k) * exp_ki_rev);

        double E_rev = std::norm(nr_ij) + std::norm(nr_jk) + std::norm(nr_ki);
        
        total_energy += 0.5 * (E_fwd + E_rev);

        // Gradient for E_fwd:
        grad(i) += (r_ki - r_ij * exp_ij_rev); 
        grad(j) += (r_ij - r_jk * exp_jk_rev);
        grad(k) += (r_jk - r_ki * exp_ki_rev);

        // Gradient for E_rev:
        grad(i) += (nr_ki - nr_ij * exp_ij_fwd);
        grad(j) += (nr_ij - nr_jk * exp_jk_fwd);
        grad(k) += (nr_jk - nr_ki * exp_ki_fwd);

    }
    // std::cout << "Current Energy: " << total_energy << std::endl;
    return total_energy;
}


void optimizePsi(
    Eigen::VectorXcd& psi,
    const Eigen::MatrixXi& F,
    const Eigen::MatrixXd& face_grad_phis_2d,
    const Eigen::VectorXd& wrinkleAmps,
    double lr, 
    int max_iters, 
    double tol) 
{
    double prev_energy = -1.0;
    Eigen::VectorXcd grad;

    for (int it = 0; it < max_iters; ++it) {
        double energy = computeEnergyAndGradient(psi, F, face_grad_phis_2d, wrinkleAmps, grad);

        // Update psi: psi -= lr * grad
        psi -= lr * grad;

        // Projection: psi /= |psi| (Projection back to unit magnitude)
        // Eigen's cwiseAbs() computes the magnitude of each complex number
        psi = psi.cwiseQuotient(psi.cwiseAbs());

        // Convergence check (energy change)
        double change = (lr*grad).norm();
        if (change < tol && it > 0) {
            std::cout << "VBD converged after " << it + 1 << " outer iterations." << std::endl;
            break;
        }

        prev_energy = energy;
    }
}

const double TOLERANCE = 1e-6;
Eigen::VectorXcd psi_current_complex;

void optimizeAndSaveZvalsGD(
    Eigen::MatrixXd& faceOmega2D,
    std::vector<Complex>& curZ,
    Eigen::VectorXd& wrinkleAmps,
    const std::vector<Complex>& prevZvals,
    const Eigen::MatrixXi& triF,
    const RenderConfig& config,
    const int i
)
{

    int max_outer_iters;
    if (i == 0) {
        max_outer_iters = config.initialIters;
    } 
    else {
        max_outer_iters = config.subsequentIters;
    }

    copyComplexStdVectorToEigen(prevZvals, psi_current_complex);

    optimizePsi(psi_current_complex, triF, faceOmega2D, wrinkleAmps, config.learningRate, max_outer_iters, TOLERANCE);

    copyComplexEigenToStdVector(psi_current_complex, curZ);

}

double computeEnergyAndGradientObj(
    const Eigen::VectorXcd& psi,
    const MeshConnectivity& mesh, // Pass the class instance
    const Eigen::VectorXd& edgeOmegas, 
    const Eigen::VectorXd& wa,
    Eigen::VectorXcd& grad) 
{
    int num_vertices = psi.size();
    int num_faces = mesh.nFaces();
    const Eigen::MatrixXi& F = mesh.faces(); // Access faces directly from class
    
    grad.setZero(num_vertices);
    double total_energy = 0.0;

    for (int f_idx = 0; f_idx < num_faces; ++f_idx) {
        int i = F(f_idx, 0);
        int j = F(f_idx, 1);
        int k = F(f_idx, 2);

        // --- Retrieve Edge Indices (Opposite Vertex Convention) ---
        // Edge ij is opposite vertex k (index 2)
        int e_ij_idx = mesh.faceEdge(f_idx, 2);
        // Edge jk is opposite vertex i (index 0)
        int e_jk_idx = mesh.faceEdge(f_idx, 0);
        // Edge ki is opposite vertex j (index 1)
        int e_ki_idx = mesh.faceEdge(f_idx, 1);

        // --- Determine Orientation ---
        // Check if global edge definition starts with the source vertex of our traversal
        double ori_ij = (mesh.edgeVertex(e_ij_idx, 0) == i) ? 1.0 : -1.0;
        double ori_jk = (mesh.edgeVertex(e_jk_idx, 0) == j) ? 1.0 : -1.0;
        double ori_ki = (mesh.edgeVertex(e_ki_idx, 0) == k) ? 1.0 : -1.0;

        // Fetch w values and correct for orientation
        double w_ij = edgeOmegas(e_ij_idx) * ori_ij;
        double w_jk = edgeOmegas(e_jk_idx) * ori_jk;
        double w_ki = edgeOmegas(e_ki_idx) * ori_ki;

        // --- Rest of the Physics Calculation (Unchanged) ---
        double s_ij = 0.5 * (wa[i] + wa[j]);
        double s_jk = 0.5 * (wa[j] + wa[k]);
        double s_ki = 0.5 * (wa[k] + wa[i]);

        Complex exp_ij_fwd = std::exp(Complex(0.0, w_ij));
        Complex exp_jk_fwd = std::exp(Complex(0.0, w_jk));
        Complex exp_ki_fwd = std::exp(Complex(0.0, w_ki));
        
        Complex r_ij = s_ij * (psi(j) - psi(i) * exp_ij_fwd);
        Complex r_jk = s_jk * (psi(k) - psi(j) * exp_jk_fwd);
        Complex r_ki = s_ki * (psi(i) - psi(k) * exp_ki_fwd);
        
        double E_fwd = std::norm(r_ij) + std::norm(r_jk) + std::norm(r_ki);

        Complex exp_ij_rev = std::conj(exp_ij_fwd);
        Complex exp_jk_rev = std::conj(exp_jk_fwd);
        Complex exp_ki_rev = std::conj(exp_ki_fwd);
        
        total_energy += E_fwd;

        grad(i) += 2.0 * (r_ki - r_ij * exp_ij_rev);
        grad(j) += 2.0 * (r_ij - r_jk * exp_jk_rev);
        grad(k) += 2.0 * (r_jk - r_ki * exp_ki_rev);
    }

    return total_energy;
}



void optimizePsiObj(
    Eigen::VectorXcd& psi,
    const MeshConnectivity& mesh, // Pass MeshConnectivity
    const Eigen::VectorXd& edgeOmegas,
    const Eigen::VectorXd& wrinkleAmps,
    double lr, 
    int max_iters, 
    double tol) 
{
    double prev_energy = -1.0;
    Eigen::VectorXcd grad;

    for (int it = 0; it < max_iters; ++it) {
        // Use updated function signature
        double energy = computeEnergyAndGradientObj(psi, mesh, edgeOmegas, wrinkleAmps, grad);
        std::cout << energy << std::endl;

        psi -= lr * grad;
        psi = psi.cwiseQuotient(psi.cwiseAbs());

        double change = (lr*grad).norm();
        if (change < tol && it > 0) {
            std::cout << "Converged after " << it + 1 << " iterations." << std::endl;
            break;
        }
        prev_energy = energy;
    }
}


void optimizeAndSaveZvalsGDObj(
    const Eigen::VectorXd& edgeOmegas, 
    std::vector<Complex>& curZ,
    Eigen::VectorXd& wrinkleAmps,
    const std::vector<Complex>& prevZvals,
    const MeshConnectivity& mesh, // Add this parameter
    const RenderConfig& config,
    const int i
)
{
    int max_outer_iters = (i == 0) ? config.initialIters : config.subsequentIters;

    copyComplexStdVectorToEigen(prevZvals, psi_current_complex);

    // Pass mesh to optimizer
    optimizePsiObj(psi_current_complex, mesh, edgeOmegas, wrinkleAmps, config.learningRate, max_outer_iters, TOLERANCE);

    copyComplexEigenToStdVector(psi_current_complex, curZ);
}


Eigen::VectorXd x;

void optimizeAndSaveZvalsBFGS(
    Eigen::MatrixXd& faceOmega2D,
    std::vector<Complex>& curZ,
    Eigen::VectorXd& wrinkleAmps,
    const std::vector<Complex>& prevZvals,
    const Eigen::MatrixXi& triF,
    const RenderConfig& config,
    const LBFGSpp::LBFGSParam<double>& param,
    const int& i
)
{
    int nverts = prevZvals.size();
    // std::cout << "\nOptimizing Frame " << i << " using Projected L-BFGS (step-by-step)..." << std::endl;

    psi_current_complex.resize(nverts);
    x = Eigen::VectorXd::Zero(2 * nverts);

    int max_outer_iters;
    if (i == 0) {
        max_outer_iters = config.initialIters;
    } 
    else {
        max_outer_iters = config.subsequentIters;
    }

    copyComplexStdVectorToEigen(prevZvals, psi_current_complex);

    auto objective_lambda = [&](const Eigen::VectorXd& current_x, Eigen::VectorXd& grad) -> double {

        Eigen::VectorXcd psi(nverts);
        if (current_x.size() != 2 * nverts) throw std::runtime_error("LBFGS vector size mismatch");
        for (int k = 0; k < nverts; ++k) psi(k) = Complex(current_x(2 * k), current_x(2 * k + 1));
        Eigen::VectorXcd complex_grad;

        double energy = computeEnergyAndGradient(psi, triF, faceOmega2D, wrinkleAmps, complex_grad); // Use your chosen gradient func

        if (complex_grad.size() != nverts) throw std::runtime_error("Complex gradient size mismatch");
        grad.resize(2 * nverts);
        for (int k = 0; k < nverts; ++k) {
            grad(2 * k)     = complex_grad(k).real();
            grad(2 * k + 1) = complex_grad(k).imag();
        }
        return energy;
    };

    Eigen::VectorXcd psi_previous_complex;
    double change = std::numeric_limits<double>::max();
    LBFGSpp::LBFGSSolver<double> solver(param); // Create solver outside loop if state is preserved

    for (int outer_iter = 0; outer_iter < max_outer_iters; ++outer_iter)
    {
        psi_previous_complex = psi_current_complex; // Store state before LBFGS step

        for(int k=0; k < nverts; ++k) {
            x(2 * k)     = psi_current_complex(k).real();
            x(2 * k + 1) = psi_current_complex(k).imag();
        }

        int niter_inner = 0; // Will be 0 or 1
        try {
            double current_fx;
            niter_inner = solver.minimize(objective_lambda, x, current_fx);

            // std::cout << "  OuterIter " << outer_iter << ": Energy after LBFGS step = " << current_fx << std::endl;

        } catch (const std::exception& e) {
            std::cerr << "L-BFGS inner step failed at outer iter " << outer_iter << " for frame " << i << ": " << e.what() << std::endl;

            psi_current_complex = psi_previous_complex;
            break;
        }

        for (int k = 0; k < nverts; ++k) {
            psi_current_complex(k) = Complex(x(2 * k), x(2 * k + 1));
        }


        for (int k = 0; k < nverts; ++k) {
            double mag = std::abs(psi_current_complex(k));
            if (mag > 1e-12) {
                psi_current_complex(k) /= mag;
            } else {
                psi_current_complex(k) = Complex(0.0, 0.0);
            }
        }

        change = (psi_current_complex - psi_previous_complex).norm();

        if (change < TOLERANCE && outer_iter > 0) { // Check only after first step
            std::cout << "GD converged after " << outer_iter + 1 << " outer iterations." << std::endl;
            break;
        }
        if (outer_iter == max_outer_iters - 1) {
            std::cout << "Warning: GD reached max outer iterations (" << max_outer_iters << "). Change: " << change << std::endl;
        }
    }

    copyComplexEigenToStdVector(psi_current_complex, curZ);
}
