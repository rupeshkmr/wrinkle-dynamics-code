#pragma once
#include "ComplexLoop.h"

class ComplexLoopZuenko : public ComplexLoop	// We modify the Loop.h
{
public:
    //void virtual Subdivide(const Eigen::VectorXd& omega, const std::vector<std::complex<double>>& zvals, Eigen::VectorXd& omegaNew, std::vector<std::complex<double>>& upZvals, int level, std::vector<SparseMatrixX> &S0List, 
	// std::vector<SparseMatrixX> &S1List,
	// std::vector<Mesh> &meshList	) override;

    std::complex<double> interpZ(const std::vector<std::complex<double>>& zList, const std::vector<Eigen::Vector3d>& gradThetaList, std::vector<double>& coords, const std::vector<Eigen::Vector3d>& pList);
    
    void interpZ(
        const std::vector<std::complex<double>>& zList1, const std::vector<Eigen::Vector3d>& gradThetaList1,
        const std::vector<std::complex<double>>& zList2, const std::vector<Eigen::Vector3d>& gradThetaList2,
        const std::vector<double>& coords, const std::vector<Eigen::Vector3d>& pList,
        std::complex<double>& z1, std::complex<double>& z2);

    std::complex<double> computeZandGradZ(const Eigen::VectorXd& omega, const std::vector<std::complex<double>>& zvals, int fid, const Eigen::Vector3d& bary, Eigen::Vector3cd* gradz, Mesh &mesh);
    Eigen::Vector3d computeGradThetaFromOmegaPerface(const Eigen::VectorXd& omega, int fid, int vInF, Mesh &mesh);
    Eigen::Vector3d computeBaryGradThetaFromOmegaPerface(const Eigen::VectorXd& omega, int fid, const Eigen::Vector3d& bary, Mesh &mesh);

    void updateLoopedZvals(const Eigen::VectorXd& omega, const std::vector<std::complex<double>>& zvals, std::vector<std::complex<double>>& upZvals, Mesh &mesh);
    
    void updateLoopedZvals(
        const Eigen::VectorXd& omega1, const std::vector<std::complex<double>>& zvals1, std::vector<std::complex<double>>& upZvals1,
        const Eigen::VectorXd& omega2, const std::vector<std::complex<double>>& zvals2, std::vector<std::complex<double>>& upZvals2,
        Mesh &mesh);
};
