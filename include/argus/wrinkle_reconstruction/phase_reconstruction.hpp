#ifndef PHASE_RECONSTRUCTION
#define PHASE_RECONSTRUCTION

#include "../common.hpp"
#include "../sheet.hpp"
/**
 * Computes the averaged 1-form (omega) for a shared edge between two faces.
 * Optimized for Line Fields: Resolves 180-degree ambiguity.
 * * Assumptions:
 * Face 1 (0,1,2): tvec1.x() is projection on edge 0-1.
 * Face 2 (0,3,1): tvec2.y() is projection on edge 0-1.
 */
namespace argus
{
	class WrinklePhase
	{
	private:
	public:
		~WrinklePhase() = default;
		WrinklePhase(argus::Sheet& sheet);
		void initVBDPhaseReconstruction(argus::Sheet& sheet);
		Real computePhaseReconstructionEnergy(argus::Sheet& sheet, const Eigen::Vector<Real, Eigen::Dynamic>& phiold, Real timestep);
		void solvePhi(argus::Sheet& sheet, Real timestep);
		std::complex<Real> getPhaseFromFaceEnergy(argus::Sheet& sheet, int vid, Real timestep, const Eigen::Vector<std::complex<Real>, Eigen::Dynamic>& x_current);
	};
}
#endif