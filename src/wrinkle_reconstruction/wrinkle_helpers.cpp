#include <argus/wrinkle_reconstruction/wrinkle_helpers.hpp>

Real argus::computeEdgeOmegaDirect(
    const Eigen::Vector2d& tvec1,
    const Eigen::Vector2d& tvec2,
    Real amp1,
    Real amp2)
{
	// 1. Extract the projections onto the shared edge (v0-v1)
	Real omega1 = tvec1.x();
	Real omega2 = tvec2.y();

	// 2. Sign Synchronization (Local Lifting)
	// Since we are dealing with line fields, we align omega2 to omega1
	// so they are in the same 'hemisphere' before averaging.
	if (omega1 * omega2 < 0.0)
	{
		omega2 = -omega2;
	}

	// 3. Amplitude-Weighted Average
	// Prevents low-signal noise from corrupting the dominant wrinkle direction.
	Real totalAmp = amp1 + amp2;
	if (totalAmp < 1e-12)
		return 0.0;

	return (amp1 * omega1 + amp2 * omega2) / totalAmp;
}
