#ifndef WRINKLE_HELPERS_HPP
#define WRINKLE_HELPERS_HPP

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
	Real computeEdgeOmegaDirect(
	    const Eigen::Vector2d& tvec1,
	    const Eigen::Vector2d& tvec2,
	    Real amp1,
	    Real amp2);

}
#endif