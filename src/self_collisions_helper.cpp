#include <argus/self_collisions_helper.hpp>
#include <argus/self_collisions.hpp>
namespace argus
{
	// self collision related stuff
	void initializeStaticWeights(PrecomputedSelfCollisionData& staticData)
	{
		// We use a local array to define the points, then copy to the aligned buffer.
		// Each row is {w0, w1, w2}
		const Real weights[15][3] = {
			// --- 3 Vertices ---
			{ 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 },

			// --- 3 Mid-edges ---
			{ 0.5, 0.5, 0.0 }, { 0.0, 0.5, 0.5 }, { 0.5, 0.0, 0.5 },

			// --- 1 Barycenter ---
			{ 0.3333333333, 0.3333333333, 0.3333333333 },

			// --- 8 Inner Points (Symmetrically distributed) ---
			{ 0.6, 0.2, 0.2 }, { 0.2, 0.6, 0.2 }, { 0.2, 0.2, 0.6 }, // Inner Ring A
			{ 0.4, 0.4, 0.2 }, { 0.2, 0.4, 0.4 }, { 0.4, 0.2, 0.4 }, // Inner Ring B
			{ 0.1, 0.45, 0.45 }, { 0.45, 0.1, 0.45 } // Fillers
		};

		// Copy into the aligned flattenedWeights buffer
		for (int i = 0; i < 15; ++i)
		{
			staticData.flattenedWeights[i * 3 + 0] = weights[i][0];
			staticData.flattenedWeights[i * 3 + 1] = weights[i][1];
			staticData.flattenedWeights[i * 3 + 2] = weights[i][2];
		}
	}
	void precomputeRestRadiiAndGrid(
	    const std::vector<trimesh::triangle_t>& faces,
	    const VectorN& Xrest,
	    PrecomputedSelfCollisionData& staticData,
	    Real safetyMultiplier, // 0.5 makes them touch, 0.55 ensures overlap
	    Real minThickness) // Absolute floor (e.g., 1mm)
	{
		const int nFaces = faces.size();
		// const int K = staticData.K; // 15

		// 1. Prepare buffers
		staticData.flattenedRadii.assign(nFaces * PrecomputedSelfCollisionData::K, 0.0);
		Real globalMaxRadius = 0.0;

// 2. Parallel loop over faces to calculate rest-distances
#pragma omp parallel
		{
			Real localMaxRadius = 0.0;

#pragma omp for schedule(static)
			for (int f = 0; f < nFaces; ++f)
			{
				const auto& face = faces[f];
				Vector3 restPos[PrecomputedSelfCollisionData::K];
				// A. Calculate world-space positions of all 15 probes in rest state
				for (int i = 0; i < PrecomputedSelfCollisionData::K; ++i)
				{
					Real w0 = staticData.flattenedWeights[i * 3 + 0];
					Real w1 = staticData.flattenedWeights[i * 3 + 1];
					Real w2 = staticData.flattenedWeights[i * 3 + 2];

					restPos[i] = w0 * Xrest.segment<3>(3 * face.v[0]) + w1 * Xrest.segment<3>(3 * face.v[1]) + w2 * Xrest.segment<3>(3 * face.v[2]);
				}

				// B. For each probe, find distance to nearest neighbor on SAME face
				for (int i = 0; i < PrecomputedSelfCollisionData::K; ++i)
				{
					Real minSqDist = std::numeric_limits<Real>::max();

					for (int j = 0; j < PrecomputedSelfCollisionData::K; ++j)
					{
						if (i == j)
							continue;
						Real d2 = (restPos[i] - restPos[j]).squaredNorm();
						if (d2 < minSqDist)
							minSqDist = d2;
					}

					// C. Calculate radius and clamp
					Real radius = std::sqrt(minSqDist) * safetyMultiplier;
					if (radius < minThickness)
						radius = minThickness;

					// D. Store and track local max
					staticData.flattenedRadii[f * PrecomputedSelfCollisionData::K + i] = radius;
					if (radius > localMaxRadius)
						localMaxRadius = radius;
				}
			}

// E. Thread-safe update of the global max radius
#pragma omp critical
			{
				if (localMaxRadius > globalMaxRadius)
					globalMaxRadius = localMaxRadius;
			}
		}

		// 3. Finalize Grid Cell Size (h)
		staticData.maxRadius = globalMaxRadius;

		// h must be at least 2 * maxRadius to ensure points in adjacent cells
		// can 'see' each other's repulsion volume.
		// 2.1 provides a small buffer for floating point safety.
		staticData.gridCellSize = globalMaxRadius * 2.1;
	}

	void initSelfCollision(const std::vector<trimesh::triangle_t>& faces, const VectorN& Xrest, PrecomputedSelfCollisionData& selfCollisionData)
	{
		// 1. Fill the 15 weights into staticData.flattenedWeights using your helper
		initializeStaticWeights(selfCollisionData);

		// 2. Compute radii per face/probe and determine the grid size (h) using your helper
		// This fills selfCollisionData.flattenedRadii and selfCollisionData.gridCellSize
		precomputeRestRadiiAndGrid(
		    faces,
		    Xrest,
		    selfCollisionData,
		    0.55, // safetyMultiplier
		    0.0001 // minThickness
		);

		// 3. One-time Setup: Initialize the Point Cloud Buffer Metadata
		// This "bakes" the static info into the 128-byte aligned structs
		// so the solver doesn't have to look them up later.
		const int nFaces = faces.size();
		const int K = selfCollisionData.K;
		selfCollisionData.pointCloudBuffer.resize(nFaces * K);

#pragma omp parallel for schedule(static)
		for (int f = 0; f < nFaces; ++f)
		{
			const auto& face = faces[f];
			for (int i = 0; i < K; ++i)
			{
				auto& pt = selfCollisionData.pointCloudBuffer[f * K + i];

				// --- Store Static Data (Never changes during simulation) ---
				pt.radius = selfCollisionData.flattenedRadii[f * K + i];
				pt.faceId = f;
				pt.vIdx[0] = face.v[0];
				pt.vIdx[1] = face.v[1];
				pt.vIdx[2] = face.v[2];

				// Store weights locally in the point for fast repulsion math
				pt.b[0] = selfCollisionData.flattenedWeights[i * 3 + 0];
				pt.b[1] = selfCollisionData.flattenedWeights[i * 3 + 1];
				pt.b[2] = selfCollisionData.flattenedWeights[i * 3 + 2];

				// Initialize position to zero (will be updated in the dynamic Gather phase)
				pt.pos.setZero();

				// Clear padding for memory alignment safety
				for (int p = 0; p < 7; ++p)
					pt.padding[p] = 0.0;
			}
		}
	}
	void resolveSelfCollisionsSubstepped(
	    const std::vector<trimesh::triangle_t>& faces,
	    const VectorN& xt, // Position at start of frame
	    VectorN& y, // Target position
	    VectorN& v, // Frame velocity
	    int maxSubsteps,
	    const VectorN& invInertia,
	    PrecomputedSelfCollisionData& sd,
	    Real stiffness,
	    Real damping,
	    Real totalDt)
	{
		Real substepDt = totalDt / static_cast<Real>(maxSubsteps);

		// 1. SYNC INITIAL VELOCITY
		// Ensure the velocity we start with matches the total displacement predicted
		v = (y - xt) / totalDt;

		VectorN currentPos = xt;

		for (int i = 0; i < maxSubsteps; ++i)
		{
			// 2. PREDICT SUBSTEP POSITION
			VectorN nextPos = currentPos + v * substepDt;

			// 3. RESOLVE (Solver now uses and updates 'v' correctly for 'substepDt')
			selfcollisionhandler::solveAllSelfCollisionsRepulsions(
			    faces,
			    currentPos, // Xold
			    nextPos, // Xnew
			    v, // V (Correctly scaled to substepDt inside)
			    invInertia,
			    sd,
			    stiffness,
			    damping,
			    substepDt // Passing the tiny DT
			);

			// 4. ADVANCE
			currentPos = nextPos;
		}

		// 5. FINAL SYNC
		y = currentPos;
	}
}