#include "../../include/ComplexLoop/ComplexLoop.h"
#include "../../include/CommonTools.h"
#include <iostream>
#include <cassert>
#include <memory>
#include <chrono>
#include <cstdlib>
#include <string>

bool
ComplexLoop::IsVertRegular(int vert) const
{
	if (_mesh.IsVertBoundary(vert)) return true;
	return (_mesh.GetVertEdges(vert).size() == 6);
}

bool
ComplexLoop::AreIrregularVertsIsolated() const
{
	for (int edge = 0; edge < _mesh.GetEdgeCount(); ++edge)
	{
		const std::vector<int>& eVerts = _mesh.GetEdgeVerts(edge);
		if (IsVertRegular(eVerts[0])) continue;
		if (IsVertRegular(eVerts[1])) continue;
		return false;
	}
	return true;
}

int
ComplexLoop::_GetVertVertIndex(int vert) const
{
	return vert;
}

int
ComplexLoop::_GetEdgeVertIndex(int edge) const
{

	return _mesh.GetVertCount() + edge;
}

int
ComplexLoop::_GetFaceVertIndex(int /*face*/) const
{
	assert(false);
	return -1;
}

int
ComplexLoop::_GetEdgeEdgeIndex(int edge, int vertInEdge) const
{
	return 2 * edge + vertInEdge;
}

int
ComplexLoop::_GetFaceEdgeIndex(int face, int edgeInFace) const
{

	assert(_mesh.IsTriangulated());
	return 2 * _mesh.GetEdgeCount() + 3 * face + edgeInFace;
}

int
ComplexLoop::_GetCentralFaceIndex(int face) const
{

	assert(_mesh.IsTriangulated());
	return 4 * face + 3;
}

int
ComplexLoop::_GetCornerFaceIndex(int face, int vertInFace) const
{

	assert(_mesh.IsTriangulated());
	return 4 * face + vertInFace;
}

// void
// ComplexLoop::GetSubdividedEdges(std::vector< std::vector<int> >& edgeToVert) const
// {

// 	assert(_mesh.IsTriangulated());

// 	int E = _mesh.GetEdgeCount();
// 	int F = _mesh.GetFaceCount();
// 	edgeToVert.resize(2 * E + 3 * F);

// 	for (int edge = 0; edge < E; ++edge)
// 	{
// 		const std::vector<int>& eVerts = _mesh.GetEdgeVerts(edge);
// 		for (int i = 0; i < 2; ++i)
// 		{
// 			int v0 = _GetVertVertIndex(eVerts[i]);
// 			int v1 = _GetEdgeVertIndex(edge);
// 			if (v0 > v1) std::swap(v0, v1);

// 			int index = _GetEdgeEdgeIndex(edge, i);
// 			edgeToVert[index].push_back(v0);
// 			edgeToVert[index].push_back(v1);
// 		}
// 	}

// 	for (int face = 0; face < F; ++face)
// 	{
// 		const std::vector<int>& fEdges = _mesh.GetFaceEdges(face);
// 		for (int i = 0; i < 3; ++i)
// 		{
// 			int v0 = _GetEdgeVertIndex(fEdges[i]);
// 			int v1 = _GetEdgeVertIndex(fEdges[(i + 1) % 3]);
// 			if (v0 > v1) std::swap(v0, v1);

// 			int index = _GetFaceEdgeIndex(face, i);
// 			edgeToVert[index].push_back(v0);
// 			edgeToVert[index].push_back(v1);
// 		}
// 	}
// }

void ComplexLoop::GetSubdividedEdges(std::vector<std::vector<int>>& edgeToVert) const
{
    assert(_mesh.IsTriangulated());

    int E = _mesh.GetEdgeCount();
    int F = _mesh.GetFaceCount();
    edgeToVert.clear();
    edgeToVert.resize(2 * E + 3 * F);

    // pre-reserve expected small size for each spot
    for (auto &v : edgeToVert) v.reserve(2);

    // Parallelize edge loop (each index is unique)
    #pragma omp parallel for schedule(static)
    for (int edge = 0; edge < E; ++edge)
    {
        const std::vector<int>& eVerts = _mesh.GetEdgeVerts(edge);
        for (int i = 0; i < 2; ++i)
        {
            int v0 = _GetVertVertIndex(eVerts[i]);
            int v1 = _GetEdgeVertIndex(edge);
            if (v0 > v1) std::swap(v0, v1);

            int index = _GetEdgeEdgeIndex(edge, i);
            // only this thread writes edgeToVert[index], so push_back is safe
            edgeToVert[index].push_back(v0);
            edgeToVert[index].push_back(v1);
        }
    }

    // Parallelize face loop too
    #pragma omp parallel for schedule(static)
    for (int face = 0; face < F; ++face)
    {
        const std::vector<int>& fEdges = _mesh.GetFaceEdges(face);
        for (int i = 0; i < 3; ++i)
        {
            int v0 = _GetEdgeVertIndex(fEdges[i]);
            int v1 = _GetEdgeVertIndex(fEdges[(i + 1) % 3]);
            if (v0 > v1) std::swap(v0, v1);

            int index = _GetFaceEdgeIndex(face, i);
            edgeToVert[index].push_back(v0);
            edgeToVert[index].push_back(v1);
        }
    }
}


// void
// ComplexLoop::GetSubdividedFaces(std::vector< std::vector<int> >& faceToVert)
// {

// 	assert(_mesh.IsTriangulated());

// 	int V = _mesh.GetVertCount();
// 	int F = _mesh.GetFaceCount();
// 	faceToVert.resize(4 * F);

// 	std::vector<int> faceFlagsNew;
// 	faceFlagsNew.resize(4 * F, 0);

// 	for (int face = 0; face < _mesh.GetFaceCount(); ++face)
// 	{
// 		int central = _GetCentralFaceIndex(face);
// 		const std::vector<int>& fVerts = _mesh.GetFaceVerts(face);
// 		const std::vector<int>& fEdges = _mesh.GetFaceEdges(face);
// 		for (int j = 0; j < 3; ++j)
// 		{
// 			// Corner face
// 			int index = _GetCornerFaceIndex(face, j);
// 			faceToVert[index].push_back(_GetVertVertIndex(fVerts[j]));
// 			faceToVert[index].push_back(_GetEdgeVertIndex(fEdges[j]));
// 			faceToVert[index].push_back(_GetEdgeVertIndex(fEdges[(j + 2) % 3]));
// 			// Central face
// 			faceToVert[central].push_back(_GetEdgeVertIndex(fEdges[j]));
// 		}
// 	}
// }

void ComplexLoop::GetSubdividedFaces(std::vector<std::vector<int>>& faceToVert)
{
    assert(_mesh.IsTriangulated());

    int V = _mesh.GetVertCount();
    int F = _mesh.GetFaceCount();
    faceToVert.clear();
    faceToVert.resize(4 * F);

    // reserve: corners need 3, central needs 3
    for (int face = 0; face < F; ++face)
    {
        int base = 4 * face;
        // corners at base + 0..2 -> each gets 3 ints
        faceToVert[base + 0].reserve(3);
        faceToVert[base + 1].reserve(3);
        faceToVert[base + 2].reserve(3);
        // central at base + 3
        faceToVert[base + 3].reserve(3);
    }

    #pragma omp parallel for schedule(static)
    for (int face = 0; face < F; ++face)
    {
        int central = _GetCentralFaceIndex(face);
        const std::vector<int>& fVerts = _mesh.GetFaceVerts(face);
        const std::vector<int>& fEdges = _mesh.GetFaceEdges(face);
        for (int j = 0; j < 3; ++j)
        {
            int cornerIndex = _GetCornerFaceIndex(face, j);
            // Corner face
            faceToVert[cornerIndex].push_back(_GetVertVertIndex(fVerts[j]));
            faceToVert[cornerIndex].push_back(_GetEdgeVertIndex(fEdges[j]));
            faceToVert[cornerIndex].push_back(_GetEdgeVertIndex(fEdges[(j + 2) % 3]));
            // Central face - each face only writes to its own central
            faceToVert[central].push_back(_GetEdgeVertIndex(fEdges[j]));
        }
    }
}


Scalar
ComplexLoop::_GetAlpha(int vert) const
{

	assert(!_mesh.IsVertBoundary(vert));
	const std::vector<int>& vEdges = _mesh.GetVertEdges(vert);

	// Fig5 left [Wang et al. 2006]
	Scalar alpha = 0.375;
	if (vEdges.size() == 3) alpha /= 2;
	else                    alpha /= vEdges.size();
	return alpha;
}

Scalar
ComplexLoop::_GetBeta(int vert) const
{

	assert(!_mesh.IsVertBoundary(vert));
	const std::vector<int>& vFaces = _mesh.GetVertFaces(vert);

	// Fig5 right [Wang et al. 2006] 
	Scalar beta = 0.;
	if (vFaces.size() >= 6) beta = 0.25;
	else if (vFaces.size() == 5) beta = 0.25 - 0.0625 * std::pow(std::sin(0.4 * M_PI), 2);
	else if (vFaces.size() == 4) beta = 0.125;
	else if (vFaces.size() == 3) beta = 0.25 / 3.;
	else assert(false);
	return beta;
}

void
ComplexLoop::BuildS0(SparseMatrixX& A) const
{

	assert(_mesh.IsTriangulated());

	std::vector<TripletX> triplets;
	int V = _mesh.GetVertCount();
	int E = _mesh.GetEdgeCount();

	// Even verts
	for (int vi = 0; vi < V; ++vi)
	{
		if (_mesh.IsVertBoundary(vi))
			_AssembleVertEvenBoundary(vi, std::back_inserter(triplets));
		else
			_AssembleVertEvenInterior(vi, std::back_inserter(triplets));
	}

	// Odd verts
	for (int edge = 0; edge < E; ++edge)
	{
		if (_mesh.IsEdgeBoundary(edge))
			_AssembleVertOddBoundary(edge, std::back_inserter(triplets));
		else
			_AssembleVertOddInterior(edge, std::back_inserter(triplets));
	}

	A.resize(V + E, V);
	A.setFromTriplets(triplets.begin(), triplets.end());
}

void
ComplexLoop::_AssembleVertEvenInterior(int vi, TripletInserter out) const
{
	// Fig5 left [Wang et al. 2006]   
	Scalar alpha = _GetAlpha(vi);
	int row = _GetVertVertIndex(vi);
	const std::vector<int>& edges = _mesh.GetVertEdges(vi);
	for (int k = 0; k < edges.size(); ++k)
	{
		int edge = edges[k];
		int viInEdge = _mesh.GetVertIndexInEdge(edge, vi);
		int vj = _mesh.GetEdgeVerts(edge)[(viInEdge + 1) % 2];
		*out++ = TripletX(row, vj, alpha);
	}
	*out++ = TripletX(row, vi, 1. - alpha * edges.size());
}

void
ComplexLoop::_AssembleVertEvenBoundary(int vi, TripletInserter out) const
{
	// Fig8 mid-top [Wang et al. 2006] 
	std::vector<int> boundary(2);
	boundary[0] = _mesh.GetVertEdges(vi).front();
	boundary[1] = _mesh.GetVertEdges(vi).back();

	int row = _GetVertVertIndex(vi);

	if(_isFixBnd)
		*out++ = TripletX(row, vi, 1.0);
	else
	{
		for (int j = 0; j < boundary.size(); ++j)
		{
			int edge = boundary[j];
			assert(_mesh.IsEdgeBoundary(edge));
			int viInEdge = _mesh.GetVertIndexInEdge(edge, vi);
			int vj = _mesh.GetEdgeVerts(edge)[(viInEdge + 1) % 2];
			*out++ = TripletX(row, vj, 0.125);
		}
		*out++ = TripletX(row, vi, 0.75);
	}
	


}

void
ComplexLoop::_AssembleVertOddInterior(int edge, TripletInserter out) const
{
	// Fig4 left-bot [Wang et al. 2006]
	for (int j = 0; j < 2; ++j)
	{
		int face = _mesh.GetEdgeFaces(edge)[j];
		int offset = _mesh.GetEdgeIndexInFace(face, edge);

		int vi = _mesh.GetFaceVerts(face)[(offset + 0) % 3];
		int vj = _mesh.GetFaceVerts(face)[(offset + 1) % 3];
		int vk = _mesh.GetFaceVerts(face)[(offset + 2) % 3];

		int row = _GetEdgeVertIndex(edge);
		*out++ = TripletX(row, vi, 0.1875);
		*out++ = TripletX(row, vj, 0.1875);
		*out++ = TripletX(row, vk, 0.125);
	}
}

void
ComplexLoop::_AssembleVertOddBoundary(int edge, TripletInserter out) const
{
	// Fig8 mid-bot [Wang et al. 2006]
	int row = _GetEdgeVertIndex(edge);
	for (int j = 0; j < 2; ++j)
	{
		int vj = _mesh.GetEdgeVerts(edge)[j];
		*out++ = TripletX(row, vj, 0.5);
	}
}

// void
// ComplexLoop::BuildS1(SparseMatrixX& A) const
// {

// 	assert(_mesh.IsTriangulated());

// 	std::vector<TripletX> triplets;
// 	int E = _mesh.GetEdgeCount();
// 	int F = _mesh.GetFaceCount();

// 	for (int edge = 0; edge < E; ++edge)
// 	{
// 		if (_mesh.IsEdgeBoundary(edge))
// 		{
// 			_AssembleEdgeEvenBoundary(edge, 0, std::back_inserter(triplets));
// 			_AssembleEdgeEvenBoundary(edge, 1, std::back_inserter(triplets));
// 		}
// 		else
// 		{
// 			const std::vector<int>& eVerts = _mesh.GetEdgeVerts(edge);
// 			for (int i = 0; i < eVerts.size(); ++i)
// 			{
// 				if (_mesh.IsVertBoundary(eVerts[i]))
// 					_AssembleEdgeEvenPartialBoundary(edge, i, std::back_inserter(triplets));
// 				else
// 					_AssembleEdgeEvenInterior(edge, i, std::back_inserter(triplets));
// 			}
// 		}
// 	}

// 	for (int face = 0; face < F; ++face)
// 	{
// 		_AssembleEdgeOdd(face, 0, std::back_inserter(triplets));
// 		_AssembleEdgeOdd(face, 1, std::back_inserter(triplets));
// 		_AssembleEdgeOdd(face, 2, std::back_inserter(triplets));
// 	}

// 	A.resize(2 * E + 3 * F, E);
// 	A.setFromTriplets(triplets.begin(), triplets.end());
// }

#include <omp.h>

// Parallel BuildS1 using per-thread triplet buckets
void
ComplexLoop::BuildS1(SparseMatrixX& A) const
{
    assert(_mesh.IsTriangulated());

    const int E = _mesh.GetEdgeCount();
    const int F = _mesh.GetFaceCount();

    // Determine number of threads we will use
    int maxThreads = 1;
    #ifdef _OPENMP
    #pragma omp parallel
    {
        #pragma omp master
        maxThreads = omp_get_num_threads();
    }
    #endif

    // Per-thread triplet buffers to avoid concurrent push_back
    std::vector<std::vector<TripletX>> threadTriplets;
    threadTriplets.resize(maxThreads);

    // Heuristic reserve: small guess to reduce reallocations.
    // (tweak if you can estimate average triplets per edge/face)
    for (int t = 0; t < maxThreads; ++t) threadTriplets[t].reserve(256);

    // Parallelize edge processing
    #pragma omp parallel for schedule(dynamic)
    for (int edge = 0; edge < E; ++edge)
    {
        int tid = 0;
        #ifdef _OPENMP
        tid = omp_get_thread_num();
        #endif

        auto &local = threadTriplets[tid];

        if (_mesh.IsEdgeBoundary(edge))
        {
            _AssembleEdgeEvenBoundary(edge, 0, std::back_inserter(local));
            _AssembleEdgeEvenBoundary(edge, 1, std::back_inserter(local));
        }
        else
        {
            const std::vector<int>& eVerts = _mesh.GetEdgeVerts(edge);
            const int vsz = static_cast<int>(eVerts.size());
            for (int i = 0; i < vsz; ++i)
            {
                if (_mesh.IsVertBoundary(eVerts[i]))
                    _AssembleEdgeEvenPartialBoundary(edge, i, std::back_inserter(local));
                else
                    _AssembleEdgeEvenInterior(edge, i, std::back_inserter(local));
            }
        }
    }

    // Parallelize face processing (each face writes 3 triplets via helpers)
    #pragma omp parallel for schedule(dynamic)
    for (int face = 0; face < F; ++face)
    {
        int tid = 0;
        #ifdef _OPENMP
        tid = omp_get_thread_num();
        #endif

        auto &local = threadTriplets[tid];

        _AssembleEdgeOdd(face, 0, std::back_inserter(local));
        _AssembleEdgeOdd(face, 1, std::back_inserter(local));
        _AssembleEdgeOdd(face, 2, std::back_inserter(local));
    }

    // Merge per-thread triplets into single vector
    size_t totalTriplets = 0;
    for (int t = 0; t < maxThreads; ++t) totalTriplets += threadTriplets[t].size();

    std::vector<TripletX> triplets;
    triplets.reserve(totalTriplets);
    for (int t = 0; t < maxThreads; ++t)
    {
        auto &local = threadTriplets[t];
        if (!local.empty())
        {
            triplets.insert(triplets.end(), local.begin(), local.end());
            // optional: free memory of local
            std::vector<TripletX>().swap(local);
        }
    }

    // Build sparse matrix (single-threaded)
    A.resize(2 * E + 3 * F, E);
    A.setFromTriplets(triplets.begin(), triplets.end());
}

void
ComplexLoop::_AssembleEdgeEvenBoundary(int edge, int vertInEdge, TripletInserter out) const
{
	int row = _GetEdgeEdgeIndex(edge, vertInEdge);
	int vert = _mesh.GetEdgeVerts(edge)[vertInEdge];
	int rSign = (_GetVertVertIndex(vert) < _GetEdgeVertIndex(edge)) ? 1 : -1;

	int nEdge = _mesh.GetVertEdges(vert).front();
	int vertInNedge = _mesh.GetVertIndexInEdge(nEdge, vert);
	int nSign = _mesh.GetVertSignInEdge(nEdge, vertInNedge);

	int pEdge = _mesh.GetVertEdges(vert).back();
	int vertInPedge = _mesh.GetVertIndexInEdge(pEdge, vert);
	int pSign = _mesh.GetVertSignInEdge(pEdge, vertInPedge);

	assert(edge == nEdge || edge == pEdge);

	if (edge == nEdge)
	{
		// Symmetric case of Fig8 right in [Wang et al. 2006]
		if(_isFixBnd)
			*out++ = TripletX(row, nEdge, (nSign == rSign) ? -0.5 : 0.5);
		else
		{
			*out++ = TripletX(row, nEdge, (nSign == rSign) ? -0.375 : 0.375);
			*out++ = TripletX(row, pEdge, (pSign == rSign) ? 0.125 : -0.125);
		}

	}
	else
	{
		// Fig8 right in [Wang et al. 2006]
		if (_isFixBnd)
			*out++ = TripletX(row, pEdge, (pSign == rSign) ? -0.5 : 0.5);
		else
		{
			*out++ = TripletX(row, pEdge, (pSign == rSign) ? -0.375 : 0.375);
			*out++ = TripletX(row, nEdge, (nSign == rSign) ? 0.125 : -0.125);
		}
	}
}

void
ComplexLoop::_InsertEdgeEdgeValue(int row, int col, int vert, int rSign, Scalar val, TripletInserter out) const
{
	// Handy function that sets the sign of val for an edge col incident to vert.
	int vertInCol = _mesh.GetVertIndexInEdge(col, vert);
	int sign = _mesh.GetVertSignInEdge(col, vertInCol);
	*out++ = TripletX(row, col, (sign == rSign) ? -val : val);
}

void
ComplexLoop::_InsertEdgeFaceValue(int row, int face, int vert, int rSign, Scalar val, TripletInserter out) const
{
	// Handy function that sets the sign of val for an edge col incident to face.
	int vertInFace = _mesh.GetVertIndexInFace(face, vert);
	int colInFace = (vertInFace + 1) % 3;
	int col = _mesh.GetFaceEdges(face)[colInFace];
	int sign = _mesh.GetEdgeSignInFace(face, colInFace);
	*out++ = TripletX(row, col, (sign == rSign) ? val : -val);
}

void
ComplexLoop::_AssembleEdgeEvenPartialBoundary(int edge, int vertInEdge, TripletInserter out) const
{
	int row = _GetEdgeEdgeIndex(edge, vertInEdge);
	int vert = _mesh.GetEdgeVerts(edge)[vertInEdge];
	int rSign = (_GetVertVertIndex(vert) < _GetEdgeVertIndex(edge)) ? 1 : -1;

	const std::vector<int>& vEdges = _mesh.GetVertEdges(vert);
	const int edgeCount = vEdges.size();

	const std::vector<int>& vFaces = _mesh.GetVertFaces(vert);
	const int faceCount = vFaces.size();

	int edgeInVert = _mesh.GetEdgeIndexInVert(vert, edge);
	assert(edgeInVert != edgeCount - 1);
	assert(edgeInVert != 0);

	std::vector< std::pair<int, Scalar> > eValues;
	std::vector< std::pair<int, Scalar> > fValues;

	assert(faceCount > 1);

	if (faceCount == 2)
	{
		// Case not covered in [Wang et al. 2006]
		if (_isFixBnd)
		{
			eValues.push_back(std::make_pair(vEdges[0], 0.53125 / 3));
			eValues.push_back(std::make_pair(vEdges[1], 0.8125 / 3));
			eValues.push_back(std::make_pair(vEdges[2], 0.53125 / 3));
		}
		else
		{
			eValues.push_back(std::make_pair(vEdges[0], 0.15625 / 3));
			eValues.push_back(std::make_pair(vEdges[1], 0.8125 / 3));
			eValues.push_back(std::make_pair(vEdges[2], 0.15625 / 3));
		}
		

		fValues.push_back(std::make_pair(vFaces[0], 0.15625 / 3));
		fValues.push_back(std::make_pair(vFaces[1], -0.15625 / 3));
	}
	else if (faceCount == 3 && edgeInVert == 1)
	{
		// Fig10 mid-left in [Wang et al. 2006]
		if (_isFixBnd)
		{
			eValues.push_back(std::make_pair(vEdges[0], 0.53125 / 3));
			eValues.push_back(std::make_pair(vEdges[1], 0.3125));
			eValues.push_back(std::make_pair(vEdges[2], 0.09375));
			eValues.push_back(std::make_pair(vEdges[3], 0.125 / 3));
		}
		else
		{
			eValues.push_back(std::make_pair(vEdges[0], 0.15625 / 3));
			eValues.push_back(std::make_pair(vEdges[1], 0.3125));
			eValues.push_back(std::make_pair(vEdges[2], 0.09375));
			eValues.push_back(std::make_pair(vEdges[3], -0.25 / 3));
		}

		fValues.push_back(std::make_pair(vFaces[0], 0.15625 / 3));
		fValues.push_back(std::make_pair(vFaces[1], -0.03125 / 3));
		fValues.push_back(std::make_pair(vFaces[2], -0.125 / 3));
	}
	else if (faceCount == 3 && edgeInVert == 2)
	{
		// Symmetric case of Fig10 mid-left in [Wang et al. 2006]
		
		if (_isFixBnd)
		{
			eValues.push_back(std::make_pair(vEdges[3], 0.53125 / 3));
			eValues.push_back(std::make_pair(vEdges[2], 0.3125));
			eValues.push_back(std::make_pair(vEdges[1], 0.09375));
			eValues.push_back(std::make_pair(vEdges[0], 0.125 / 3));
		}
		else
		{
			eValues.push_back(std::make_pair(vEdges[3], 0.15625 / 3));
			eValues.push_back(std::make_pair(vEdges[2], 0.3125));
			eValues.push_back(std::make_pair(vEdges[1], 0.09375));
			eValues.push_back(std::make_pair(vEdges[0], -0.25 / 3));
		}

		fValues.push_back(std::make_pair(vFaces[2], -0.15625 / 3));
		fValues.push_back(std::make_pair(vFaces[1], 0.03125 / 3));
		fValues.push_back(std::make_pair(vFaces[0], 0.125 / 3));
	}
	else if (faceCount == 4 && edgeInVert == 2)
	{
		// Fig10 mid-right in [Wang et al. 2006]
		
		if (_isFixBnd)
		{
			eValues.push_back(std::make_pair(vEdges[0], 0.03125));
			eValues.push_back(std::make_pair(vEdges[1], 0.125));
			eValues.push_back(std::make_pair(vEdges[2], 0.3125)); // typo in [Wang et al. 2006]
			eValues.push_back(std::make_pair(vEdges[3], 0.125));
			eValues.push_back(std::make_pair(vEdges[4], 0.03125));
		}
		else
		{
			eValues.push_back(std::make_pair(vEdges[0], -0.09375));
			eValues.push_back(std::make_pair(vEdges[1], 0.125));
			eValues.push_back(std::make_pair(vEdges[2], 0.3125)); // typo in [Wang et al. 2006]
			eValues.push_back(std::make_pair(vEdges[3], 0.125));
			eValues.push_back(std::make_pair(vEdges[4], -0.09375));
		}

		fValues.push_back(std::make_pair(vFaces[0], 0.03125));
		fValues.push_back(std::make_pair(vFaces[1], 0.03125));
		fValues.push_back(std::make_pair(vFaces[2], -0.03125));
		fValues.push_back(std::make_pair(vFaces[3], -0.03125));
	}
	else if (edgeInVert == 1)
	{
		// Fig.10 mid-mid in [Wang et al. 2006]

		if (_isFixBnd)
		{
			eValues.push_back(std::make_pair(vEdges[0], 0.53125 / 3));
			eValues.push_back(std::make_pair(vEdges[1], 0.3125)); // typo in [Wang et al. 2006]
			eValues.push_back(std::make_pair(vEdges[2], 0.09375));
			eValues.push_back(std::make_pair(vEdges[3], 0.125 / 3));
		}
		else
		{
			eValues.push_back(std::make_pair(vEdges[0], 0.15625 / 3));
			eValues.push_back(std::make_pair(vEdges[1], 0.3125)); // typo in [Wang et al. 2006]
			eValues.push_back(std::make_pair(vEdges[2], 0.09375));
			eValues.push_back(std::make_pair(vEdges[3], 0.125 / 3));
			eValues.push_back(std::make_pair(vEdges.back(), -0.125));
		}

		fValues.push_back(std::make_pair(vFaces[0], 0.15625 / 3));
		fValues.push_back(std::make_pair(vFaces[1], -0.03125 / 3));
		fValues.push_back(std::make_pair(vFaces[2], -0.125 / 3));
	}
	else if (edgeInVert == edgeCount - 2)
	{
		// Symmetric case of Fig.10 mid-mid in [Wang et al. 2006]
		
		if (_isFixBnd)
		{
			eValues.push_back(std::make_pair(vEdges[edgeCount - 1], 0.53125 / 3));
			eValues.push_back(std::make_pair(vEdges[edgeCount - 2], 0.3125)); // typo in [Wang et al. 2006]
			eValues.push_back(std::make_pair(vEdges[edgeCount - 3], 0.09375));
			eValues.push_back(std::make_pair(vEdges[edgeCount - 4], 0.125 / 3));
		}
		else
		{
			eValues.push_back(std::make_pair(vEdges[edgeCount - 1], 0.15625 / 3));
			eValues.push_back(std::make_pair(vEdges[edgeCount - 2], 0.3125)); // typo in [Wang et al. 2006]
			eValues.push_back(std::make_pair(vEdges[edgeCount - 3], 0.09375));
			eValues.push_back(std::make_pair(vEdges[edgeCount - 4], 0.125 / 3));
			eValues.push_back(std::make_pair(vEdges.front(), -0.125));
		}

		fValues.push_back(std::make_pair(vFaces[faceCount - 1], -0.15625 / 3));
		fValues.push_back(std::make_pair(vFaces[faceCount - 2], 0.03125 / 3));
		fValues.push_back(std::make_pair(vFaces[faceCount - 3], 0.125 / 3));
	}
	else if (edgeInVert == 2)
	{
		// Fig10 bot-left in [Wang et al. 2006]

		if (_isFixBnd)
		{
			eValues.push_back(std::make_pair(vEdges[0], 0.03125));
			eValues.push_back(std::make_pair(vEdges[1], 0.125));
			eValues.push_back(std::make_pair(vEdges[2], 0.3125)); // typo in [Wang et al. 2006]
			eValues.push_back(std::make_pair(vEdges[3], 0.125));
			eValues.push_back(std::make_pair(vEdges[4], 0.03125));
		}
		else
		{
			eValues.push_back(std::make_pair(vEdges[0], -0.09375));
			eValues.push_back(std::make_pair(vEdges[1], 0.125));
			eValues.push_back(std::make_pair(vEdges[2], 0.3125)); // typo in [Wang et al. 2006]
			eValues.push_back(std::make_pair(vEdges[3], 0.125));
			eValues.push_back(std::make_pair(vEdges[4], 0.03125));
			eValues.push_back(std::make_pair(vEdges[edgeCount - 1], -0.125));
		}
		

		fValues.push_back(std::make_pair(vFaces[0], 0.03125));
		fValues.push_back(std::make_pair(vFaces[1], 0.03125));
		fValues.push_back(std::make_pair(vFaces[2], -0.03125));
		fValues.push_back(std::make_pair(vFaces[3], -0.03125));
	}
	else if (edgeInVert == edgeCount - 3)
	{
		// Symmetric case of Fig10 bot-left in [Wang et al. 2006]
		if (_isFixBnd)
		{
			eValues.push_back(std::make_pair(vEdges[edgeCount - 1], 0.03125));
			eValues.push_back(std::make_pair(vEdges[edgeCount - 2], 0.125));
			eValues.push_back(std::make_pair(vEdges[edgeCount - 3], 0.3125)); // typo in [Wang et al. 2006]
			eValues.push_back(std::make_pair(vEdges[edgeCount - 4], 0.125));
			eValues.push_back(std::make_pair(vEdges[edgeCount - 5], 0.03125));
		}
		else
		{
			eValues.push_back(std::make_pair(vEdges[edgeCount - 1], -0.09375));
			eValues.push_back(std::make_pair(vEdges[edgeCount - 2], 0.125));
			eValues.push_back(std::make_pair(vEdges[edgeCount - 3], 0.3125)); // typo in [Wang et al. 2006]
			eValues.push_back(std::make_pair(vEdges[edgeCount - 4], 0.125));
			eValues.push_back(std::make_pair(vEdges[edgeCount - 5], 0.03125));
			eValues.push_back(std::make_pair(vEdges[0], -0.125));
		}

		

		fValues.push_back(std::make_pair(vFaces[faceCount - 1], -0.03125));
		fValues.push_back(std::make_pair(vFaces[faceCount - 2], -0.03125));
		fValues.push_back(std::make_pair(vFaces[faceCount - 3], 0.03125));
		fValues.push_back(std::make_pair(vFaces[faceCount - 4], 0.03125));
	}
	else
	{
		// Fig10 bot-mid in [Wang et al. 2006]
		if (_isFixBnd)
		{
			eValues.push_back(std::make_pair(vEdges[edgeInVert - 2], 0.03125));
			eValues.push_back(std::make_pair(vEdges[edgeInVert - 1], 0.125));
			eValues.push_back(std::make_pair(vEdges[edgeInVert], 0.3125)); // typo in [Wang et al. 2006]
			eValues.push_back(std::make_pair(vEdges[edgeInVert + 1], 0.125));
			eValues.push_back(std::make_pair(vEdges[edgeInVert + 2], 0.03125));
		}
		else
		{
			eValues.push_back(std::make_pair(vEdges[edgeInVert - 2], 0.03125));
			eValues.push_back(std::make_pair(vEdges[edgeInVert - 1], 0.125));
			eValues.push_back(std::make_pair(vEdges[edgeInVert], 0.3125)); // typo in [Wang et al. 2006]
			eValues.push_back(std::make_pair(vEdges[edgeInVert + 1], 0.125));
			eValues.push_back(std::make_pair(vEdges[edgeInVert + 2], 0.03125));
			eValues.push_back(std::make_pair(vEdges.front(), -0.125));
			eValues.push_back(std::make_pair(vEdges.back(), -0.125));
		}
		

		fValues.push_back(std::make_pair(vFaces[edgeInVert - 2], 0.03125));
		fValues.push_back(std::make_pair(vFaces[edgeInVert - 1], 0.03125));
		fValues.push_back(std::make_pair(vFaces[edgeInVert], -0.03125));
		fValues.push_back(std::make_pair(vFaces[edgeInVert + 1], -0.03125));
	}

	for (size_t i = 0; i < eValues.size(); ++i)
	{
		_InsertEdgeEdgeValue(row, eValues[i].first, vert, rSign, eValues[i].second, out);
	}

	for (size_t i = 0; i < fValues.size(); ++i)
	{
		_InsertEdgeFaceValue(row, fValues[i].first, vert, rSign, fValues[i].second, out);
	}
}

void
ComplexLoop::_AssembleEdgeEvenInterior(int edge, int vertInEdge, TripletInserter out) const
{
	int vert = _mesh.GetEdgeVerts(edge)[vertInEdge];
	int edgeInVert = _mesh.GetEdgeIndexInVert(vert, edge);

	int row = _GetEdgeEdgeIndex(edge, vertInEdge);
	int rSign = (_GetVertVertIndex(vert) < _GetEdgeVertIndex(edge)) ? 1 : -1;

	const std::vector<int>& vEdges = _mesh.GetVertEdges(vert);
	const std::vector<int>& vFaces = _mesh.GetVertFaces(vert);

	std::vector< std::pair<int, Scalar> > eValues;
	std::vector< std::pair<int, Scalar> > fValues;

	const int count = vEdges.size();
	assert(count == vFaces.size());

	Scalar alpha = _GetAlpha(vert);
	Scalar beta = _GetBeta(vert);

	if (count == 3)
	{
		// Fig6 bot [Wang et al. 2006]
		eValues.push_back(std::make_pair(vEdges[(edgeInVert + 0) % count], 0.375 - alpha - 0.25 * beta));
		eValues.push_back(std::make_pair(vEdges[(edgeInVert + 1) % count], 0.125 * (1. + beta) - alpha));
		eValues.push_back(std::make_pair(vEdges[(edgeInVert + 2) % count], 0.125 * (1. + beta) - alpha));
	}
	else if (count == 4)
	{
		// Case not covered in [Wang et al. 2006]
		eValues.push_back(std::make_pair(vEdges[(edgeInVert + 0) % count], 0.375 - alpha - 0.25 * beta));
		eValues.push_back(std::make_pair(vEdges[(edgeInVert + 1) % count], 0.125 - alpha));
		eValues.push_back(std::make_pair(vEdges[(edgeInVert + 2) % count], 0.25 * beta - alpha));
		eValues.push_back(std::make_pair(vEdges[(edgeInVert + 3) % count], 0.125 - alpha));
	}
	else
	{
		// Fig6 top [Wang et al. 2006]
		for (int i = 0; i < count; ++i)
		{
			Scalar val = 0.;
			if (i == 0) val = 0.375 - alpha - 0.25 * beta;
			else if (i == 1 || i == count - 1) val = 0.125 - alpha;
			else if (i == 2 || i == count - 2) val = 0.125 * beta - alpha;
			else val = -alpha;
			eValues.push_back(std::make_pair(vEdges[(edgeInVert + i) % count], val));
		}
	}

	if (count == 3)
	{
		// Fig6 bot [Wang et al. 2006]
		fValues.push_back(std::make_pair(vFaces[(edgeInVert + 0) % count], -0.125 * beta));
		fValues.push_back(std::make_pair(vFaces[(edgeInVert + 2) % count], 0.125 * beta));
	}
	else
	{
		// Fig6 top [Wang et al. 2006]
		fValues.push_back(std::make_pair(vFaces[(edgeInVert + 0) % count], -0.125 * beta));
		fValues.push_back(std::make_pair(vFaces[(edgeInVert + 1) % count], -0.125 * beta));
		fValues.push_back(std::make_pair(vFaces[(edgeInVert + count - 1) % count], 0.125 * beta));
		fValues.push_back(std::make_pair(vFaces[(edgeInVert + count - 2) % count], 0.125 * beta));
	}

	for (size_t i = 0; i < eValues.size(); ++i)
	{
		_InsertEdgeEdgeValue(row, eValues[i].first, vert, rSign, eValues[i].second, out);
	}

	for (size_t i = 0; i < fValues.size(); ++i)
	{
		_InsertEdgeFaceValue(row, fValues[i].first, vert, rSign, fValues[i].second, out);
	}
}

void
ComplexLoop::_AssembleEdgeOdd(int face, int edgeInFace, TripletInserter out) const
{
	int row = _GetFaceEdgeIndex(face, edgeInFace);

	int vertInFace = (edgeInFace + 1) % 3;
	int vert = _mesh.GetFaceVerts(face)[vertInFace];

	int nEdge = _mesh.GetFaceEdges(face)[vertInFace];
	int nSign = _mesh.GetEdgeSignInFace(face, vertInFace);

	int oEdge = _mesh.GetFaceEdges(face)[(vertInFace + 1) % 3];
	int oSign = _mesh.GetEdgeSignInFace(face, (vertInFace + 1) % 3);

	int pEdge = _mesh.GetFaceEdges(face)[(vertInFace + 2) % 3];
	int pSign = _mesh.GetEdgeSignInFace(face, (vertInFace + 2) % 3);

	int rSign = (_GetEdgeVertIndex(nEdge) < _GetEdgeVertIndex(pEdge)) ? 1 : -1;

	bool nBdry = _mesh.IsEdgeBoundary(nEdge);
	bool pBdry = _mesh.IsEdgeBoundary(pEdge);

	if (nBdry && pBdry)
	{
		// Ear case not covered in [Wang et al. 2006]
		*out++ = TripletX(row, oEdge, (oSign == rSign) ? 0.25 : -0.25);
		*out++ = TripletX(row, pEdge, (pSign == rSign) ? -0.25 : 0.25);
		*out++ = TripletX(row, nEdge, (nSign == rSign) ? -0.25 : 0.25);
	}
	else if (nBdry)
	{
		// Symmetric case of Fig10 top-left in [Wang et al. 2006]
		*out++ = TripletX(row, oEdge, (oSign == rSign) ? 0.21875 : -0.21875);
		*out++ = TripletX(row, pEdge, (pSign == rSign) ? -0.1875 : 0.1875);
		*out++ = TripletX(row, nEdge, (nSign == rSign) ? -0.15625 : 0.15625);
	}
	else if (pBdry)
	{
		// Fig10 top-left in [Wang et al. 2006]
		*out++ = TripletX(row, oEdge, (oSign == rSign) ? 0.21875 : -0.21875);
		*out++ = TripletX(row, pEdge, (pSign == rSign) ? -0.15625 : 0.15625);
		*out++ = TripletX(row, nEdge, (nSign == rSign) ? -0.1875 : 0.1875);
	}
	else
	{
		// Fig4 mid-bot [Wang et al. 2006]
		*out++ = TripletX(row, oEdge, (oSign == rSign) ? 0.1875 : -0.1875);
		*out++ = TripletX(row, pEdge, (pSign == rSign) ? -0.09375 : 0.09375);
		*out++ = TripletX(row, nEdge, (nSign == rSign) ? -0.09375 : 0.09375);
	}

	std::vector<int> flaps;
	flaps.push_back(pEdge);
	flaps.push_back(nEdge);
	for (int i = 0; i < flaps.size(); ++i)
	{
		int edge = flaps[i];
		if (_mesh.IsEdgeBoundary(edge)) continue;

		int faceInEdge = _mesh.GetFaceIndexInEdge(edge, face);
		int flap = _mesh.GetEdgeFaces(edge)[(faceInEdge + 1) % 2];
		int edgeInFlap = _mesh.GetEdgeIndexInFace(flap, edge);

		for (int j = 1; j <= 2; ++j)
		{
			// Fig4 mid-bot [Wang et al. 2006]
			Scalar val = 0.09375;
			int cEdge = _mesh.GetFaceEdges(flap)[(edgeInFlap + j) % 3];
			int cSign = _mesh.GetEdgeSignInFace(flap, (edgeInFlap + j) % 3);
			if (_mesh.GetVertIndexInEdge(cEdge, vert) == -1) val = -0.03125;
			*out++ = TripletX(row, cEdge, (cSign == rSign) ? -val : val);
		}
	}
}

void ComplexLoop::meshSubdivide(int level, std::vector<SparseMatrixX> &S0List, std::vector<SparseMatrixX> &S1List, std::vector<Mesh> &meshList)
{
	static const bool timingEnabled = []() {
		const char* v = std::getenv("CWF_TIMING");
		return (v != nullptr) && (std::string(v) != "0");
	}();

	using Clock = std::chrono::high_resolution_clock;
	auto totalStart = Clock::now();

	int nverts = _mesh.GetVertCount();
	int startVerts = _mesh.GetVertCount();
	int startEdges = _mesh.GetEdgeCount();
	int startFaces = _mesh.GetFaceCount();

	MatrixX X;
	_mesh.GetPos(X);
	
	for (int l = 0; l < level; ++l)
	{
		auto iterStart = Clock::now();
		auto t0 = iterStart;

		//SparseMatrixX tmpS0 = S0List[l];
		auto t1 = Clock::now();
		//SparseMatrixX tmpS1 = S1List[l];
		auto t2 = Clock::now();

		X = S0List[l] * X;
		auto t3 = Clock::now();

		std::vector<Vector3> points;
		ConvertToVector3(X, points);
		auto t4 = Clock::now();

		// std::vector< std::vector<int> > edgeToVert;
		// GetSubdividedEdges(edgeToVert);
		auto t5 = Clock::now();

		// std::vector< std::vector<int> > faceToVert;
		// GetSubdividedFaces(faceToVert);
		auto t6 = Clock::now();

		_mesh.Populate(points);
		auto t7 = Clock::now();

		if (timingEnabled)
		{
			auto ms = [](const Clock::time_point& a, const Clock::time_point& b) {
				return std::chrono::duration_cast<std::chrono::milliseconds>(b - a).count();
			};
			std::cout << "[CWF_TIMING][ComplexLoop::meshSubdivide] level " << l << "/" << level
			          << " | BuildS0=" << ms(t0, t1) << "ms"
			          << " BuildS1=" << ms(t1, t2) << "ms"
			          << " S0*X=" << ms(t2, t3) << "ms"
			          << " ConvertToVector3=" << ms(t3, t4) << "ms"
			          << " GetSubdividedEdges=" << ms(t4, t5) << "ms"
			          << " GetSubdividedFaces=" << ms(t5, t6) << "ms"
			          << " Populate=" << ms(t6, t7) << "ms"
			          << " iterTotal=" << ms(iterStart, t7) << "ms"
			          << std::endl;
		}
	}

	if (timingEnabled)
	{
		auto totalEnd = Clock::now();
		auto totalMs = std::chrono::duration_cast<std::chrono::milliseconds>(totalEnd - totalStart).count();
		int endVerts = _mesh.GetVertCount();
		int endEdges = _mesh.GetEdgeCount();
		int endFaces = _mesh.GetFaceCount();
		std::cout << "[CWF_TIMING][ComplexLoop::meshSubdivide] total=" << totalMs << "ms"
		          << " levels=" << level
		          << " V " << startVerts << "->" << endVerts
		          << " E " << startEdges << "->" << endEdges
		          << " F " << startFaces << "->" << endFaces
		          << std::endl;
	}
}