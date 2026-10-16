#include "Mesh.h"

#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

void
Mesh::Copy(const Mesh& mesh)
{
    _points = mesh._points;
    _faceToVert = mesh._faceToVert;
    _faceToEdge = mesh._faceToEdge;
    _edgeToVert = mesh._edgeToVert;
    _edgeToFace = mesh._edgeToFace;
    _vertToFace = mesh._vertToFace;
    _vertToEdge = mesh._vertToEdge;
    _isTriangulated = mesh._isTriangulated;
    _isQuadrangulated = mesh._isQuadrangulated;
}

void 
Mesh::Clear()
{
    _points.clear();
    _faceToVert.clear();
    _faceToEdge.clear();
    _edgeToVert.clear();
    _edgeToFace.clear();
    _vertToFace.clear();
    _vertToEdge.clear();
    _isTriangulated = false;
    _isQuadrangulated = false;
}

int 
Mesh::GetEdgeSignInFace(int face, int edgeInFace) const
{
    assert(edgeInFace < GetFaceEdges(face).size());
    int fVert = GetFaceVerts(face)[edgeInFace];
    int edge  = GetFaceEdges(face)[edgeInFace]; 
    int eVert = GetEdgeVerts(edge)[0];
    return (fVert == eVert)? 1 : -1;
}

int 
Mesh::GetVertSignInEdge(int edge, int vertInEdge) const
{
    assert(vertInEdge < GetEdgeVerts(edge).size());
    return (vertInEdge == 1)? 1 : -1;
}

bool
Mesh::HasBoundary() const
{
    for (int edge = 0; edge < GetEdgeCount(); ++edge) 
    {
        if (IsEdgeBoundary(edge)) return true;
    }
    return false;
}

void
Mesh::GetBoundaryVerts(std::vector< std::vector<int> > &verts) const
{
    std::vector<bool> isDone(GetVertCount(), false);
    for (int vert = 0; vert < GetVertCount(); ++vert) 
    {
        if (isDone[vert]) continue;
        isDone[vert] = true;

        std::vector<int> boundary;
        GetBoundaryFromVert(vert, boundary);
        
        if (boundary.empty()) continue;
        verts.push_back(boundary);

        for (int i = 0; i < boundary.size(); ++i)
        {
            isDone[boundary[i]] = true;
        }
    }
}

void 
Mesh::GetBoundaryFromVert(int vert, std::vector<int>& boundary) const
{
    if (!IsVertBoundary(vert)) return;
    int next = vert;
    do
    {
        assert(IsVertBoundary(next));
        boundary.push_back(next);

        int edge = GetVertEdges(next)[0];
        assert(IsEdgeBoundary(edge));

        int vertInEdge = GetVertIndexInEdge(edge, next);
        next = GetEdgeVerts(edge)[(vertInEdge+1)%2];
    } while (next != vert);
}


void 
Mesh::GetBoundaryVerts(std::vector<int> &verts) const
{
    for (int vert = 0; vert < GetVertCount(); ++vert) 
    {
        if (IsVertBoundary(vert)) 
        {
            verts.push_back(vert);
        }
    }
}

void
Mesh::GetNonBoundaryVerts(std::vector<int> &verts) const
{
    for (int vert = 0; vert < GetVertCount(); ++vert) 
    {
        if (!IsVertBoundary(vert)) 
        {
            verts.push_back(vert);
        }
    }
}

void 
Mesh::GetBoundaryEdges(std::vector<int> &edges) const
{
    for (int edge = 0; edge < GetEdgeCount(); ++edge) 
    {
        if (IsEdgeBoundary(edge)) 
        {
            edges.push_back(edge);
        }
    }
}

void 
Mesh::GetNonBoundaryEdges(std::vector<int> &edges) const
{
    for (int edge = 0; edge < GetEdgeCount(); ++edge) 
    {
        if (!IsEdgeBoundary(edge)) 
        {
            edges.push_back(edge);
        }
    }
}

void
Mesh::BuildD(int i, SparseMatrixX& D) const
{
    if      (i == 0) BuildD0(D);
    else if (i == 1) BuildD1(D);
    else assert(false);
}

void 
Mesh::BuildD0(SparseMatrixX &D) const
{
    int V = GetVertCount();
    int E = GetEdgeCount();

    std::vector<TripletX> coef;
    coef.reserve(2*E);

    for (int edge = 0; edge < E; ++edge) 
    {
        const std::vector<int>& eVerts = GetEdgeVerts(edge);
        for (size_t i = 0; i < eVerts.size(); ++i)
        {
            coef.push_back(TripletX(edge, eVerts[i], GetVertSignInEdge(edge, i)));
        }
    }

    D.resize(E, V);
    D.setFromTriplets(coef.begin(), coef.end());
}

void 
Mesh::BuildD1(SparseMatrixX &D) const
{
    int E = GetEdgeCount();
    int F = GetFaceCount();

    std::vector<TripletX> coef;
    coef.reserve(2*E);

    for (int face = 0; face < F; ++face) 
    {
        const std::vector<int>& fEdges = GetFaceEdges(face);
        for (size_t i = 0; i < fEdges.size(); ++i)
        {
            coef.push_back(TripletX(face, fEdges[i], GetEdgeSignInFace(face, i)));
        }
    }

    D.resize(F, E);
    D.setFromTriplets(coef.begin(), coef.end());
}

void
Mesh::Populate(
    const std::vector<Vector3> &points,
    const std::vector<int> &faceVertCount,
    const std::vector<int> &faceToVert)
{
    size_t offset = 0;
    std::vector< std::vector<int> > faceVerts(faceVertCount.size());
    for (size_t i = 0; i < faceVerts.size(); ++i)
    {
        for (size_t j = 0; j < faceVertCount[i]; ++j)
        {
            faceVerts[i].push_back(faceToVert[offset+j]);
        }
        offset += faceVertCount[i];
    }
    Populate(points, faceVerts);
}

void
Mesh::Populate(
    const std::vector<Vector3> &points,
    const std::vector< std::vector<int> > &faceToVert)
{
    std::vector< std::vector<int> > edgeToVert;
    std::map< std::pair<int,int>, int > heToEdge;
    for (int face = 0; face < faceToVert.size(); ++face) 
    {
        const std::vector<int>& fVerts = faceToVert[face];
        for (int i = 0; i < fVerts.size(); ++i) 
        {
            int vi = fVerts[i];
            int vj = fVerts[(i+1) % fVerts.size()];
            assert(vi != vj);

            std::pair<int, int> he = std::make_pair(vi, vj);
            if (he.first > he.second) std::swap(he.first, he.second);
            if (heToEdge.find(he) != heToEdge.end()) continue;

            heToEdge[he] = edgeToVert.size();
            edgeToVert.push_back(std::vector<int>(2));
            edgeToVert.back()[0] = he.first;
            edgeToVert.back()[1] = he.second;
        }
    }
    Populate(points, faceToVert, edgeToVert);
}

void
Mesh::Populate(
    const std::vector<Vector3> &points)
{
    //Clear();

    _points = std::move(points);

}

// void
// Mesh::Populate(
//     const std::vector<Vector3> &points,
//     const std::vector< std::vector<int> > &faceToVert,
//     const std::vector< std::vector<int> > &edgeToVert)
// {
//     Clear();

//     _points = points;

//     _faceToVert = faceToVert;
//     _faceToEdge.resize(_faceToVert.size());
    
//     _edgeToVert = edgeToVert;
//     _edgeToFace.resize(_edgeToVert.size());
    
//     _vertToEdge.resize(_points.size());
//     _vertToFace.resize(_points.size());

//     std::map< std::pair<int,int>, int > heToEdge;
//     for (int edge = 0; edge < _edgeToVert.size(); ++edge)
//     {
//         const std::vector<int>& eVert = _edgeToVert[edge];
//         assert(eVert.size() == 2);

//         std::pair<int, int> he = std::make_pair(eVert[0], eVert[1]);
//         assert(he.first < he.second);

//         heToEdge[he] = edge;
//         _vertToEdge[he.first ].push_back(edge); 
//         _vertToEdge[he.second].push_back(edge); 
//     }

//     for (int face = 0; face < _faceToVert.size(); ++face) 
//     {
//         const std::vector<int>& fVerts = _faceToVert[face];
//         for (int i = 0; i < fVerts.size(); ++i) 
//         {
//             int vi = fVerts[i];
//             int vj = fVerts[(i+1) % fVerts.size()];
//             assert(vi != vj);

//             std::pair<int, int> he = std::make_pair(vi, vj);
//             if (he.first > he.second) std::swap(he.first, he.second);

//             assert(heToEdge.find(he) != heToEdge.end());
//             int edge = heToEdge[he];

//             _vertToFace[vi].push_back(face);
//             _faceToEdge[face].push_back(edge);
//             _edgeToFace[edge].push_back(face);
//         }
//     }

//     _CheckValidity();

//     _SortOneRings();

//     _isTriangulated = _AreFacesTris();
//     _isQuadrangulated = _AreFacesQuads();
// }

#include <unordered_map>
#include <cstdint>
#include <cassert>
#include <atomic>
#include <omp.h> // ensure your build enables OpenMP

void Mesh::Populate(
    const std::vector<Vector3> &points,
    const std::vector<std::vector<int>> &faceToVert,
    const std::vector<std::vector<int>> &edgeToVert)
{
    Clear();

    _points = points;
    _faceToVert = faceToVert;
    _edgeToVert = edgeToVert;

    const int V = static_cast<int>(_points.size());
    const int F = static_cast<int>(_faceToVert.size());
    const int E = static_cast<int>(_edgeToVert.size());

    _faceToEdge.clear(); _faceToEdge.resize(F);
    _edgeToFace.clear(); _edgeToFace.resize(E);
    _vertToEdge.clear(); _vertToEdge.resize(V);
    _vertToFace.clear(); _vertToFace.resize(V);

    // Build heToEdge map (sequential) -- no changes
    std::unordered_map<uint64_t, int> heToEdge;
    heToEdge.reserve(static_cast<size_t>(E) * 2);
    for (int edge = 0; edge < E; ++edge)
    {
        const auto &ev = _edgeToVert[edge];
        assert(ev.size() == 2);
        int a = ev[0], b = ev[1];
        if (a > b) std::swap(a, b);
        uint64_t key = (static_cast<uint64_t>(static_cast<uint32_t>(a)) << 32)
                     | static_cast<uint32_t>(b);
        heToEdge[key] = edge;
    }

    // --- Counting phase (parallel, using Thread-Local Storage) ---
    // These are the final, global counters (no longer atomic)
    std::vector<int> vEdgeCount(V, 0);
    std::vector<int> vFaceCount(V, 0);
    std::vector<int> eFaceCount(E, 0);
    std::vector<int> faceEdgeCount(F, 0); // deterministic

    #pragma omp parallel
    {
        // --- Thread-Local Storage (TLS) ---
        // Each thread gets its own private (non-atomic) counter.
        // This completely eliminates atomic contention during counting.
        std::vector<int> tls_vEdge(V, 0);
        std::vector<int> tls_vFace(V, 0);
        std::vector<int> tls_eFace(E, 0);

        // Count edges per vertex (parallel)
        #pragma omp for schedule(static) nowait
        for (int edge = 0; edge < E; ++edge)
        {
            const auto &e = _edgeToVert[edge];
            tls_vEdge[e[0]]++; // No atomics!
            tls_vEdge[e[1]]++;
        }

        // Count face incidences (parallel)
        #pragma omp for schedule(dynamic)
        for (int face = 0; face < F; ++face)
        {
            const auto &fv = _faceToVert[face];
            const int n = static_cast<int>(fv.size());

            // This is safe to write from multiple threads
            // because they will all write the same value.
            faceEdgeCount[face] = n;

            for (int i = 0; i < n; ++i)
            {
                int vi = fv[i];
                int vj = fv[(i + 1) % n];
                int a = vi, b = vj;
                if (a > b) std::swap(a, b);
                uint64_t key = (static_cast<uint64_t>(static_cast<uint32_t>(a)) << 32)
                             | static_cast<uint32_t>(b);

                auto it = heToEdge.find(key);
                assert(it != heToEdge.end());
                int eidx = it->second;

                tls_eFace[eidx]++; // No atomics!
                tls_vFace[vi]++;
            }
        }

        // --- Reduction Step ---
        // Sum all thread-local counts into the global counts.
        // Use critical sections for simplicity; this is fast as it
        // only runs once per thread, *after* all the counting.
        #pragma omp critical
        for (int i = 0; i < V; ++i) {
            vEdgeCount[i] += tls_vEdge[i];
            vFaceCount[i] += tls_vFace[i];
        }
        
        #pragma omp critical
        for (int i = 0; i < E; ++i) {
            eFaceCount[i] += tls_eFace[i];
        }

    } // --- end of parallel region ---


    // --- Resize inner vectors (now parallel) ---
    // We parallelize the V+E+F independent allocations.
    #pragma omp parallel for schedule(static)
    for (int v = 0; v < V; ++v) {
        if (vEdgeCount[v] > 0) _vertToEdge[v].assign(vEdgeCount[v], -1);
        if (vFaceCount[v] > 0) _vertToFace[v].assign(vFaceCount[v], -1);
    }
    
    #pragma omp parallel for schedule(static)
    for (int e = 0; e < E; ++e) {
        if (eFaceCount[e] > 0) _edgeToFace[e].assign(eFaceCount[e], -1);
    }
    
    #pragma omp parallel for schedule(static)
    for (int f = 0; f < F; ++f) {
        if (faceEdgeCount[f] > 0) _faceToEdge[f].assign(faceEdgeCount[f], -1);
    }

    // --- Fill phase (parallel) ---
    // This logic is unchanged. Using atomic cursors to write into the
    // pre-sized vectors is the correct approach for this data structure.
    std::vector<std::atomic<int>> vertEdgePos(V), vertFacePos(V), edgeFacePos(E);
    for (int i = 0; i < V; ++i) { vertEdgePos[i].store(0, std::memory_order_relaxed); vertFacePos[i].store(0, std::memory_order_relaxed); }
    for (int i = 0; i < E; ++i) edgeFacePos[i].store(0, std::memory_order_relaxed);

    // Fill vert->edge (parallel)
    #pragma omp parallel for schedule(static)
    for (int edge = 0; edge < E; ++edge)
    {
        const auto &e = _edgeToVert[edge];
        int a = e[0], b = e[1];

        int pa = vertEdgePos[a].fetch_add(1, std::memory_order_relaxed);
        _vertToEdge[a][pa] = edge;

        int pb = vertEdgePos[b].fetch_add(1, std::memory_order_relaxed);
        _vertToEdge[b][pb] = edge;
    }

    // Fill face->edge, vert->face, edge->face (parallel over faces)
    #pragma omp parallel for schedule(dynamic)
    for (int face = 0; face < F; ++face)
    {
        const auto &fv = _faceToVert[face];
        const int n = static_cast<int>(fv.size());

        for (int i = 0; i < n; ++i)
        {
            int vi = fv[i];
            int vj = fv[(i + 1) % n];
            int a = vi, b = vj;
            if (a > b) std::swap(a, b);
            uint64_t key = (static_cast<uint64_t>(static_cast<uint32_t>(a)) << 32)
                         | static_cast<uint32_t>(b);

            auto it = heToEdge.find(key);
            assert(it != heToEdge.end());
            int eidx = it->second;

            // face -> edge: we reserved faceEdgeCount[face] == n and now write at i
            _faceToEdge[face][i] = eidx;

            // vert -> face: write at pos fetched atomically
            int pv = vertFacePos[vi].fetch_add(1, std::memory_order_relaxed);
            _vertToFace[vi][pv] = face;

            // edge -> face
            int pe = edgeFacePos[eidx].fetch_add(1, std::memory_order_relaxed);
            _edgeToFace[eidx][pe] = face;
        }
    }

    // Post-processing & validation
    _CheckValidity();
    _SortOneRings();

    _isTriangulated = _AreFacesTris();
    _isQuadrangulated = _AreFacesQuads();
}


void
Mesh::_CheckValidity() const
{
    for (int vert = 0; vert < _vertToFace.size(); ++vert)
    {
        assert(_vertToFace[vert].size());
        assert(_vertToEdge[vert].size());
        assert(
            (_vertToFace[vert].size()   == _vertToEdge[vert].size()) || 
            (_vertToFace[vert].size()+1 == _vertToEdge[vert].size()));
    }

    for (int edge = 0; edge < _edgeToVert.size(); ++edge)
    {
        assert(_edgeToVert[edge].size() == 2);
        assert(
            (_edgeToFace[edge].size() == 2) ||
            (_edgeToFace[edge].size() == 1));
    }

    for (int face = 0; face < _faceToVert.size(); ++face) 
    {
        assert(_faceToVert[face].size() > 2);
        assert(_faceToVert[face].size() == _faceToEdge[face].size());
    }
}

void 
Mesh::Triangulate()
{
    if (_isTriangulated) return;

    std::vector<Vector3> points = _points;
    std::vector< std::vector<int> > faceVerts;
    std::vector< std::vector<int> > edgeVerts = _edgeToVert;

    std::map< std::pair<int,int>, int > heToEdge;
    for (int edge = 0; edge < edgeVerts.size(); ++edge)
    {
        const std::vector<int>& eVert = edgeVerts[edge];
        std::pair<int, int> he = std::make_pair(eVert[0], eVert[1]);
        heToEdge[he] = edge;
    }

    // Triangulate with fan of triangles
    for (int face = 0; face < GetFaceCount(); ++face) 
    {
        const std::vector<int> &fVert = GetFaceVerts(face);
        for (size_t j = 2; j < fVert.size(); ++j) 
        {
            std::vector<int> t(3);
            t[0] = fVert[0  ];
            t[1] = fVert[j-1];
            t[2] = fVert[j  ];
            faceVerts.push_back(t);

            for (size_t k = 0; k < 3; ++k)
            {
                int v0 = fVert[k];
                int v1 = fVert[(k+1)%3];

                std::pair<int, int> he = std::make_pair(v0, v1);
                if (he.first > he.second) std::swap(he.first, he.second);
                if (heToEdge.find(he) != heToEdge.end()) continue;

                heToEdge[he] = edgeVerts.size();
                edgeVerts.push_back(std::vector<int>(2));
                edgeVerts.back()[0] = he.first;
                edgeVerts.back()[1] = he.second;
            }
        }
    }

    Populate(points, faceVerts, edgeVerts);
}

void
Mesh::_SortOneRings()
{
    for (size_t vert = 0; vert < GetVertCount(); ++vert)
    {
        int firstFace = GetVertFaces(vert)[0];
        if (IsVertBoundary(vert))
        {
            const std::vector<int>& vFaces = GetVertFaces(vert);
            for (int i = 0; i < vFaces.size(); ++i)
            {
                int face = vFaces[i];
                int index = GetVertIndexInFace(face, vert);
                int edge = GetFaceEdges(face)[index];
                if (!IsEdgeBoundary(edge)) continue;
                firstFace = face;
                break;
            } 
        }

        std::vector<int> vFaces, vEdges;
        vFaces.reserve(GetVertFaces(vert).size());
        vEdges.reserve(GetVertEdges(vert).size());

        int face = firstFace;
        do
        {
            int index = GetVertIndexInFace(face, vert);
            int edge = GetFaceEdges(face)[index];

            vFaces.push_back(face);
            vEdges.push_back(edge);

            if (index == 0) index = GetFaceEdges(face).size();
            edge = GetFaceEdges(face)[--index];

            if (IsEdgeBoundary(edge)) 
            {
                vEdges.push_back(edge);
                break;
            }
            
            index = GetFaceIndexInEdge(edge, face);
            face = GetEdgeFaces(edge)[(index+1)%2];
        }
        while (face != firstFace);

        if (IsVertBoundary(vert)) assert(vFaces.size()+1 == vEdges.size());
        else                      assert(vFaces.size()   == vEdges.size());

        _vertToFace[vert] = vFaces;
        _vertToEdge[vert] = vEdges;
    }
}

bool
Mesh::_AreFacesTris() const
{
    for (int face = 0; face < _faceToVert.size(); ++face) 
    {
        if (_faceToVert[face].size() != 3) return false;
    }
    return true;
}

bool 
Mesh::_AreFacesQuads() const
{
    for (int face = 0; face < _faceToVert.size(); ++face) 
    {
        if (_faceToVert[face].size() != 4) return false;
    }
    return true;
}
