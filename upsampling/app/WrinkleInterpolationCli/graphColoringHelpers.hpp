#pragma once
#include <Eigen/Core>
#include <vector>
#include <set>
#include <algorithm>
#include <iostream>

// -----------------------------------------------------------------------------
// 1. Build vertex->vertex adjacency (for graph coloring)
// -----------------------------------------------------------------------------
inline void build_vertex_to_vertex_adjacency(
    const Eigen::MatrixXi& F,
    int nverts,
    std::vector<std::vector<int>>& vert_adj)
{
    std::vector<std::set<int>> temp(nverts);

    for (int f = 0; f < F.rows(); ++f) {
        for (int j = 0; j < 3; ++j) {
            int v0 = F(f, j);
            int v1 = F(f, (j + 1) % 3); // next vertex in the triangle

            if (v0 >= 0 && v0 < nverts && v1 >= 0 && v1 < nverts) {
                temp[v0].insert(v1);
                temp[v1].insert(v0);
            }
        }
    }

    vert_adj.assign(nverts, {});
    for (int v = 0; v < nverts; ++v) {
        vert_adj[v].assign(temp[v].begin(), temp[v].end());
    }
}

// -----------------------------------------------------------------------------
// 2. Build vertex->face adjacency (for energy / VBD)
// -----------------------------------------------------------------------------
inline void build_vertex_to_face_adjacency(
    const Eigen::MatrixXi& F,
    int nverts,
    std::vector<std::vector<int>>& vert_to_faces)
{
    vert_to_faces.clear();
    vert_to_faces.resize(nverts);

    for (int f_idx = 0; f_idx < F.rows(); ++f_idx) {
        int i = F(f_idx, 0);
        int j = F(f_idx, 1);
        int k = F(f_idx, 2);

        if (i >= 0 && i < nverts) vert_to_faces[i].push_back(f_idx);
        if (j >= 0 && j < nverts) vert_to_faces[j].push_back(f_idx);
        if (k >= 0 && k < nverts) vert_to_faces[k].push_back(f_idx);
    }

#ifndef NDEBUG
    // Sanity check: each (v, f_idx) must satisfy v ∈ {F(f_idx,0), F(f_idx,1), F(f_idx,2)}
    for (int v = 0; v < nverts; ++v) {
        for (int f_idx : vert_to_faces[v]) {
            int a = F(f_idx, 0);
            int b = F(f_idx, 1);
            int c = F(f_idx, 2);
            if (v != a && v != b && v != c) {
                std::cerr << "[build_vertex_to_face_adjacency] ERROR: vertex "
                          << v << " lists face " << f_idx
                          << " but F(f) = (" << a << "," << b << "," << c << ")\n";
            }
        }
    }
#endif
}

// -----------------------------------------------------------------------------
// 3. Maximum Cardinality Search ordering (for coloring)
// -----------------------------------------------------------------------------
inline std::vector<int> mcs_ordering(const std::vector<std::vector<int>>& graph)
{
    const int n = static_cast<int>(graph.size());
    std::vector<int> weight(n, 0);
    std::vector<int> order(n, -1);
    std::vector<char> selected(n, 0);

    // Fill order from the end: order[k] is the k-th vertex in MCS order
    for (int k = n - 1; k >= 0; --k) {
        int max_v = -1;
        int max_w = -1;

        // pick unselected vertex with maximum weight
        for (int v = 0; v < n; ++v) {
            if (!selected[v] && weight[v] >= max_w) {
                max_w = weight[v];
                max_v = v;
            }
        }

        if (max_v == -1) {
            std::cerr << "[mcs_ordering] ERROR: no unselected vertex found.\n";
            break;
        }

        order[k] = max_v;
        selected[max_v] = 1;

        // bump weights of unselected neighbors
        for (int nb : graph[max_v]) {
            if (!selected[nb]) {
                weight[nb] += 1;
            }
        }
    }

    return order;
}

// -----------------------------------------------------------------------------
// 4. Greedy coloring using a given ordering
// -----------------------------------------------------------------------------
inline std::vector<int> greedy_coloring_with_order(
    const std::vector<std::vector<int>>& graph,
    const std::vector<int>& order)
{
    const int n = static_cast<int>(graph.size());
    std::vector<int> color(n, -1);
    int maxColor = -1;

    // We need enough room for "forbidden" colors; worst-case n colors
    std::vector<char> forbidden(n, 0);

    for (int idx = 0; idx < n; ++idx) {
        int v = order[idx];

        // mark colors used by neighbors
        std::fill(forbidden.begin(), forbidden.end(), 0);
        for (int nb : graph[v]) {
            int c = color[nb];
            if (c >= 0) {
                forbidden[c] = 1;
            }
        }

        // pick the smallest non-forbidden color
        int c = 0;
        while (c <= maxColor && forbidden[c]) {
            ++c;
        }

        color[v] = c;
        if (c > maxColor) maxColor = c;
    }

    return color;
}

// -----------------------------------------------------------------------------
// 5. Public API: vertex coloring for your mesh (for GS batching)
// -----------------------------------------------------------------------------
inline void vertex_coloring_mcs(
    const Eigen::MatrixXi& F,
    int nverts,
    std::vector<std::vector<int>>& V_by_color)
{
    // 1. Build adjacency graph: vertex -> neighbors
    std::vector<std::vector<int>> vert_graph;
    build_vertex_to_vertex_adjacency(F, nverts, vert_graph);

    // 2. Compute MCS ordering
    std::vector<int> order = mcs_ordering(vert_graph);

    // 3. Greedy color in that ordering
    std::vector<int> vertex_colors = greedy_coloring_with_order(vert_graph, order);

    // 4. Sanity check: ensure coloring is valid (no adjacent same-colored verts)
    bool ok = true;
    for (int v = 0; v < nverts; ++v) {
        int c = vertex_colors[v];
        if (c < 0) {
            std::cerr << "[vertex_coloring_mcs] ERROR: vertex " << v << " uncolored.\n";
            ok = false;
        }
        for (int nb : vert_graph[v]) {
            if (vertex_colors[nb] == c) {
                std::cerr << "[vertex_coloring_mcs] ERROR: vertices " << v
                          << " and " << nb << " share color " << c << ".\n";
                ok = false;
            }
        }
    }

    if (!ok) {
        std::cerr << "[vertex_coloring_mcs] WARNING: coloring is invalid!\n";
        // you could std::abort() here if you want to hard-fail
    }

    // 5. Bucket vertices by color
    int num_colors = 0;
    for (int c : vertex_colors) {
        if (c + 1 > num_colors) num_colors = c + 1;
    }

    V_by_color.clear();
    V_by_color.resize(num_colors);
    for (int v = 0; v < nverts; ++v) {
        int c = vertex_colors[v];
        if (c >= 0 && c < num_colors) {
            V_by_color[c].push_back(v);
        }
    }

    std::cout << "MCS coloring complete. Num colors = " << num_colors << "\n";
    for (int c = 0; c < num_colors; ++c) {
        std::cout << "  Color " << c << ": " << V_by_color[c].size() << " vertices\n";
    }
}

// -----------------------------------------------------------------------------
// 6. Verify graph coloring against mesh faces
// -----------------------------------------------------------------------------
inline bool verify_graph_coloring(
    const Eigen::MatrixXi& F,                  // Face connectivity
    int nverts,                                // Total number of vertices
    const std::vector<std::vector<int>>& V_by_color // The output to verify
)
{
    std::cout << "--- Verifying Graph Coloring ---" << std::endl;

    // Flatten: vertex -> color
    std::vector<int> vertex_colors(nverts, -1);
    int colored_count = 0;

    for (size_t c = 0; c < V_by_color.size(); ++c) {
        for (int v : V_by_color[c]) {
            if (v < 0 || v >= nverts) {
                std::cerr << "[Error] Vertex index " << v << " is out of bounds (0-" << nverts - 1 << ")\n";
                return false;
            }
            if (vertex_colors[v] != -1) {
                std::cerr << "[Error] Vertex " << v << " is assigned to multiple color buckets (Color " 
                          << vertex_colors[v] << " and " << c << ")\n";
                return false;
            }
            vertex_colors[v] = static_cast<int>(c);
            colored_count++;
        }
    }

    if (colored_count != nverts) {
        std::cerr << "[Error] Coloring is incomplete. " << colored_count 
                  << " vertices colored out of " << nverts << ".\n";
        for (int v = 0; v < nverts; ++v) {
            if (vertex_colors[v] == -1) {
                std::cerr << "   - First missing vertex: " << v << "\n";
                break;
            }
        }
        return false;
    }

    // Check adjacency via faces
    bool valid = true;
    int conflict_count = 0;

    for (int f = 0; f < F.rows(); ++f) {
        for (int edge = 0; edge < 3; ++edge) {
            int v1 = F(f, edge);
            int v2 = F(f, (edge + 1) % 3);

            if (v1 < 0 || v1 >= nverts || v2 < 0 || v2 >= nverts) continue;

            if (vertex_colors[v1] == vertex_colors[v2]) {
                if (conflict_count < 5) {
                    std::cerr << "[Conflict] Neighboring vertices " << v1 << " and " << v2 
                              << " share Color " << vertex_colors[v1] << " (Face " << f << ")\n";
                }
                valid = false;
                conflict_count++;
            }
        }
    }

    if (!valid) {
        std::cerr << "[Failed] Found " << conflict_count << " edge conflicts.\n";
    } else {
        std::cout << "[Success] Graph coloring is valid.\n";
    }

    return valid;
}
