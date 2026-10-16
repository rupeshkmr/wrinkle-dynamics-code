#pragma once
#include "commonHeader.hpp"


bool readFaces(std::string &facesFile, Eigen::MatrixXi& triF) 
{
    std::ifstream file(facesFile);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open faces file: " << facesFile << std::endl;
        return false;
    }

    std::vector<int> face_indices;
    int num_faces = 0;
    std::string line;

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        std::stringstream ss(line);
        std::string cell;
        int count = 0;
        
        while (std::getline(ss, cell, ',')) {
            try {
                face_indices.push_back(std::stoi(cell));
                count++;
            } catch (const std::exception& e) {
                std::cerr << "Error reading face index: " << e.what() << std::endl;
                return false;
            }
        }
        
        if (count != 3) {
            std::cerr << "Error: Face line does not contain 3 indices." << std::endl;
            // Depending on strictness, you might return false or skip.
            // Returning false for strict compliance.
            return false;
        }
        num_faces++;
    }

    if (num_faces == 0) {
        triF.resize(0, 3);
        return true;
    }

    // Resize and map the flattened vector to the Eigen matrix
    triF.resize(num_faces, 3);
    for (int i = 0; i < num_faces; ++i) {
        triF(i, 0) = face_indices[i * 3];
        triF(i, 1) = face_indices[i * 3 + 1];
        triF(i, 2) = face_indices[i * 3 + 2];
    }
    
    //std::cout << "Successfully read " << num_faces << " faces." << std::endl;
    return true;
}

bool readVertices(std::string &verticesFile, Eigen::MatrixXd &triVframe)
{
    std::ifstream file(verticesFile);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open vertices file: " << verticesFile << std::endl;
        return false;
    }

    std::vector<std::vector<double>> coords;
    coords.resize(3);
    
    int num_lines = 0;
    std::string line;
    int vertex_count = -1;

    while (std::getline(file, line) && num_lines < 3) {
        if (line.empty()) continue;
        
        std::stringstream ss(line);
        std::string cell;
        int current_count = 0;
        
        while (std::getline(ss, cell, ',')) {
            try {
                // Remove potential carriage return (\r) if reading files created on Windows
                if (!cell.empty() && cell.back() == '\r') {
                    cell.pop_back();
                }
                coords[num_lines].push_back(std::stod(cell));
                current_count++;
            } catch (const std::exception& e) {
                std::cerr << "Error reading coordinate: " << e.what() << std::endl;
                return false;
            }
        }
        
        if (vertex_count == -1) {
            vertex_count = current_count;
        } else if (vertex_count != current_count) {
            std::cerr << "Error: Coordinate lines have inconsistent vertex counts." << std::endl;
            return false;
        }
        
        num_lines++;
    }

    if (num_lines != 3) {
        std::cerr << "Error: Vertices file must contain exactly 3 coordinate lines (X, Y, Z)." << std::endl;
        return false;
    }

    if (vertex_count == 0) {
        triVframe.resize(0, 3);
        return true;
    }

    triVframe.resize(vertex_count, 3); 
    for (int j = 0; j < vertex_count; ++j) {
        triVframe(j, 0) = coords[0][j]; // X
        triVframe(j, 1) = coords[1][j]; // Y
        triVframe(j, 2) = coords[2][j]; // Z
    }
    
    //std::cout << "Successfully read " << vertex_count << " vertices." << std::endl;
    return true;
}

bool loadFaceOmegas(const std::string& filePath, const int& nfaces, Eigen::MatrixXd& faceOmegaMat)
{
    std::ifstream file(filePath);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open face omegas file: " << filePath << std::endl;
        return false;
    }

    // Resize the output matrix to the standard format: Nfaces x 2 components
    faceOmegaMat.resize(nfaces, 2);
    std::string line;
    int row_count = 0;

    // We expect exactly two lines (one for each component)
    while (std::getline(file, line) && row_count < 2) 
    {
        if (line.empty()) continue;

        std::stringstream ss(line);
        std::string cell;
        int col = 0;
        
        // Read values for the current component across all faces
        while (std::getline(ss, cell, ',')) 
        {
            // If we've read more faces than expected, something is wrong
            if (col >= nfaces) {
                std::cerr << "Error: Row " << row_count + 1 << " has more than " << nfaces << " components." << std::endl;
                return false;
            }

            try {
                // Store the value:
                // faceOmegaMat(face_index, component_index)
                // Current face index is 'col'. Current component index is 'row_count'.
                faceOmegaMat(col, row_count) = std::stod(cell);
            } catch (const std::exception& e) {
                std::cerr << "Error converting face omega value in row " << row_count + 1 << ", column " << col + 1 << ": " << e.what() << std::endl;
                return false;
            }
            col++;
        }
        
        // Check if the number of faces read matches the expected count
        if (col != nfaces) {
            std::cerr << "Error: Row " << row_count + 1 << " read " << col << " components, expected " << nfaces << "." << std::endl;
            return false;
        }

        row_count++;
    }

    // Final check to ensure exactly 2 rows were processed
    if (row_count != 2) {
        std::cerr << "Error: Expected 2 rows in the file but read " << row_count << "." << std::endl;
        return false;
    }

    //std::cout << "Successfully read and transposed " << nfaces << " face omegas." << std::endl;
    return true;
}


#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iostream>
#include <regex>

// Helper: read Only UV faces from OBJ into Eigen matrix.
// - Must match the triangulation of readObjVerticesAndFaces.
static bool getUVFaces(
    const std::string& objPath,
    Eigen::MatrixXi& uvFaces
) {
    std::ifstream in(objPath);
    if (!in.is_open()) {
        std::cerr << "getUVFaces: failed to open " << objPath << std::endl;
        return false;
    }

    int vt_count = 0; 
    std::vector<Eigen::RowVector3i> faces;

    std::string line;
    std::regex slash_re(R"(/)");

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        size_t p = 0;
        while (p < line.size() && isspace((unsigned char)line[p])) ++p;
        if (p >= line.size()) continue;
        if (line[p] == '#') continue;

        // Track vt count for negative indices
        if (line.compare(p, 3, "vt ") == 0 || line.compare(p, 3, "vt\t") == 0) {
            vt_count++;
        }
        else if (line.compare(p, 2, "f ") == 0 || line.compare(p, 2, "f\t") == 0) {
            std::istringstream iss(line.substr(p + 1));
            std::vector<int> idxs;
            std::string tok;
            while (iss >> tok) {
                // split by '/'
                std::vector<std::string> parts;
                std::sregex_token_iterator it(tok.begin(), tok.end(), slash_re, -1);
                std::sregex_token_iterator end;
                for (; it != end; ++it) parts.push_back(it->str());

                int vt_idx = 0; 
                // v/vt/vn -> parts[1] is vt
                if (parts.size() >= 2 && !parts[1].empty()) {
                     try {
                        vt_idx = std::stoi(parts[1]);
                    } catch (...) {
                        vt_idx = 0;
                    }
                }
                
                if (vt_idx != 0) {
                    if (vt_idx < 0) {
                         vt_idx = vt_count + vt_idx + 1;
                    }
                    vt_idx -= 1; // 0-based
                }
                
                idxs.push_back(vt_idx);
            }

            // triangulate if face has >3 vertices using fan triangulation
            if (idxs.size() >= 3) {
                for (size_t t = 1; t + 1 < idxs.size(); ++t) {
                    faces.emplace_back(idxs[0], idxs[t], idxs[t+1]);
                }
            }
        }
    }

    if (faces.empty()) {
        uvFaces.resize(0, 3);
        return true;
    }

    uvFaces.resize((int)faces.size(), 3);
    for (size_t fi = 0; fi < faces.size(); ++fi) {
        uvFaces.row(fi) = faces[fi];
    }
    
    return true;
}

// Helper: read OBJ vertices and triangular faces into Eigen matrices.
// - Supports face lines like "f v1 v2 v3" or "f v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3".
// - Produces 0-based indices in triF.
static bool readObjVerticesAndFaces(
    const std::string& objPath,
    const int& uvPresent,
    Eigen::MatrixXd& V,       // Nx3
    Eigen::MatrixXi& F,        // Mx3 (0-based indices)
    Eigen::MatrixXd& uvCoordinates
) {
    std::ifstream in(objPath);
    if (!in.is_open()) {
        std::cerr << "readObjVerticesAndFaces: failed to open " << objPath << std::endl;
        return false;
    }

    std::vector<Eigen::RowVector3d> verts;
    std::vector<Eigen::RowVector2d> uvs;  // Store UV coordinates
    std::vector<Eigen::RowVector3i> faces;

    std::string line;
    // Regex to split face tokens like "v", "v/vt", "v//vn", "v/vt/vn"
    std::regex slash_re(R"(/)");
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        // Trim leading spaces
        size_t p = 0;
        while (p < line.size() && isspace((unsigned char)line[p])) ++p;
        if (p >= line.size()) continue;
        if (line[p] == '#') continue;

        // vertex
        if (line.compare(p, 2, "v ") == 0 || line.compare(p, 2, "v\t") == 0) {
            std::istringstream iss(line.substr(p + 1));
            double x, y, z;
            if (!(iss >> x >> y >> z)) {
                // maybe only 2 coords, handle gracefully by adding 0 for z
                iss.clear();
                iss.str(line.substr(p + 1));
                if (!(iss >> x >> y)) continue;
                z = 0.0;
            }
            verts.emplace_back(x, y, z);
        }
        // texture coordinate (vt u v [w])
        else if (uvPresent == 1 && (line.compare(p, 3, "vt ") == 0 || line.compare(p, 3, "vt\t") == 0)) {
            std::istringstream iss(line.substr(p + 2));
            double u, v;
            if (!(iss >> u >> v)) continue;
            uvs.emplace_back(u, v);
        }
        // face
        else if (line.compare(p, 2, "f ") == 0 || line.compare(p, 2, "f\t") == 0) {
            std::istringstream iss(line.substr(p + 1));
            std::vector<int> idxs;
            std::string tok;
            while (iss >> tok) {
                // split by '/'
                std::vector<std::string> parts;
                std::sregex_token_iterator it(tok.begin(), tok.end(), slash_re, -1);
                std::sregex_token_iterator end;
                for (; it != end; ++it) parts.push_back(it->str());

                if (parts.size() == 0) continue;
                // the first part is the vertex index (1-based). OBJ uses 1-based indexing and negative indices allowed.
                int vi = 0;
                try {
                    vi = std::stoi(parts[0]);
                } catch (...) {
                    vi = 0;
                }
                // convert negative indices: -1 means last vertex
                if (vi < 0) vi = (int)verts.size() + vi + 1; // still 1-based after this
                // ensure within bounds (we'll check later)
                idxs.push_back(vi - 1); // convert to 0-based
            }
            // triangulate if face has >3 vertices using fan triangulation
            if (idxs.size() >= 3) {
                for (size_t t = 1; t + 1 < idxs.size(); ++t) {
                    faces.emplace_back(idxs[0], idxs[t], idxs[t+1]);
                }
            }
        }
    }

    // Convert vectors to Eigen matrices
    if (verts.empty()) {
        std::cerr << "readObjVerticesAndFaces: no vertices found in " << objPath << std::endl;
        return false;
    }
    V.resize((int)verts.size(), 3);
    for (size_t vi = 0; vi < verts.size(); ++vi) {
        V.row((int)vi) = verts[vi];
    }

    // Convert UV coordinates if present
    // Map (u, v) from OBJ to (u, 0, v) in uvCoordinates
    if (!uvs.empty() && uvPresent) {
        uvCoordinates.resize((int)uvs.size(), 3);
        for (size_t ui = 0; ui < uvs.size(); ++ui) {
            uvCoordinates(ui, 0) = uvs[ui][0];  // u -> first coordinate
            uvCoordinates(ui, 1) = 0.0;         // 0 -> second coordinate
            uvCoordinates(ui, 2) = uvs[ui][1];  // v -> third coordinate
        }
    } else {
        uvCoordinates.resize(0, 3);
    }

    if (faces.empty()) {
        std::cerr << "readObjVerticesAndFaces: no faces found in " << objPath << std::endl;
        return false;
    }
    F.resize((int)faces.size(), 3);
    for (size_t fi = 0; fi < faces.size(); ++fi) {
        // Basic bounds checking
        int a = faces[fi][0], b = faces[fi][1], c = faces[fi][2];
        if (a < 0 || b < 0 || c < 0 || a >= (int)verts.size() || b >= (int)verts.size() || c >= (int)verts.size()) {
            std::cerr << "readObjVerticesAndFaces: face index out of range in " << objPath << " face #" << fi << std::endl;
            return false;
        }
        F(fi, 0) = a;
        F(fi, 1) = b;
        F(fi, 2) = c;
    }

    return true;
}

void readSeamMap(
    const std::string& seamFile,
    std::map<int, int>& seamMap,
    int uvpresent,
    const Eigen::MatrixXi& F,     // required if uvpresent != 0
    const Eigen::MatrixXi& FUV    // required if uvpresent != 0
) {
    seamMap.clear();

    // ------------------------------------------------------------
    // Case 1: UVs are present → build seam map automatically
    // ------------------------------------------------------------
    if (uvpresent) {
        std::cout << "UVs detected — building seam map from (F, FUV)" << std::endl;

        // Map (orig_vertex, uv_vertex) → cut_vertex
        std::map<std::pair<int,int>, int> pairToCut;
        int nextCutIndex = 0;

        for (int f = 0; f < F.rows(); ++f) {
            for (int c = 0; c < 3; ++c) {
                int v_orig = F(f, c);
                int v_uv   = FUV(f, c);

                std::pair<int,int> key(v_orig, v_uv);

                if (pairToCut.find(key) == pairToCut.end()) {
                    pairToCut[key] = nextCutIndex;
                    seamMap[nextCutIndex] = v_orig;
                    nextCutIndex++;
                }
            }
        }

        std::cout << "Constructed seam map with "
                  << seamMap.size()
                  << " cut vertices" << std::endl;

        return;
    }
    else {
        // ------------------------------------------------------------
        // Case 2: No UVs → read seam map from file (original behavior)
        // ------------------------------------------------------------
        std::cout << "Attempting to read seam map from: " << seamFile << std::endl;

        std::ifstream in(seamFile);
        if (!in.is_open()) {
            throw std::runtime_error("Failed to open seam map file: " + seamFile);
        }

        std::cout << "File opened successfully" << std::endl;

        std::string line;
        int lineNo = 0;
        int parsedCount = 0;

        while (std::getline(in, line)) {
            ++lineNo;

            size_t end = line.find_last_not_of(" \t\r\n");
            if (end != std::string::npos)
                line = line.substr(0, end + 1);
            else
                continue;

            if (line.empty())
                continue;

            std::istringstream iss(line);
            int a, b;

            if (!(iss >> a >> b)) {
                std::cerr << "Warning: Failed to parse line "
                        << lineNo << ": '" << line << "'" << std::endl;
                continue;
            }

            seamMap[a] = b;
            parsedCount++;

            if (parsedCount <= 5) {
                std::cout << "  Parsed: " << a << " -> " << b << std::endl;
            }
        }

        std::cout << "Total lines read: " << lineNo << std::endl;
        std::cout << "Loaded " << seamMap.size()
                << " seam vertex pairs from file" << std::endl;
    }
}

void buildUVto3DMapping(
    const Eigen::MatrixXi& F,
    const Eigen::MatrixXi& F_uv,
    int num3DVerts,
    int numUVVerts,
    std::vector<int>& uvTo3D,
    std::vector<std::vector<int>>& threeDtoUV
    )
{
    uvTo3D.assign(numUVVerts, -1);
    threeDtoUV.assign(num3DVerts, std::vector<int>());

    for (int i = 0; i < F.rows(); ++i) {
        for (int j = 0; j < 3; ++j) {
            int v3d = F(i, j);
            int vuv = F_uv(i, j);

            if (vuv >= 0 && vuv < numUVVerts) {
                uvTo3D[vuv] = v3d;
            }
            
            if (v3d >= 0 && v3d < num3DVerts) {
                // Check if already added
                bool found = false;
                for (int u : threeDtoUV[v3d]) {
                    if (u == vuv) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    threeDtoUV[v3d].push_back(vuv);
                }
            }
        }
    }
}

void map3DVertexDataToUV(
    const Eigen::VectorXd& data3D,
    const std::vector<int>& uvTo3D,
    Eigen::VectorXd& dataUV
) {
    int numUV = uvTo3D.size();
    dataUV.resize(numUV);
    for (int i = 0; i < numUV; ++i) {
        int v3d = uvTo3D[i];
        if (v3d >= 0 && v3d < data3D.size()) {
            dataUV[i] = data3D[v3d];
        } else {
             dataUV[i] = 0.0;
        }
    }
}

