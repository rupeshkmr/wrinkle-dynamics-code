#include "viz.hpp"
#include "polyscope/surface_mesh.h"

void VizCallback::updateAnimation() {
  if (m_Mesh.verts.size() >= m_frame) {
    polyscope::registerSurfaceMesh("Mesh", m_Mesh.verts[m_frame],
                                   m_Mesh.faces[m_frame]);
  } else {
    m_isPlaying = false;
    m_frame = m_Mesh.verts.size() - 1;
  }
}

void VizCallback::updateMesh(const Eigen::Matrix<double, Eigen::Dynamic, 3> &V,
                             const Eigen::Matrix<int, Eigen::Dynamic, 3> &F) {

  m_Mesh.faces.push_back(F);
  m_Mesh.verts.push_back(V);
}