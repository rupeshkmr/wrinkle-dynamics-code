#ifndef VIZ_HPP
#define VIZ_HPP
#include "polyscope/messages.h"
#include "polyscope/pick.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/view.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
struct meshData {
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, 3>> verts;
  std::vector<Eigen::Matrix<int, Eigen::Dynamic, 3>> faces;
};

class VizCallback {
  int m_frame;
  int m_totalFrames;
  bool m_isPlaying;
  struct meshData m_Mesh;

public:
  VizCallback(int total_frames) : m_totalFrames(total_frames) {
    m_isPlaying = false;
    m_frame = 0;
  }

  int sanitizeFrame(int f) {
    if (f < 0)
      return m_totalFrames - 1;
    if (f >= m_totalFrames)
      return 0;
    return f;
  }
  void setFrame(int f) { m_frame = sanitizeFrame(f); }
  void incrementFrame() { m_frame = sanitizeFrame(m_frame + 1); }
  void decrementFrame() { m_frame = sanitizeFrame(m_frame - 1); }
  void callback() {
    // 1. Wrap in a dedicated UI window
    ImGui::Begin("Simulation Control");

    // 2. Play/Pause Logic
    if (m_isPlaying) {
      if (ImGui::Button("Pause"))
        m_isPlaying = false;
      else {
        incrementFrame();
        updateAnimation();
      }
    } else {
      if (ImGui::Button("Play")) {
        m_isPlaying = true;
        incrementFrame();
        updateAnimation();
      }
    }

    ImGui::SameLine();

    // 3. Manual Step Buttons
    if (ImGui::Button("<<")) {
      decrementFrame();
      updateAnimation();
    }
    ImGui::SameLine();
    if (ImGui::Button(">>")) {
      incrementFrame();
      updateAnimation();
    }

    // 4. Frame Slider (Very helpful for shell research)
    if (ImGui::SliderInt("Frame", &m_frame, 0, m_totalFrames - 1)) {
      updateAnimation();
    }

    // 5. The "Playing" engine
    if (m_isPlaying) {
      incrementFrame();
      updateAnimation();
    }

    ImGui::End();
  }
  void updateAnimation();
  void updateMesh(const Eigen::Matrix<double, Eigen::Dynamic, 3> &V,
                  const Eigen::Matrix<int, Eigen::Dynamic, 3> &F);
};
#endif