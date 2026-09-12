# Dynamic Wrinkling on Coarsely Meshed Cloth

Implementation for the paper:  
**Dynamic Wrinkling on Coarsely Meshed Cloth**  
R. Kumar, S. Khurana, and R. Narain  
*Computer Graphics Forum (Proc. EUROGRAPHICS)*, 2026.  
DOI: [10.1111/cgf.70564](https://doi.org/10.1111/cgf.70564)

---

## Overview

This repository provides an implementation of our dynamic wrinkling simulation framework, which extends tension-field wrinkles (TFW) to capture high-frequency dynamic wrinkling effects efficiently on coarse cloth meshes.

---

## Dependencies

### System Requirements & Libraries
- **C++17** compliant compiler (GCC 9+, Clang 10+)
- **CMake** (version 3.16 or higher)
- **Eigen3**
- **SuiteSparse**
- **fcl** (Flexible Collision Library)
- **glm**
- **nlohmann_json**

On Ubuntu/Debian, the core dependencies can be installed via:
```bash
sudo apt-get update
sudo apt-get install build-essential cmake \
    libsuitesparse-dev libfcl-dev libglm-dev nlohmann-json3-dev libeigen3-dev
```

### Bundled & Automatically Managed Dependencies
- **Boost** (headers and filesystem library) is bundled under `vendor/boost/` so system Boost packages are not required.
- **libigl** is fetched automatically via CMake `FetchContent` during the initial configuration.
- Third-party modules for geometry processing, optimization, and visualization (`polyscope`, `TinyAD`, `LBFGSpp`, `SecondFundamentalForm`, `MeshLib`, `halfedge`) are included in `vendor/`.

---

## Building

Configure and build using CMake:

```bash
# Configure with visualization (Polyscope GUI) enabled
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DARGUS_ENABLE_VISUALIZATION=ON

# Build the targets
cmake --build build -j$(nproc)
```

The main simulation executable will be produced at `build/bin/wrinkle_dynamics` (and the shared library at `build/lib/libargus.so`).

### Build Options
- `-DARGUS_ENABLE_VISUALIZATION=ON/OFF` (default `OFF`): Enable Polyscope GUI support.
- `-DARGUS_CHECKPOINT=ON/OFF` (default `OFF`): Enable simulation state checkpoint dumps (positions, velocities, wrinkle amplitudes, and frequencies).

---

## Running Simulations

The driver executable takes a JSON configuration file and execution flags:

```bash
# Run headless (no GUI display needed)
./build/bin/wrinkle_dynamics configs/<config_file>.json --nogui

# Run with interactive Polyscope GUI
./build/bin/wrinkle_dynamics configs/<config_file>.json --gui

# Run and save checkpoints for post-processing/analysis
./build/bin/wrinkle_dynamics configs/<config_file>.json --nogui --checkpoint
```

### Example Configurations
All sample configs are located in the `configs/` directory.

- **Cloth Drape over Sphere (16x16 coarse mesh)**:
  ```bash
  ./build/bin/wrinkle_dynamics configs/1773844743.json --gui
  ```
- **Two-Edges Clamped Square Cloth**:
  ```bash
  ./build/bin/wrinkle_dynamics configs/1779268717.json --gui
  ```
- **Rotating Cylinder Twist**:
  ```bash
  ./build/bin/wrinkle_dynamics configs/1776251448.json --gui
  ```
- **Rotating Arm Model**:
  ```bash
  ./build/bin/wrinkle_dynamics configs/1780661372.json --gui
  ```
- **Fast Test / Comparison Run (50 frames)**:
  ```bash
  ./build/bin/wrinkle_dynamics configs/test_compare.json --nogui
  ```

---

## Output Checkpoints

When running with `--checkpoint`, output files are saved into `checkpoints/<scenario_name>/`. The state dumps include:
- `positions`: Nodal 3D coordinates per frame
- `velocities`: Vertex velocities per frame
- `amplitudes`: Wrinkle amplitude values per face
- `dphisPerFace`: Wrinkle direction and spatial frequency vectors per face
- `energies`: System energy components per frame

---

## Citation

```bibtex
@article{kumar2026dynamic,
  author    = {Kumar, R. and Khurana, S. and Narain, R.},
  title     = {Dynamic Wrinkling on Coarsely Meshed Cloth},
  journal   = {Computer Graphics Forum},
  year      = {2026},
  doi       = {10.1111/cgf.70564}
}
```
