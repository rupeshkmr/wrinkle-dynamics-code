# Dynamic Wrinkling on Coarsely Meshed Cloth

Implementation for the paper:  
**Dynamic Wrinkling on Coarsely Meshed Cloth**  
R. Kumar, S. Khurana, and R. Narain  
*Computer Graphics Forum (Proc. EUROGRAPHICS)*, 2026.  
DOI: [10.1111/cgf.70564](https://doi.org/10.1111/cgf.70564)

---

## Overview

This repository provides an implementation of our dynamic wrinkling simulation framework, which extends tension-field wrinkles (TFW) to capture high-frequency dynamic wrinkling effects efficiently on coarse cloth meshes.

The repository also contains the wrinkle upsampler (`upsampling/`, based on [Complex Wrinkle Field Evolution](https://github.com/zhenchen-jay/Complex-Wrinkle-Field)), which turns the coarse simulation output (mesh, per-face wrinkle amplitudes and frequencies) into high-resolution wrinkled meshes. A single script runs the whole pipeline, simulation followed by upsampling (see [Full Pipeline](#full-pipeline-simulate-and-upsample)).

---

## Getting the Code

Clone the repository with the `--recursive` flag to automatically fetch all required submodules in `vendor/`:

```bash
git clone --recursive https://github.com/rupeshkmr/wrinkle-dynamics-code.git
cd wrinkle-dynamics-code
```

If you have already cloned the repository without `--recursive`, initialize and update the submodules manually:

```bash
git submodule update --init --recursive
```

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

The upsampler (optional, `-DARGUS_BUILD_UPSAMPLER=ON`) additionally needs **OpenMP** and **SuiteSparse** from the system; its other dependencies (polyscope, libigl, TBB, Spectra, geometry-central, CLI11, LBFGSpp, Discregrid) are fetched automatically during the build. Running the upsampling scripts needs **Python 3** (standard library only).

### Bundled & Automatically Managed Dependencies
- **Boost** (headers and filesystem library) is bundled under `vendor/boost/` so system Boost packages are not required.
- **libigl** is fetched automatically via CMake `FetchContent` during the initial configuration.
- Third-party modules for geometry processing, optimization, and visualization (`polyscope`, `TinyAD`, `LBFGSpp`, `SecondFundamentalForm`, `MeshLib`, `halfedge`, `tmd`, `Discregrid`, `dlib`) are managed in `vendor/`.

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

To also build the upsampler (produced at `build/bin/upsample`; this also turns on `ARGUS_CHECKPOINT`, which the upsampler needs):

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DARGUS_BUILD_UPSAMPLER=ON
cmake --build build -j$(nproc)
```

With CMake 4.x, add `-DCMAKE_POLICY_VERSION_MINIMUM=3.5` to the configure command (some vendored dependencies declare an older minimum CMake version).

### Build Options
- `-DARGUS_ENABLE_VISUALIZATION=ON/OFF` (default `OFF`): Enable Polyscope GUI support.
- `-DARGUS_CHECKPOINT=ON/OFF` (default `OFF`): Enable simulation state checkpoint dumps (positions, velocities, wrinkle amplitudes, and frequencies); required for `--checkpoint`.
- `-DARGUS_BUILD_UPSAMPLER=ON/OFF` (default `OFF`): Build the wrinkle upsampler in `upsampling/` as `build/bin/upsample`.

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

- **Cloth Drape over Sphere (16x16 coarse mesh, symplectic Euler)**:
  ```bash
  ./build/bin/wrinkle_dynamics configs/sphere_drape.json --gui
  ```
- **Two Adjacent Edges Clamped Square Cloth (32x32 coarse mesh)**:
  ```bash
  ./build/bin/wrinkle_dynamics configs/cloth_drape.json --gui
  ```
- **Rotating Cylinder Twist**:
  ```bash
  ./build/bin/wrinkle_dynamics configs/cylinder_twist.json --gui
  ```
- **Rotating Arm Model**:
  ```bash
  ./build/bin/wrinkle_dynamics configs/sleeve.json --gui
  ```
- **Sphere Drape with Hybrid Stretching (stretching-feedback ablation, `tfwstvkhybridstretch`; TFW bending retained)**:
  ```bash
  ./build/bin/wrinkle_dynamics configs/1790250839.json --gui
  ```
Each example's solver is selected by its `timestepper` field. For symplectic runs, the per-face frequency
solve currently uses a hard-coded tolerance and iteration limit in
`Sheet::optimizeParallelEnergies`; the generic `optimizer_*` fields do not set
that solve's stopping criteria.

---

## Output Checkpoints

When running with `--checkpoint`, output files are saved into `checkpoints/<scenario_name>/`. The state dumps include:
- `positions`: Nodal 3D coordinates per frame
- `velocities`: Vertex velocities per frame
- `amplitudes`: Wrinkle amplitude values per face
- `dphisPerFace`: Wrinkle direction and spatial frequency vectors per face
- `energies`: System energy components per frame

---

## Upsampling

`scripts/upsample.py` upsamples one or more checkpoints into wrinkled high-resolution meshes:

```bash
# upsample every frame
python3 scripts/upsample.py checkpoints/<checkpoint_id>

# keep every 20th frame (and the last one); several checkpoints run in parallel
python3 scripts/upsample.py checkpoints/<id1> checkpoints/<id2> --every 20
```

Results are written to `output/upsample/<checkpoint_id>/render/`:
- `wrinkledMesh/wrinkledMesh_<frame>.obj`: the upsampled wrinkled mesh
- `upsampledAmplitudes/`, `upsampledPhases/`: the upsampled wrinkle amplitude and phase per vertex

The upsampler needs a reference (rest/UV) mesh on which the wrinkle field is defined. It is chosen in this order:
1. `--ref <mesh.obj> --seam <seamVerts.txt>` given on the command line;
2. `<mesh>_reference.obj` and `<mesh>_seamVerts.txt` (or `_seamVert.txt`) next to the cloth mesh of the configuration, e.g. `examples/assets/arm/arm_1000verts_reference.obj` for the arm;
3. otherwise the rest mesh of the checkpoint without seams, which is what flat sheets use.

For the twisted cylinder with seams, pass `--ref examples/assets/cylinder_rotating/cylinder_16x16_withSeams_Reference.obj --seam examples/assets/cylinder_rotating/cylinder_16x16_withSeams_seamVert.txt`.

Other options: `--frames N` (only the first N frames), `--param key=value` (override an upsampler parameter, e.g. `--param upsampledTimes=2`), `--out <dir>`, `--bin <upsample executable>`, `--ram-gb` (memory budget for parallel runs; the upsampler keeps all upsampled frames in memory). Run `python3 scripts/upsample.py -h` for details.

### Full Pipeline (Simulate and Upsample)

```bash
scripts/run_pipeline.sh configs/<config_file>.json --every 20
```

This runs `wrinkle_dynamics` headless with checkpoints and then upsamples the new checkpoint. Checkpoint directories are named after the start time in seconds, so do not start two simulations within the same second (run several pipelines one after another, or a second apart). Any options after the config are passed to `scripts/upsample.py`. The build directory is found automatically (`build/`, or the newest `build/*/` containing both executables); set `BUILD_DIR` to choose one.

---

## Acknowledgements

- Portions of the code and simulation formulation build upon the reference implementation of **Wrinkled Tension Fields** ([zhenchen-jay/WrinkledTensionFields](https://github.com/zhenchen-jay/WrinkledTensionFields)). We gratefully acknowledge the authors for open-sourcing their codebase.
- The upsampler in `upsampling/` is adapted from **Complex Wrinkle Field Evolution** by Chen et al. ([zhenchen-jay/Complex-Wrinkle-Field](https://github.com/zhenchen-jay/Complex-Wrinkle-Field), MIT license, see `upsampling/LICENSE`).
- Several libraries in the `vendor/` directory (such as `halfedge`, `MeshLib`, `SecondFundamentalForm`, `tmd`, and `LBFGSpp`) are modified adaptations derived directly from or tailored for the [WrinkledTensionFields](https://github.com/zhenchen-jay/WrinkledTensionFields) framework.

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
