# CUDA Setup (Ubuntu/Linux)

## 1. Install NVIDIA driver

Check GPU:

```bash
lspci | grep -i nvidia
```

Install recommended driver (example):

```bash
sudo ubuntu-drivers autoinstall
sudo reboot
```

Verify:

```bash
nvidia-smi
```

## 2. Install CUDA Toolkit

Option A: Ubuntu package (quick)

```bash
sudo apt update
sudo apt install -y nvidia-cuda-toolkit
```

Option B: NVIDIA official repo/toolkit (newer; recommended for recent GPUs)

Follow: https://developer.nvidia.com/cuda-downloads

Verify:

```bash
nvcc --version
```

## 3. Build integration in this repo

This repository includes:

- `include/ComplexLoop/ComplexLoopZuenkoCUDA.h`
- `src/ComplexLoop/ComplexLoopZuenkoCUDA.cu`

`app/WrinkleInterpolationCli/CMakeLists.txt` is wired to compile/link CUDA for CLI when CMake finds a CUDA compiler (`nvcc`).

## 4. Configure and build

```bash
cd build
cmake ..
make -j
```

If CMake cannot find CUDA, pass:

```bash
cmake .. -DCUDAToolkit_ROOT=/usr/local/cuda
```

If needed, also set the CUDA compiler explicitly:

```bash
cmake .. -DCUDAToolkit_ROOT=/usr/local/cuda -DCMAKE_CUDA_COMPILER=/usr/local/cuda/bin/nvcc
```

## 5. Runtime backend selection

`app/WrinkleInterpolationCli/main.cpp` has a selector:

```cpp
constexpr UpsamplingBackend kUpsamplingBackend = UpsamplingBackend::CUDA;
```

Switch to CPU by setting:

```cpp
constexpr UpsamplingBackend kUpsamplingBackend = UpsamplingBackend::CPU;
```
