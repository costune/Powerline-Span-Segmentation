# PCGrid - Pointcloud Grid for Powerline Span Segmentation

A tool for extracting correct powerline span with input line and tower pointclouds, combining catenary model regression.


## Table of Contents

- [Introduction](#introduction)
- [Project Structure](#project-structure)
- [Installation](#installation)
- [Usage](#usage)
- [Examples](#examples)
- [BibTeX](#bibtex)

## Introduction

The project is to identify correct span infomation in airbone Lidar powerline pointclouds. With the input of line pcd and tower pcd with noise, the algorithm can output checked main-line tower positions and mian-line powerline pcd.

The algorithm support both `C++ binary file` and `Python binding with Pybind`. For Further infomation, refering the paper [Enhanced Automatic Span Segmentation of Airborne LiDAR Powerline Point Clouds: Mitigating Adjacent Powerline Interference](https://www.mdpi.com/1424-8220/25/20/6448).

## Project Structure

```
├── fit_line.py         # Powerline fitting with output span
└── PCGrid
    ├── CMakeLists.txt  # for c++ entry
    ├── ext.cpp         # bindings
    ├── pcgrid          # python package __init__
    ├── setup.py        # setup file
    ├── src             # source code of PCGrid
    └── third_party     # third-party dependencies
```

## Installation

### 1. Clone

```bash
git clone https://github.com/costune/Powerline-Span-Segmentation.git --recursive
```

### 2. Using PCGrid in C++

Due to using [libLAS](https://github.com/libLAS/libLAS) to read and write las file, and the repo is no longer maintained. The `libboost` in newer system may not compatible with libLAS. For adaptability, We use `Boost-1.65.0` locally. If you are on Unix system, you can download tar file from [Boost-1.65.0](https://www.boost.org/releases/1.65.0/), and extract the lib folder into `PCGrid/third_party`. 

```bash
cd PCGrid/third_party
curl -O https://archives.boost.io/release/1.65.0/source/boost_1_65_0.tar.gz
tar -xf boost_1_65_0.tar.gz

cd boost_1_65_0
chmod +x ./bootstrap.sh
./b2 -j$(nproc)
```

> If you use other versions of Boost, you may modify the Boost_INCLUDE_DIR and Boost_LIBRARY_DIRS in CMakeLists. If you install the Boost into system, you need to set Boost_INCLUDE_DIR=/usr/local/include and Boost_LIBRARY_DIRS=/usr/local/lib.

We use `CMake` to manage the code.

```bash
mkdir -p PCGrid/build && cd PCGrid/build
cmake ..
make -j$(nproc)
```

The bin file is under `PCGrid/build/bin`.

> If encounted libboost_system.so.1.65.0 not find, add lib path of Boost into env LD_LIBRARY_PATH.

### 3. Using PCGrid in Python

```bash
conda create -n pcgrid python=3.8
conda activate pcgrid
pip install laspy numpy scipy scikit-learn pyshp pyproj tqdm pybind11
pip install PCGrid/
```

## Usage

### 1. C++ bin file

```bash
./PCGrid/build/bin/Las2PowerLine lineLasFile towerLasFile outDir [--verbose]
```

### 2. Python

```Python
from pcgrid import SpanSegment
```

## Examples

We provide a python script for fitting powerline catenary arguments. Using `-h` or `--help` for flag helping.

```bash
python fit_line.py -h
```


## BibTeX

If you find our project help in your research, please cite:

```
@Article{ma2025Enhanced,
    AUTHOR = {Ma, Yi and Wang, Guofang and Liu, Tianle and Wang, Yifan and Geng, Hao and Jiang, Wanshou},
    TITLE = {Enhanced Automatic Span Segmentation of Airborne LiDAR Powerline Point Clouds: Mitigating Adjacent Powerline Interference},
    JOURNAL = {Sensors},
    VOLUME = {25},
    YEAR = {2025},
    NUMBER = {20},
    ARTICLE-NUMBER = {6448},
}
```
