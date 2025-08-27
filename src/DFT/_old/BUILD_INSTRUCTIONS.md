# Building LAMMPS with DFT Package

## Prerequisites

Before building, ensure you have:
- libxc library installed
- libint2 library installed  
- Eigen3 library installed
- nlohmann_json library installed

## Build Instructions

```bash
cd ~/github/lammps-alphataubio
mkdir build-dft
cd build-dft

# Configure with DFT package enabled
cmake ../cmake \
  -DPKG_DFT=yes \
  -DCMAKE_PREFIX_PATH="/path/to/libxc;/path/to/libint2;/path/to/eigen3" \
  -DCMAKE_BUILD_TYPE=Release

# Build
make -j4
```

## Required Libraries Installation (macOS)

```bash
# Using Homebrew
brew install libxc
brew install eigen
brew install nlohmann-json

# libint2 may need to be built from source:
git clone https://github.com/evaleev/libint.git
cd libint
./autogen.sh
./configure --prefix=/usr/local
make -j4
make install
```

## Usage

```lammps
# LAMMPS input script example
units real
atom_style atomic

# Define pair style
# pair_style dft <functional> <basis.json> [options]
pair_style dft B3LYP sto-3g.json grid 50000 tol 1e-8

# Set pair coefficients (atom_type atom_type charge vdw_radius)
pair_coeff 1 1 1.0 1.2  # H atom
pair_coeff 1 2 1.0 1.2  # H-C interaction  
pair_coeff 2 2 6.0 1.7  # C atom

# Run simulation
run 100
```

## Supported Functionals

Any functional from libxc library, including:
- LDA: LDA, PW92, VWN
- GGA: PBE, BLYP, BP86, PW91
- Hybrid: B3LYP, PBE0, HSE06
- Meta-GGA: TPSS, M06, SCAN
- Range-separated: wB97, wB97X, CAM-B3LYP
- Special: wB97M-V (custom implementation)

## Basis Sets

Basis sets should be in JSON format from Basis Set Exchange:
```bash
# Download basis set
curl -X GET "https://www.basissetexchange.org/api/basis/sto-3g/format/json" > sto-3g.json
```

## Troubleshooting

If you encounter linking errors, ensure:
1. All required libraries are in your library path
2. CMake can find all dependencies
3. The DFT source files are properly included in the build

For macOS ARM64 (M1/M2), you may need:
```bash
export CXXFLAGS="-arch arm64"
export LDFLAGS="-L/opt/homebrew/lib"
export CPPFLAGS="-I/opt/homebrew/include"
```
