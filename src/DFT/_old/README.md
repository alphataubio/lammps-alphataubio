# LAMMPS DFT Package

## Overview

The DFT package provides a `pair_style dft` that enables density functional theory (DFT) calculations within LAMMPS molecular dynamics simulations. This implementation interfaces with LibXC for exchange-correlation functionals and uses libint2 for efficient integral evaluation, providing access to a wide range of DFT methods including LDA, GGA, meta-GGA, and hybrid functionals.

## Key Features

- **Extensive Functional Support**: Access to 600+ functionals through LibXC
- **Special Implementation of wB97M-V**: Full support for the highly accurate wB97M-V functional
- **Efficient Integral Evaluation**: Uses libint2 for fast and accurate molecular integrals
- **Hybrid Functional Support**: Handles exact exchange mixing for hybrid functionals
- **Range-Separated Functionals**: Supports long-range corrected functionals
- **Dispersion Corrections**: D3, D3BJ, and D4 dispersion corrections available
- **JSON Basis Sets**: Compatible with Basis Set Exchange JSON format
- **Parallel Execution**: Optimized for MPI and OpenMP parallelization

## Building LAMMPS with DFT Package

### Prerequisites

1. **C++ Compiler** with C++14 support
2. **CMake** >= 3.16
3. **MPI** (optional but recommended)

The following dependencies will be automatically downloaded and built if not found:
- **LibXC** (>= 5.0.0)
- **libint2** (>= 2.6.0)
- **Eigen3** (>= 3.3)

### Build Instructions

```bash
cd lammps-alphataubio
mkdir build
cd build

# Basic build with DFT package
cmake ../cmake -DPKG_DFT=ON

# Build with additional optimizations
cmake ../cmake -DPKG_DFT=ON -DBUILD_OMP=ON -DBUILD_MPI=ON

# Compile
make -j$(nproc)
```

#### Using System Libraries (Optional)

If you have LibXC and libint2 already installed:

```bash
cmake ../cmake -DPKG_DFT=ON \
  -DLIBXC_ROOT=/path/to/libxc \
  -DLibint2_DIR=/path/to/libint2/lib/cmake/libint2
```

## Usage

### Basic Syntax

```lammps
pair_style dft <functional> <basis.json> [options]
```

**Parameters:**
- `functional`: Name of the XC functional (e.g., "PBE", "B3LYP", "wB97M-V")
- `basis.json`: Path to JSON file containing basis set definition
- `options`:
  - `grid <N>`: Number of grid points for numerical integration (default: 50000)
  - `grid_type <type>`: Grid type - "Lebedev", "Becke", "MultiExp" (default: "Lebedev")
  - `tol <value>`: SCF convergence tolerance (default: 1e-8)
  - `maxiter <N>`: Maximum SCF iterations (default: 100)
  - `dispersion <type>`: Dispersion correction - "D3", "D3BJ", "D4" (optional)

### Pair Coefficients

```lammps
pair_coeff <i> <j> <nuclear_charge> <vdw_radius> [cutoff]
```

**Parameters:**
- `i`, `j`: Atom type indices
- `nuclear_charge`: Nuclear charge (number of protons)
- `vdw_radius`: van der Waals radius (Angstroms)
- `cutoff`: Optional pairwise cutoff distance

### Example Input Scripts

#### Water Molecule with wB97M-V

```lammps
units           real
atom_style      full

# Create water molecule
region          box block -10 10 -10 10 -10 10
create_box      2 box

create_atoms    2 single 0.0 0.0 0.0      # Oxygen
create_atoms    1 single 0.757 0.586 0.0  # Hydrogen
create_atoms    1 single -0.757 0.586 0.0 # Hydrogen

mass            1 1.008   # H
mass            2 15.999  # O

# DFT with wB97M-V functional
pair_style      dft wB97M-V def2-tzvp.H_O.json grid 75000 tol 1e-9

pair_coeff      1 1 1.0 1.20   # H-H
pair_coeff      1 2 1.0 1.20   # H-O
pair_coeff      2 2 8.0 1.52   # O-O

minimize        1e-10 1e-12 1000 10000
```

## Supported Functionals

### Local Density Approximation (LDA)
- `LDA`: Slater exchange + VWN correlation
- `LDA_X+LDA_C_PW`: Slater exchange + Perdew-Wang correlation
- `LDA_X+LDA_C_PZ`: Slater exchange + Perdew-Zunger correlation

### Generalized Gradient Approximation (GGA)
- `PBE`: Perdew-Burke-Ernzerhof
- `BLYP`: Becke88 exchange + Lee-Yang-Parr correlation
- `PW91`: Perdew-Wang 1991
- `RPBE`: Revised PBE
- `revPBE`: Revised PBE by Zhang-Yang
- `OPTPBE`: Optimized PBE
- `OPTB88`: Optimized B88

### Meta-GGA Functionals
- `TPSS`: Tao-Perdew-Staroverov-Scuseria
- `revTPSS`: Revised TPSS
- `SCAN`: Strongly Constrained and Appropriately Normed
- `r2SCAN`: Regularized SCAN
- `M06L`: Minnesota 06 local functional

### Hybrid Functionals
- `B3LYP`: Becke 3-parameter Lee-Yang-Parr (20% HF exchange)
- `PBE0`: PBE hybrid (25% HF exchange)
- `HSE06`: Heyd-Scuseria-Ernzerhof (25% short-range HF exchange)
- `M062X`: Minnesota 06 2X (54% HF exchange)

### Range-Separated Functionals
- `wB97M-V`: ωB97M-V with VV10 non-local correlation
- `wB97X-V`: ωB97X-V with VV10 non-local correlation
- `CAM-B3LYP`: Coulomb-attenuating method B3LYP
- `LC-PBE`: Long-range corrected PBE
- `LC-wPBE`: Long-range corrected ωPBE

## Basis Set Format

The package accepts basis sets in JSON format compatible with the Basis Set Exchange (BSE). Example structure:

```json
{
  "molssi_bse_schema": {
    "schema_type": "complete",
    "schema_version": "0.1"
  },
  "name": "def2-TZVP",
  "elements": {
    "1": {
      "electron_shells": [
        {
          "function_type": "gto",
          "angular_momentum": [0],
          "exponents": ["34.0613410", "5.1235746", "1.1646626"],
          "coefficients": [["0.60251978E-02", "0.45021094E-01", "0.20189726"]]
        }
      ]
    }
  }
}
```

You can download basis sets from: https://www.basissetexchange.org/

## Performance Considerations

### Grid Size
- **Coarse** (25000 points): Quick calculations, qualitative results
- **Medium** (50000 points): Good balance of speed and accuracy
- **Fine** (75000 points): Production calculations
- **Ultrafine** (100000+ points): High accuracy, benchmark quality

### Functional Choice Impact
- **LDA**: Fastest, least accurate
- **GGA**: 2-3x slower than LDA, good accuracy
- **Meta-GGA**: 3-5x slower than LDA, better accuracy
- **Hybrid**: 10-50x slower than LDA, high accuracy
- **Range-separated**: Similar to hybrid

### Parallelization
- **MPI**: Distributes atoms and grid points across processors
- **OpenMP**: Parallelizes integral evaluation and grid operations
- **GPU**: Not yet supported (planned for future release)

### Memory Requirements
- Scales as O(N²) with number of basis functions
- Typical requirements:
  - Small molecule (< 50 atoms): 1-2 GB
  - Medium system (50-200 atoms): 4-8 GB
  - Large system (> 200 atoms): 16+ GB

## Advanced Features

### Self-Consistent Field (SCF) Control

```lammps
# Tight convergence for accurate forces
pair_style dft PBE basis.json tol 1e-10 maxiter 200

# Loose convergence for rapid screening
pair_style dft PBE basis.json tol 1e-6 maxiter 50
```

### Dispersion Corrections

```lammps
# D3 dispersion with Becke-Johnson damping
pair_style dft PBE basis.json dispersion D3BJ

# D4 dispersion (most accurate, slower)
pair_style dft PBE basis.json dispersion D4
```

### Grid Types

```lammps
# Lebedev grid (spherical, recommended)
pair_style dft PBE basis.json grid_type Lebedev

# Becke grid (atomic partitioning)
pair_style dft PBE basis.json grid_type Becke

# MultiExp grid (multi-exponential, experimental)
pair_style dft PBE basis.json grid_type MultiExp
```

## Limitations

1. **Periodic Boundary Conditions**: Currently limited support for PBC
2. **Open-Shell Systems**: Only closed-shell (spin-restricted) calculations
3. **Relativistic Effects**: No relativistic corrections available
4. **Solvation**: No implicit solvation models
5. **Excited States**: Ground state only

## Troubleshooting

### SCF Convergence Issues
- Increase max iterations: `maxiter 200`
- Use density mixing: Automatically applied
- Try different initial guess
- Check molecular geometry for unrealistic bonds

### Memory Errors
- Reduce grid size
- Use smaller basis set
- Enable memory-efficient algorithms
- Increase system RAM or use fewer MPI ranks

### Incorrect Energies
- Verify basis set file format
- Check nuclear charges in pair_coeff
- Ensure appropriate functional for system
- Increase grid density for meta-GGA functionals

## Citations

When using the DFT package, please cite:

1. **LAMMPS**: Thompson, A. P., et al. "LAMMPS - a flexible simulation tool for particle-based materials modeling at the atomic, meso, and continuum scales." Comp. Phys. Comm. 271, 108171 (2022).

2. **LibXC**: Lehtola, S., et al. "Recent developments in libxc — A comprehensive library of functionals for density functional theory." SoftwareX 7, 1-5 (2018).

3. **libint2**: Valeev, E. F. "Libint: A library for the evaluation of molecular integrals of many-body operators over Gaussian functions." http://libint.valeyev.net/ (2020).

4. **wB97M-V** (if used): Mardirossian, N. & Head-Gordon, M. "ωB97M-V: A combinatorially optimized, range-separated hybrid, meta-GGA density functional with VV10 nonlocal correlation." J. Chem. Phys. 144, 214110 (2016).

## Contributing

Contributions are welcome! Areas for improvement:
- GPU acceleration (CUDA/HIP)
- Open-shell (unrestricted) calculations
- Periodic boundary conditions with k-points
- Additional basis set formats
- Implicit solvation models
- Excited state methods (TDDFT)

Submit pull requests or issues to: https://github.com/lammps/lammps

## License

The DFT package is distributed under the same license as LAMMPS (GPL v2).

## Contact

For questions and support:
- LAMMPS mailing list: https://www.lammps.org/mail.html
- GitHub issues: https://github.com/lammps/lammps/issues
