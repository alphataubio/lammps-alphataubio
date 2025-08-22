# LAMMPS DFT Package

## Overview

The DFT package provides a `pair_style libxc` that interfaces with the LibXC library to enable density functional theory (DFT) calculations within LAMMPS molecular dynamics simulations. This allows for the use of any exchange-correlation functional available in LibXC, including LDA, GGA, meta-GGA, and hybrid functionals.

## Features

- **Wide Range of Functionals**: Access to over 600 exchange-correlation functionals through LibXC
- **Flexible Functional Specification**: Use separate exchange and correlation functionals or combined XC functionals
- **Efficient Implementation**: Optimized for parallel execution in LAMMPS
- **Customizable Parameters**: Control grid spacing and integration parameters

## Building with the DFT Package

### Prerequisites

1. **LibXC Library** (version >= 5.0.0)
   - Install from package manager: `sudo apt-get install libxc-dev` (Ubuntu/Debian)
   - Or build from source: https://www.tddft.org/programs/libxc/

2. **LAMMPS Build Requirements**
   - CMake >= 3.16
   - C++ compiler with C++11 support
   - MPI (optional, for parallel builds)

### Building LAMMPS with DFT Package

```bash
cd ~/github/lammps-alphataubio
mkdir build
cd build
cmake ../cmake -DPKG_DFT=ON
make -j$(nproc)
```

If LibXC is installed in a non-standard location:
```bash
cmake ../cmake -DPKG_DFT=ON -DLIBXC_ROOT=/path/to/libxc
```

## Usage

### Pair Style Syntax

```lammps
pair_style libxc <functional> <cutoff> [options]
```

**Parameters:**
- `functional`: LibXC functional name or combination (e.g., "PBE", "LDA_X+LDA_C_PW")
- `cutoff`: Global cutoff distance
- `options`:
  - `ngrid <N>`: Number of grid points for integration (default: 100)
  - `spacing <S>`: Grid spacing for integration (default: 0.01)

### Pair Coefficients

```lammps
pair_coeff <i> <j> <rho0> <decay_length> <atomic_volume> <cutoff>
```

**Parameters:**
- `i`, `j`: Atom type indices
- `rho0`: Reference electron density
- `decay_length`: Characteristic decay length for electron density
- `atomic_volume`: Effective atomic volume
- `cutoff`: Pairwise cutoff distance

### Examples

#### LDA Functional
```lammps
pair_style libxc LDA_X+LDA_C_VWN 12.0
pair_coeff 1 1 0.15 2.5 11.8 12.0
```

#### PBE GGA Functional
```lammps
pair_style libxc PBE 12.0 ngrid 150
pair_coeff 1 1 0.15 2.5 11.8 12.0
```

#### B3LYP Hybrid Functional
```lammps
pair_style libxc B3LYP 15.0 ngrid 200 spacing 0.005
pair_coeff 1 1 0.15 2.5 11.8 15.0
```

## Supported Functionals

### Common LDA Functionals
- `LDA_X`: Slater exchange
- `LDA_C_VWN`: Vosko-Wilk-Nusair correlation
- `LDA_C_PW`: Perdew-Wang correlation
- `LDA_C_PZ`: Perdew-Zunger correlation

### Common GGA Functionals
- `PBE`: Perdew-Burke-Ernzerhof
- `BLYP`: Becke88 exchange + Lee-Yang-Parr correlation
- `GGA_X_PBE+GGA_C_PBE`: PBE exchange and correlation (separate)
- `PW91`: Perdew-Wang 1991

### Hybrid Functionals
- `B3LYP`: Becke 3-parameter Lee-Yang-Parr
- `PBE0`: PBE hybrid
- `HSE06`: Heyd-Scuseria-Ernzerhof

### Meta-GGA Functionals
- `TPSS`: Tao-Perdew-Staroverov-Scuseria
- `SCAN`: Strongly Constrained and Appropriately Normed

For a complete list, see: https://www.tddft.org/programs/libxc/functionals/

## Parameter Fitting

The pair coefficients (`rho0`, `decay_length`, `atomic_volume`) should be fitted to reproduce:
1. DFT-calculated cohesive energies
2. Lattice parameters
3. Elastic constants
4. Surface energies

A fitting procedure can be implemented using:
1. Reference DFT calculations (VASP, Quantum ESPRESSO, etc.)
2. Force matching or energy matching methods
3. Optimization algorithms (genetic algorithms, simulated annealing)

## Performance Considerations

1. **Functional Choice**: 
   - LDA functionals are fastest
   - GGA functionals are ~2x slower than LDA
   - Hybrid functionals are significantly slower (10-100x)

2. **Grid Parameters**:
   - Increasing `ngrid` improves accuracy but increases cost
   - Grid spacing affects convergence

3. **Parallelization**:
   - The implementation is MPI-parallel
   - Scales well up to hundreds of processors

## Limitations

1. Currently implements a simplified electron density model
2. Self-consistent field iterations not included (uses approximate densities)
3. Spin polarization not yet supported
4. van der Waals functionals require additional implementation

## Testing

Example input scripts are provided in `examples/DFT/`:
- `in.libxc`: Basic LDA calculation
- `in.libxc_advanced`: GGA calculation with binary alloy

Run tests:
```bash
cd examples/DFT
../../build/lmp -in in.libxc
```

## References

1. LibXC: Lehtola, S., et al. "Recent developments in libxc — A comprehensive library of functionals for density functional theory." SoftwareX 7, 1-5 (2018).

2. DFT in MD: Marx, D., and Hutter, J. "Ab initio molecular dynamics: basic theory and advanced methods." Cambridge University Press (2009).

3. LAMMPS: Thompson, A. P., et al. "LAMMPS - a flexible simulation tool for particle-based materials modeling at the atomic, meso, and continuum scales." Comp. Phys. Comm. 271, 108171 (2022).

## Contributing

Contributions are welcome! Please submit pull requests or issues to:
- Add new density models
- Implement spin polarization
- Add van der Waals corrections
- Improve performance

## License

This package is distributed under the same license as LAMMPS (GPL v2).
