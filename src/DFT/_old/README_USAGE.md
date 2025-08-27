# DFT Package for LAMMPS

This package implements a full density functional theory (DFT) pair style for LAMMPS, enabling quantum mechanical calculations within molecular dynamics simulations.

## Features

- **Extensive Functional Support**: Any exchange-correlation functional from the libxc library
- **Flexible Basis Sets**: JSON format basis sets from Basis Set Exchange
- **Efficient Integrals**: Powered by libint2 for fast integral evaluation
- **Advanced SCF**: DIIS acceleration for rapid convergence
- **Complete Forces**: Both Hellmann-Feynman and Pulay forces (essential for Gaussian basis sets!)
- **Dispersion Corrections**: D3, D3BJ, and D4 methods
- **Special Functionals**: Custom implementation of wB97M-V

## Quick Start

### 1. Download a basis set:

```bash
python3 download_basis.py def2-svp H
# This creates def2-svp.H.json
```

### 2. Create a LAMMPS input script:

```lammps
# Simple H2 calculation
units real
atom_style atomic

# Create H2 molecule
region box block -10 10 -10 10 -10 10
create_box 1 box
create_atoms 1 single 0.0 0.0 -0.35
create_atoms 1 single 0.0 0.0  0.35
mass 1 1.00794

# Use DFT with PBE functional
pair_style dft PBE def2-svp.H.json
pair_coeff 1 1 1.0 0.5  # H: charge=1, vdw_radius=0.5

# Run
minimize 1e-8 1e-10 1000 10000
```

### 3. Run LAMMPS:

```bash
lmp -in in.h2
```

## Supported Functionals

The package automatically maps common functional names to libxc format:

| Common Name | LibXC Format | Type |
|------------|--------------|------|
| PBE | GGA_X_PBE+GGA_C_PBE | GGA |
| B3LYP | HYB_GGA_XC_B3LYP | Hybrid GGA |
| LDA | LDA_X+LDA_C_PW | LDA |
| BLYP | GGA_X_B88+GGA_C_LYP | GGA |
| PBE0 | HYB_GGA_XC_PBEH | Hybrid GGA |
| TPSS | MGGA_X_TPSS+MGGA_C_TPSS | Meta-GGA |
| SCAN | MGGA_X_SCAN+MGGA_C_SCAN | Meta-GGA |
| M06 | HYB_MGGA_XC_M06 | Hybrid Meta-GGA |
| wB97M-V | (custom) | Range-separated hybrid |

You can also use any libxc functional directly by its full name (e.g., GGA_X_PBE+GGA_C_PBE).

## Basis Sets

Download basis sets from Basis Set Exchange using the provided script:

```bash
# Download common basis sets
python3 download_basis.py sto-3g H C N O
python3 download_basis.py def2-svp H C N O  
python3 download_basis.py 6-31g H C N O
python3 download_basis.py cc-pvdz H C N O
```

## pair_style Syntax

```lammps
pair_style dft <functional> <basis.json> [options]
```

Options:
- `grid <N>` - Number of grid points for XC integration (default: 50000)
- `tol <tol>` - SCF energy convergence tolerance (default: 1e-8)
- `maxiter <N>` - Maximum SCF iterations (default: 100)
- `dispersion <type>` - Add dispersion correction: D3, D3BJ, or D4

## pair_coeff Syntax

```lammps
pair_coeff i j charge vdw_radius [cutoff]
```

- `i j` - Atom types
- `charge` - Nuclear charge (e.g., 1 for H, 6 for C)
- `vdw_radius` - van der Waals radius in Angstroms
- `cutoff` - Optional pair cutoff (default: 20 Angstroms)

## Examples

### Water molecule with PBE:

```lammps
units real
atom_style atomic

# Create water molecule
region box block -10 10 -10 10 -10 10
create_box 2 box
create_atoms 2 single 0.0 0.0 0.0       # O
create_atoms 1 single 0.757 0.586 0.0   # H1
create_atoms 1 single -0.757 0.586 0.0  # H2

mass 1 1.00794   # H
mass 2 15.9994   # O

pair_style dft PBE def2-svp.HO.json grid 75000
pair_coeff 1 1 1.0 0.5   # H-H
pair_coeff 1 2 1.0 0.5   # H-O (mixed)
pair_coeff 2 2 8.0 1.5   # O-O

minimize 1e-9 1e-11 1000 10000
```

### Methane with B3LYP and dispersion:

```lammps
pair_style dft B3LYP 6-31g.CH.json dispersion D3BJ
pair_coeff 1 1 1.0 0.5   # H
pair_coeff 2 2 6.0 1.7   # C
```

## Building with CMake

```bash
cd lammps
mkdir build && cd build
cmake ../cmake -DPKG_DFT=yes \
  -DCMAKE_PREFIX_PATH="/path/to/libxc;/path/to/libint2;/path/to/eigen3"
make -j4
```

## Requirements

- libxc (>= 5.0)
- libint2 (>= 2.6)
- Eigen3 (>= 3.3)
- nlohmann_json (>= 3.0)
- BLAS/LAPACK

## Limitations

- Currently supports closed-shell systems only
- No periodic boundary conditions for DFT (molecules only)
- Grid integration may need tuning for heavy elements

## Troubleshooting

### "Unknown XC functional" error
- Check that the functional name is correct
- Use the mapped common names (PBE, B3LYP, etc.) or full libxc names
- Verify libxc is properly installed

### Basis set loading errors
- Ensure the JSON file exists and is properly formatted
- Download directly from Basis Set Exchange using the provided script
- Check that all elements in your system are included in the basis file

### SCF convergence issues
- Increase grid size for better XC integration
- Adjust mixing parameter in the code
- Try different initial guess strategies

## Citation

If you use this package, please cite:
- LAMMPS: https://lammps.org
- libxc: https://www.tddft.org/programs/libxc/
- libint2: https://github.com/evaleev/libint
- Basis Set Exchange: https://www.basissetexchange.org/

## License

This package is distributed under the same license as LAMMPS (GPL v2).
