# Fix for "Unknown functional: PBE" Error

## Problem
The initial implementation tried to map functional names like "PBE" to LibXC string format (e.g., "GGA_X_PBE+GGA_C_PBE") and then parse those strings with `xc_functional_get_number()`. This approach failed because LibXC doesn't recognize these string formats directly.

## Solution
Implemented direct mapping of common functional names to LibXC integer IDs:

### Supported Common Functionals

| Name | Type | LibXC IDs | Description |
|------|------|-----------|-------------|
| **PBE** | GGA | X:101, C:130 | Perdew-Burke-Ernzerhof |
| **LDA** | LDA | X:1, C:12 | Local Density Approximation (Slater + PW) |
| **B3LYP** | Hybrid | XC:402 | Becke 3-parameter Lee-Yang-Parr (20% HF) |
| **PBE0** | Hybrid | XC:406 | PBE hybrid (25% HF) |
| **BLYP** | GGA | X:106, C:131 | Becke 88 + Lee-Yang-Parr |
| **BP86** | GGA | X:106, C:132 | Becke 88 + Perdew 86 |
| **TPSS** | Meta-GGA | X:202, C:231 | Tao-Perdew-Staroverov-Scuseria |
| **SCAN** | Meta-GGA | X:263, C:267 | Strongly Constrained and Appropriately Normed |
| **wB97M-V** | Special | Custom | Range-separated meta-GGA with VV10 NLC |

## Implementation Details

1. **Direct ID Mapping**: Each common functional name is mapped directly to its LibXC integer ID(s)
2. **Separate X/C Handling**: For functionals like PBE, exchange and correlation are initialized separately
3. **Combined XC Handling**: For functionals like B3LYP, the combined XC functional is used
4. **Fallback**: If not a common name, tries to parse as LibXC format string

## Usage Examples

```lammps
# Common name usage (recommended)
pair_style dft PBE def2-svp.json
pair_style dft B3LYP 6-31g.json
pair_style dft SCAN cc-pvtz.json

# LibXC format (for less common functionals)
pair_style dft HYB_GGA_XC_B3LYP basis.json
pair_style dft GGA_X_PBE+GGA_C_PBE basis.json
```

## Testing

To test if a functional works:

```lammps
# Minimal test
units real
atom_style atomic
region box block -10 10 -10 10 -10 10
create_box 1 box
create_atoms 1 single 0 0 0
mass 1 1.00794

# Test functional initialization
pair_style dft PBE sto-3g.H.json
pair_coeff 1 1 1.0 0.5
run 0
```

## Error Messages

If an unknown functional is used, the error message now shows:
```
Unknown functional: XYZ
Supported common names: PBE, LDA, B3LYP, BLYP, BP86, PBE0, TPSS, SCAN, wB97M-V
Or use LibXC format: HYB_GGA_XC_B3LYP, GGA_X_PBE+GGA_C_PBE, etc.
```

## Files Modified

- `src/DFT/pair_dft.cpp`: Added direct functional mapping and LibXC constants

## Build Requirements

- LibXC library (version 5.0+)
- C++ compiler with C++11 support
- LAMMPS with DFT package enabled

## Notes

- The fix uses LibXC integer constants which are stable across versions
- Meta-GGA functionals (TPSS, SCAN) require kinetic energy density evaluation
- Hybrid functionals (B3LYP, PBE0) include exact exchange and are slower
- wB97M-V has a custom implementation due to its complexity
