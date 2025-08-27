## Compilation Fixes Applied

### Fixed libint2 API Issues

1. **Shell Construction**: Fixed to use proper libint2::Shell constructor
   - Changed from direct member assignment to constructor with Contraction struct
   - Properly set angular momentum and coefficients

2. **Engine Initialization**: Updated to compute max_l and max_nprim correctly
   - Removed dependency on non-existent BasisSet methods
   - Manually compute maximum values from shell data

3. **Include Headers**: Added necessary libint2 headers
   - Added libint2/basis.h, libint2/shell.h, libint2/engine.h
   - Added algorithm header for std::max

4. **Variable Ordering**: Fixed variable declaration order in loops
   - Ensured n1, n2 are declared before use
   - Consistent pattern across all compute functions

### To Compile:

```bash
cd ~/github/lammps-alphataubio/src/DFT
chmod +x test_compile.sh
./test_compile.sh
```

### If you still get errors:

1. Check libint2 version:
```bash
pkg-config --modversion libint2
```

2. Verify library paths:
```bash
pkg-config --cflags libint2
pkg-config --libs libint2
```

3. For minimal test without LAMMPS:
```bash
clang++ -std=c++17 -I/opt/homebrew/include -I/opt/homebrew/include/eigen3 \
        -c integral_engine.cpp -o integral_engine.o
```

### Notes:
- The Shell API in libint2 uses boost::container::small_vector internally
- Contraction structure contains angular momentum and coefficients
- Shell constructor takes (exponents, contractions, origin)
