.. index:: pair_style dft

pair_style dft command
======================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style dft functional basis_file keyword value ...

* functional = name of exchange-correlation functional

  * Common names: *PBE*, *B3LYP*, *LDA*, *BLYP*, *BP86*, *PBE0*, *HSE06*, *TPSS*, *SCAN*, *M06*, *wB97M-V*
  * LibXC format: *GGA_X_PBE+GGA_C_PBE*, *HYB_GGA_XC_B3LYP*, etc.
  
* basis_file = JSON file containing basis set from Basis Set Exchange
* zero or more keyword/value pairs may be appended
* keyword = *grid* or *tol* or *maxiter* or *dispersion* or *grid_type*

  .. parsed-literal::

     *grid* value = N
       N = number of grid points for XC integration (default: 50000)
     *tol* value = tolerance
       tolerance = SCF energy convergence tolerance (default: 1e-8)
     *maxiter* value = N
       N = maximum number of SCF iterations (default: 100)
     *dispersion* value = type
       type = *D3* or *D3BJ* or *D4* = dispersion correction type
     *grid_type* value = type
       type = *Lebedev* or *Becke* or *MultiExp* = integration grid type (default: Lebedev)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style dft PBE def2-svp.json
   pair_coeff 1 1 1.0 0.5
   
   pair_style dft B3LYP 6-31g.json grid 75000 tol 1e-9
   pair_coeff * * 6.0 1.7
   
   pair_style dft GGA_X_PBE+GGA_C_PBE cc-pvdz.json dispersion D3BJ
   pair_coeff 1 1 1.0 0.5   # H
   pair_coeff 2 2 6.0 1.7   # C
   
   pair_style dft SCAN def2-tzvp.json grid 100000 maxiter 200
   pair_coeff 1 2 1.0 0.5   # H-C interaction
   
   pair_style dft wB97M-V def2-qzvpp.json grid_type Becke
   pair_coeff * * 8.0 1.5   # O

Description
"""""""""""

Style *dft* performs self-consistent field (SCF) density functional theory 
calculations using atom-centered Gaussian basis sets and exchange-correlation 
(XC) functionals. This allows for quantum mechanical treatment of electronic 
structure within LAMMPS molecular dynamics simulations.

The total DFT energy is computed as:

.. math::

   E_{DFT} = T_s[\rho] + V_{ne}[\rho] + J[\rho] + E_{xc}[\rho] + E_{nn}

where :math:`T_s` is the kinetic energy of non-interacting electrons,
:math:`V_{ne}` is the electron-nuclear attraction, :math:`J` is the classical
electron-electron repulsion (Coulomb), :math:`E_{xc}` is the exchange-correlation
energy, and :math:`E_{nn}` is the nuclear-nuclear repulsion.

The electronic structure is solved using the Roothaan-Hall equations:

.. math::

   \mathbf{F}\mathbf{C} = \mathbf{S}\mathbf{C}\boldsymbol{\epsilon}

where **F** is the Fock matrix, **C** contains molecular orbital coefficients,
**S** is the overlap matrix, and :math:`\boldsymbol{\epsilon}` contains orbital energies.

Forces on atoms include both Hellmann-Feynman and Pulay contributions:

.. math::

   \mathbf{F}_{total} = \mathbf{F}_{HF} + \mathbf{F}_{Pulay}

The Pulay force term is essential for atom-centered Gaussian basis sets because
the basis functions move with the atoms.

----------

**Supported Functionals**

The *dft* style supports all exchange-correlation functionals from the 
`LibXC library <https://www.tddft.org/programs/libxc/>`_. Common functional 
names are automatically mapped to their LibXC equivalents:

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Common Name
     - LibXC Format
     - Type
   * - PBE
     - GGA_X_PBE+GGA_C_PBE
     - GGA
   * - B3LYP
     - HYB_GGA_XC_B3LYP
     - Hybrid GGA (20% HF exchange)
   * - LDA
     - LDA_X+LDA_C_PW
     - LDA
   * - BLYP
     - GGA_X_B88+GGA_C_LYP
     - GGA
   * - BP86
     - GGA_X_B88+GGA_C_P86
     - GGA
   * - PBE0
     - HYB_GGA_XC_PBEH
     - Hybrid GGA (25% HF exchange)
   * - HSE06
     - HYB_GGA_XC_HSE06
     - Range-separated hybrid
   * - TPSS
     - MGGA_X_TPSS+MGGA_C_TPSS
     - Meta-GGA
   * - SCAN
     - MGGA_X_SCAN+MGGA_C_SCAN
     - Meta-GGA
   * - M06
     - HYB_MGGA_XC_M06
     - Hybrid meta-GGA (27% HF)
   * - M06-2X
     - HYB_MGGA_XC_M06_2X
     - Hybrid meta-GGA (54% HF)
   * - wB97M-V
     - (custom implementation)
     - Range-separated meta-GGA with VV10 NLC

You can also specify any LibXC functional directly using its full name.
Separate exchange and correlation functionals can be combined using the 
"+" notation (e.g., GGA_X_PBE+GGA_C_PBE).

----------

**Basis Sets**

Basis sets must be provided in JSON format as obtained from the
`Basis Set Exchange <https://www.basissetexchange.org/>`_. A helper
script is provided to download basis sets:

.. code-block:: bash

   python download_basis.py def2-svp H C N O
   # Creates def2-svp.HCNO.json

Common basis sets include:

* Minimal: STO-3G, STO-6G
* Pople: 3-21G, 6-31G, 6-31G*, 6-311G**
* Dunning: cc-pVDZ, cc-pVTZ, cc-pVQZ, aug-cc-pVDZ
* Karlsruhe: def2-SVP, def2-TZVP, def2-QZVP
* Ahlrichs: def2-TZVPP, def2-QZVPP

The basis set file must contain basis functions for all element types
present in your simulation.

----------

**Grid Integration**

Exchange-correlation energies for GGA and meta-GGA functionals are 
evaluated numerically on a grid. The package uses:

* **Radial grid**: Chebyshev-Gauss quadrature with Becke transformation
* **Angular grid**: Lebedev quadrature on the unit sphere
* **Partitioning**: Becke partitioning for multi-center integration

The number of grid points can be adjusted with the *grid* keyword:

* 30000-50000: Coarse grid, suitable for geometry optimization
* 50000-75000: Medium grid, default for most calculations
* 75000-100000: Fine grid, recommended for accurate energies
* 100000+: Ultra-fine grid, for benchmark calculations

----------

**SCF Convergence**

The self-consistent field calculation uses:

* **DIIS acceleration**: Direct Inversion of Iterative Subspace for faster convergence
* **Density mixing**: Linear mixing with adjustable parameter
* **Initial guess**: Core Hamiltonian (H_core = T + V_ne)

Convergence is achieved when both energy and density changes fall below
the specified tolerance.

----------

**Dispersion Corrections**

Long-range van der Waals interactions can be included via Grimme's
DFT-D methods:

* **D3**: Original D3 dispersion with zero-damping
* **D3BJ**: D3 with Becke-Johnson damping (recommended)
* **D4**: Latest version with charge-dependent C6 coefficients

These corrections are essential for accurate description of non-covalent
interactions.

----------

**Implementation Details**

The *dft* style uses:

* `LibXC <https://www.tddft.org/programs/libxc/>`_ for exchange-correlation functionals
* `Libint2 <https://github.com/evaleev/libint>`_ for molecular integrals
* `Eigen3 <http://eigen.tuxfamily.org/>`_ for linear algebra
* Custom implementation of wB97M-V functional with VV10 non-local correlation

The implementation includes:

#. **One-electron integrals**: Overlap (S), kinetic (T), nuclear attraction (V)
#. **Two-electron integrals**: Electron repulsion integrals (ERIs) with screening
#. **Fock matrix construction**: F = H_core + 2J - K + V_xc
#. **SCF solver**: Roothaan-Hall equations with DIIS
#. **Force calculation**: Hellmann-Feynman + Pulay forces
#. **Parallel efficiency**: MPI parallelization over atoms

----------

**Performance Considerations**

DFT calculations are computationally intensive. For better performance:

* Use smaller basis sets for initial equilibration (e.g., STO-3G or 3-21G)
* Reduce grid size for geometry optimization (30000-50000 points)
* Increase grid size for production runs and property calculations
* Consider using pure functionals (PBE, BLYP) over hybrids (B3LYP, PBE0)
* Meta-GGA functionals (SCAN, TPSS) are more expensive than GGA

Computational scaling:

* Basis functions: O(N^4) for exact exchange, O(N^3) for pure DFT
* Grid points: O(N_grid × N_basis^2)
* SCF iterations: Typically 10-30 for well-behaved systems

----------

The pair_coeff command for style *dft* specifies the nuclear charge
and van der Waals radius for each atom type:

.. code-block:: LAMMPS

   pair_coeff I J charge radius [cutoff]

* I,J = atom types
* charge = nuclear charge (1 for H, 6 for C, 7 for N, 8 for O, etc.)
* radius = van der Waals radius in Angstroms (used for dispersion)
* cutoff = optional pair cutoff (default: 20.0 Angstroms)

Examples:

.. code-block:: LAMMPS

   pair_coeff 1 1 1.0 0.5    # H: Z=1, rvdw=0.5 Å
   pair_coeff 2 2 6.0 1.7    # C: Z=6, rvdw=1.7 Å
   pair_coeff 3 3 7.0 1.55   # N: Z=7, rvdw=1.55 Å
   pair_coeff 4 4 8.0 1.52   # O: Z=8, rvdw=1.52 Å

For cross-interactions (I ≠ J), geometric mean mixing rules are applied
automatically.

----------

**Output Information**

During SCF iterations, the following information is printed:

.. code-block:: none

   ================================================
              DFT SCF CALCULATION
   ================================================
   Functional: PBE
   Basis functions: 28
   Electrons: 10
   Grid points: 50000
   ------------------------------------------------
    Iter    Energy         Delta E      RMS Dens
   ------------------------------------------------
      1 -76.234567890  0.00000e+00  1.23456e-01
      2 -76.345678901  1.11111e-01  2.34567e-02
      3 -76.356789012  1.11101e-02  3.45678e-03
   ------------------------------------------------
   SCF CONVERGED
   
   Energy Components:
     Kinetic:          75.123456789
     Nuclear:         -198.234567890
     Electron-Nuclear: -123.345678901
     Coulomb:          45.678901234
     XC:              -12.345678901
     Dispersion:       -0.123456789
     TOTAL:           -76.356789012
   ================================================

----------

**Special Features**

#. **wB97M-V functional**: Custom implementation of the ωB97M-V 
   range-separated hybrid meta-GGA functional with VV10 non-local 
   correlation, one of the most accurate functionals for general chemistry.

#. **Automatic functional mapping**: Common functional names are 
   automatically converted to LibXC format.

#. **Flexible basis set support**: Any basis set from Basis Set Exchange 
   can be used directly.

#. **Complete force implementation**: Both Hellmann-Feynman and Pulay 
   forces are computed for accurate dynamics with Gaussian basis sets.

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The *dft* style does not support the :doc:`pair_modify <pair_modify>`
mix, shift, table, and tail options.

The *dft* style does not write basis set or functional information
to :doc:`binary restart files <restart>`. You must re-specify the
pair_style and pair_coeff commands in an input script that reads a
restart file.

The *dft* style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command. It does not support the
*inner*, *middle*, *outer* keywords.

----------

Restrictions
""""""""""""

* The DFT package must be enabled when building LAMMPS
* Requires LibXC library (version 5.0 or higher)
* Requires Libint2 library (version 2.6 or higher)
* Requires Eigen3 library (version 3.3 or higher)
* Requires nlohmann_json library
* Currently supports only closed-shell systems (even number of electrons)
* Does not support periodic boundary conditions (molecular systems only)
* Cannot be used with atom_style that doesn't support charges

----------

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`pair_style <pair_style>`,
:doc:`units <units>`

----------

Default
"""""""

* grid = 50000
* tol = 1e-8
* maxiter = 100
* grid_type = Lebedev
* dispersion = none

----------

References
""""""""""

.. [Kohn1965] W. Kohn and L. J. Sham, 
   "Self-Consistent Equations Including Exchange and Correlation Effects",
   Phys. Rev. 140, A1133 (1965).

.. [Mardirossian2016] N. Mardirossian and M. Head-Gordon,
   "ωB97M-V: A combinatorially optimized, range-separated hybrid, meta-GGA 
   density functional with VV10 nonlocal correlation",
   J. Chem. Phys. 144, 214110 (2016).

.. [Grimme2010] S. Grimme, J. Antony, S. Ehrlich, and H. Krieg,
   "A consistent and accurate ab initio parametrization of density functional 
   dispersion correction (DFT-D) for the 94 elements H-Pu",
   J. Chem. Phys. 132, 154104 (2010).

.. [LibXC] M. A. L. Marques, M. J. T. Oliveira, and T. Burnus,
   "LibXC: A library of exchange and correlation functionals for density functional theory",
   Comput. Phys. Commun. 183, 2272 (2012).
   https://www.tddft.org/programs/libxc/

.. [Libint2] E. F. Valeev,
   "Libint: A library for the evaluation of molecular integrals of many-body operators over Gaussian functions",
   https://github.com/evaleev/libint

.. [BSE] B. P. Pritchard, D. Altarawy, B. Didier, T. D. Gibson, and T. L. Windus,
   "A New Basis Set Exchange: An Open, Up-to-date Resource for the Molecular Sciences Community",
   J. Chem. Inf. Model. 59, 4814 (2019).
   https://www.basissetexchange.org/

.. [Pulay1969] P. Pulay,
   "Ab initio calculation of force constants and equilibrium geometries in polyatomic molecules",
   Mol. Phys. 17, 197 (1969).

.. [DIIS] P. Pulay,
   "Convergence acceleration of iterative sequences. The case of SCF iteration",
   Chem. Phys. Lett. 73, 393 (1980).

.. [Becke1988] A. D. Becke,
   "A multicenter numerical integration scheme for polyatomic molecules",
   J. Chem. Phys. 88, 2547 (1988).

.. [Lebedev1999] V. I. Lebedev and D. N. Laikov,
   "A quadrature formula for the sphere of the 131st algebraic order of accuracy",
   Doklady Mathematics 59, 477 (1999).
