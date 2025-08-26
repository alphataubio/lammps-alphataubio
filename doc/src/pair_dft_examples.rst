.. _pair_dft_examples:

pair_style dft examples
========================

This page provides detailed examples of using the *pair_style dft* command
for various quantum mechanical simulations in LAMMPS.

----------

Simple H2 molecule
------------------

Hydrogen molecule with PBE functional:

.. code-block:: LAMMPS

   # H2 molecule with PBE/def2-svp
   units real
   atom_style atomic
   
   # Create box and atoms
   region box block -10 10 -10 10 -10 10
   create_box 1 box
   create_atoms 1 single 0.0 0.0 -0.37
   create_atoms 1 single 0.0 0.0  0.37
   mass 1 1.00794
   
   # DFT settings
   pair_style dft PBE def2-svp.H.json grid 50000 tol 1e-9
   pair_coeff 1 1 1.0 0.5  # H: charge=1, vdw_radius=0.5
   
   # Optimize geometry
   minimize 1e-8 1e-10 1000 10000
   
   # Output
   write_data h2_optimized.data

----------

Water molecule
--------------

Water with hybrid B3LYP functional:

.. code-block:: LAMMPS

   # H2O with B3LYP/6-31G*
   units real
   atom_style atomic
   
   # Create water molecule
   region box block -10 10 -10 10 -10 10
   create_box 2 box
   
   # O at origin, H atoms at typical positions
   create_atoms 2 single 0.0 0.0 0.0
   create_atoms 1 single 0.757 0.586 0.0
   create_atoms 1 single -0.757 0.586 0.0
   
   mass 1 1.00794   # H
   mass 2 15.9994   # O
   
   # DFT with hybrid functional
   pair_style dft B3LYP 6-31gs.HO.json grid 75000 tol 1e-10
   pair_coeff 1 1 1.0 0.5   # H
   pair_coeff 1 2 1.0 0.5   # H-O mixed
   pair_coeff 2 2 8.0 1.52  # O
   
   # Optimize
   minimize 1e-9 1e-11 2000 20000
   
   # Run short MD at 300K
   velocity all create 300.0 12345
   fix 1 all nvt temp 300.0 300.0 100.0
   timestep 0.5
   run 1000

----------

Methane with dispersion
------------------------

CH4 with PBE-D3BJ:

.. code-block:: LAMMPS

   # Methane with dispersion correction
   units real
   atom_style atomic
   
   # Create CH4 (tetrahedral geometry)
   region box block -10 10 -10 10 -10 10
   create_box 2 box
   
   # Carbon at center
   create_atoms 2 single 0.0 0.0 0.0
   
   # Hydrogen atoms at tetrahedral positions
   create_atoms 1 single  0.629  0.629  0.629
   create_atoms 1 single -0.629 -0.629  0.629
   create_atoms 1 single -0.629  0.629 -0.629
   create_atoms 1 single  0.629 -0.629 -0.629
   
   mass 1 1.00794   # H
   mass 2 12.0107   # C
   
   # PBE with D3BJ dispersion
   pair_style dft PBE def2-svp.CH.json dispersion D3BJ
   pair_coeff 1 1 1.0 0.5   # H
   pair_coeff 1 2 1.0 0.5   # H-C
   pair_coeff 2 2 6.0 1.7   # C
   
   # Optimize
   minimize 1e-9 1e-11 1000 10000

----------

Benzene with meta-GGA
----------------------

Benzene ring with SCAN functional:

.. code-block:: LAMMPS

   # Benzene with SCAN meta-GGA
   units real
   atom_style atomic
   boundary p p p
   
   # Read benzene coordinates
   read_data benzene.data
   
   # Or create manually (example positions)
   # region box block -15 15 -15 15 -15 15
   # create_box 2 box
   # create_atoms 2 single 1.396 0.000 0.0  # C atoms in hexagon
   # create_atoms 2 single 0.698 1.209 0.0
   # create_atoms 2 single -0.698 1.209 0.0
   # create_atoms 2 single -1.396 0.000 0.0
   # create_atoms 2 single -0.698 -1.209 0.0
   # create_atoms 2 single 0.698 -1.209 0.0
   # # Add H atoms...
   
   mass 1 1.00794   # H
   mass 2 12.0107   # C
   
   # Meta-GGA functional with fine grid
   pair_style dft SCAN def2-tzvp.CH.json grid 100000 maxiter 200
   pair_coeff 1 1 1.0 0.5   # H
   pair_coeff 2 2 6.0 1.7   # C
   
   # Optimize with tight convergence
   min_style cg
   minimize 1e-10 1e-12 5000 50000

----------

Ammonia with range-separated functional
----------------------------------------

NH3 with wB97M-V:

.. code-block:: LAMMPS

   # Ammonia with wB97M-V range-separated functional
   units real
   atom_style atomic
   
   # Create NH3 (pyramidal geometry)
   region box block -10 10 -10 10 -10 10
   create_box 2 box
   
   # N at origin
   create_atoms 2 single 0.0 0.0 0.0
   
   # H atoms in pyramidal arrangement
   create_atoms 1 single 0.94 0.0 -0.37
   create_atoms 1 single -0.47 0.82 -0.37
   create_atoms 1 single -0.47 -0.82 -0.37
   
   mass 1 1.00794   # H
   mass 2 14.0067   # N
   
   # wB97M-V with ultra-fine grid
   pair_style dft wB97M-V def2-qzvpp.HN.json grid 150000 tol 1e-11
   pair_coeff 1 1 1.0 0.5   # H
   pair_coeff 2 2 7.0 1.55  # N
   
   # Very tight optimization
   minimize 1e-11 1e-13 10000 100000
   
   # Compute vibrational frequencies
   dynamical_matrix all eskm 0.01 file dynmat.dat

----------

CO2 molecule
------------

Carbon dioxide with PBE0 hybrid:

.. code-block:: LAMMPS

   # CO2 linear molecule
   units real
   atom_style atomic
   
   # Create linear CO2
   region box block -10 10 -10 10 -10 10
   create_box 2 box
   
   create_atoms 1 single 0.0 0.0 0.0      # C
   create_atoms 2 single 1.16 0.0 0.0     # O
   create_atoms 2 single -1.16 0.0 0.0    # O
   
   mass 1 12.0107   # C
   mass 2 15.9994   # O
   
   # Hybrid functional PBE0
   pair_style dft PBE0 cc-pvdz.CO.json grid 75000
   pair_coeff 1 1 6.0 1.7   # C
   pair_coeff 2 2 8.0 1.52  # O
   
   # Optimize
   minimize 1e-9 1e-11 1000 10000

----------

Ethanol molecule
----------------

Ethanol with basis set comparison:

.. code-block:: LAMMPS

   # Ethanol - comparing basis sets
   units real
   atom_style atomic
   
   # Read ethanol structure
   read_data ethanol.data
   
   mass 1 1.00794   # H
   mass 2 12.0107   # C
   mass 3 15.9994   # O
   
   # Small basis for quick optimization
   pair_style dft PBE sto-3g.CHO.json grid 30000 tol 1e-7
   pair_coeff 1 1 1.0 0.5   # H
   pair_coeff 2 2 6.0 1.7   # C
   pair_coeff 3 3 8.0 1.52  # O
   
   minimize 1e-7 1e-9 100 1000
   write_data ethanol_sto3g.data
   
   # Larger basis for accurate energy
   pair_style dft PBE def2-tzvp.CHO.json grid 100000 tol 1e-10
   pair_coeff 1 1 1.0 0.5   # H
   pair_coeff 2 2 6.0 1.7   # C
   pair_coeff 3 3 8.0 1.52  # O
   
   minimize 1e-10 1e-12 1000 10000
   write_data ethanol_def2tzvp.data

----------

Hydrogen bonding in water dimer
--------------------------------

Water dimer with dispersion:

.. code-block:: LAMMPS

   # Water dimer - hydrogen bonding
   units real
   atom_style atomic
   
   region box block -15 15 -15 15 -15 15
   create_box 2 box
   
   # First water molecule
   create_atoms 2 single 0.0 0.0 0.0
   create_atoms 1 single 0.757 0.586 0.0
   create_atoms 1 single -0.757 0.586 0.0
   
   # Second water molecule (hydrogen bonded)
   create_atoms 2 single 2.8 0.0 0.0
   create_atoms 1 single 3.557 0.586 0.0
   create_atoms 1 single 2.043 0.586 0.0
   
   mass 1 1.00794   # H
   mass 2 15.9994   # O
   
   # BLYP with D3 dispersion for H-bonding
   pair_style dft BLYP def2-tzvp.HO.json dispersion D3 grid 100000
   pair_coeff 1 1 1.0 0.5   # H
   pair_coeff 2 2 8.0 1.52  # O
   
   # Optimize dimer
   minimize 1e-10 1e-12 2000 20000
   
   # Calculate binding energy
   compute pe all pe
   variable ebind equal "c_pe"
   print "Binding energy: ${ebind} kcal/mol"

----------

Tips for DFT calculations
-------------------------

**Basis set selection:**

* Use minimal basis (STO-3G) for initial structures
* Use double-zeta (def2-SVP, 6-31G) for optimization
* Use triple-zeta (def2-TZVP, cc-pVTZ) for accurate energies
* Add polarization for better geometries (6-31G*, def2-TZVPP)
* Add diffuse functions for anions (aug-cc-pVDZ)

**Functional selection:**

* LDA: Very fast but less accurate (metals)
* PBE: Good general-purpose GGA functional
* B3LYP: Popular hybrid, good for organic molecules
* PBE0: Robust hybrid functional
* SCAN: Accurate meta-GGA, good for diverse systems
* wB97M-V: Most accurate for thermochemistry

**Grid settings:**

* 30000: Quick and dirty calculations
* 50000: Default, good for optimization
* 75000: Production calculations
* 100000+: High accuracy, benchmarks

**SCF convergence:**

* tol 1e-7: Quick optimization
* tol 1e-8: Default
* tol 1e-9: Production
* tol 1e-10: High accuracy
* tol 1e-11: Benchmark quality

**Performance optimization:**

* Start with small basis and coarse grid
* Use pure functionals (PBE) over hybrids (B3LYP) when possible
* Increase maxiter if SCF struggles (default 100)
* Use DIIS (automatic) for better convergence
* Consider dispersion only for non-covalent interactions

----------

Downloading basis sets
----------------------

Use the provided Python script:

.. code-block:: bash

   # Download specific basis for elements
   python download_basis.py def2-svp H C N O
   
   # Common organic molecule basis sets
   python download_basis.py 6-31g H C N O S P
   python download_basis.py 6-31gs H C N O F Cl
   python download_basis.py cc-pvdz H C N O F
   
   # Transition metals
   python download_basis.py def2-tzvp Fe Co Ni Cu Zn
   python download_basis.py lanl2dz Ag Au Pd Pt
   
   # List all available functionals
   python list_functionals.py

----------

Common errors and solutions
---------------------------

**"Unknown XC functional"**

* Check functional name spelling
* Use common name (PBE) or full LibXC name (GGA_X_PBE+GGA_C_PBE)
* Verify LibXC installation

**"Cannot open basis set file"**

* Ensure JSON file exists in working directory
* Download from Basis Set Exchange
* Check all elements are included in basis file

**"SCF not converged"**

* Increase maxiter (e.g., maxiter 200)
* Reduce mixing parameter in code
* Start with smaller basis set
* Check molecular geometry is reasonable

**"Electron density exceeds tabulated range"**

* Atoms too close together
* Increase box size or adjust initial positions
* Check units (should be real for Angstroms)

**Memory issues**

* Reduce basis set size
* Decrease grid points
* Use fewer atoms or smaller system
