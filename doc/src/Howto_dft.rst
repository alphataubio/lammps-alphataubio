Using DFT in LAMMPS simulations
================================

The :doc:`pair_style dft <pair_dft>` command enables quantum mechanical 
density functional theory (DFT) calculations within LAMMPS molecular dynamics
simulations. This allows you to compute electronic structure, energies, and
forces at the quantum level while leveraging LAMMPS's powerful MD capabilities.
This document explains how to effectively use DFT in your simulations, from
basic setup to advanced techniques.

Overview of DFT in molecular simulations
-----------------------------------------

DFT provides a quantum mechanical description of electronic structure by solving
the Kohn-Sham equations self-consistently. Unlike classical force fields that
use predetermined parameters, DFT computes energies and forces from first 
principles based on electron density. This makes DFT particularly valuable for:

* Systems where electronic structure changes during dynamics (bond breaking/forming)
* Novel materials where force field parameters are unavailable
* High-accuracy calculations of molecular properties
* Benchmarking and validating classical force fields
* Studying chemical reactions and catalysis

The trade-off is computational cost - DFT calculations are typically 
1000-10000 times more expensive than classical force fields. Therefore,
careful consideration of system size, basis sets, and functionals is essential.

Setting up a basic DFT calculation
-----------------------------------

To perform a DFT calculation in LAMMPS, you need three key components:

1. **A functional** - Defines how to compute exchange-correlation energy
2. **A basis set** - Mathematical functions to represent molecular orbitals  
3. **Atomic charges** - Nuclear charges for each atom type

Let's start with a simple hydrogen molecule calculation to illustrate
the basic setup:

.. code-block:: LAMMPS

   # Simple H2 molecule with DFT
   units real
   atom_style atomic
   
   # Create simulation box and atoms
   region box block -10 10 -10 10 -10 10
   create_box 1 box
   create_atoms 1 single 0.0 0.0 -0.37  # H atom 1
   create_atoms 1 single 0.0 0.0  0.37  # H atom 2
   mass 1 1.00794
   
   # Setup DFT calculation
   pair_style dft PBE def2-svp.H.json
   pair_coeff 1 1 1.0 0.5  # H: nuclear charge=1, vdW radius=0.5
   
   # Optimize geometry
   minimize 1e-8 1e-10 1000 10000
   
   # Output results
   write_data h2_optimized.data

This example uses the PBE functional (a popular GGA functional) with the
def2-SVP basis set. The pair_coeff command specifies that hydrogen atoms
have nuclear charge 1 and a van der Waals radius of 0.5 Å (used for
dispersion corrections if enabled).

Obtaining basis sets
--------------------

Basis sets define the mathematical functions used to represent molecular
orbitals. The DFT package requires basis sets in JSON format from the
`Basis Set Exchange <https://www.basissetexchange.org/>`_. A Python script
is provided to download basis sets:

.. code-block:: bash

   # Download def2-SVP basis for H, C, N, O atoms
   python download_basis.py def2-svp H C N O
   
   # This creates def2-svp.HCNO.json

Common basis sets and their uses:

* **Minimal basis sets** (STO-3G): Very fast, useful for initial structures
* **Double-zeta** (def2-SVP, 6-31G): Good balance of speed and accuracy
* **Triple-zeta** (def2-TZVP, cc-pVTZ): High accuracy for production
* **Polarized** (6-31G\*, def2-TZVPP): Better for molecular geometries
* **Diffuse functions** (aug-cc-pVDZ): Important for anions and weak interactions

Choosing the right functional
------------------------------

The choice of exchange-correlation functional significantly impacts both
accuracy and computational cost. Different functionals excel at different
properties:

**Local Density Approximation (LDA)**
   Fastest but least accurate. Suitable for metals and quick estimates:

   .. code-block:: LAMMPS

      pair_style dft LDA sto-3g.json

**Generalized Gradient Approximation (GGA)**
   Good balance of speed and accuracy. PBE and BLYP are popular choices:

   .. code-block:: LAMMPS

      pair_style dft PBE def2-svp.json      # General purpose
      pair_style dft BLYP def2-svp.json     # Often better for organics

**Hybrid functionals**
   Include exact exchange, more accurate but 3-5x slower. B3LYP is the
   most popular for organic molecules:

   .. code-block:: LAMMPS

      pair_style dft B3LYP 6-31gs.json      # 20% exact exchange
      pair_style dft PBE0 def2-svp.json     # 25% exact exchange

**Meta-GGA functionals**
   Include kinetic energy density, very accurate for diverse systems:

   .. code-block:: LAMMPS

      pair_style dft SCAN def2-tzvp.json    # Excellent for many properties
      pair_style dft TPSS def2-svp.json     # Good and relatively fast

**Range-separated functionals**
   Best for charge-transfer and excited states:

   .. code-block:: LAMMPS

      pair_style dft wB97M-V def2-qzvpp.json  # State-of-the-art accuracy

Understanding computational cost
---------------------------------

DFT calculations scale steeply with system size and basis set:

* **Basis functions**: Scales as O(N⁴) for hybrids, O(N³) for pure functionals
* **Grid points**: Linear scaling but large prefactor
* **SCF iterations**: Typically 10-30 cycles to convergence

To manage computational cost:

1. **Start with small basis sets** for structure optimization
2. **Use pure functionals** (PBE, BLYP) when hybrids aren't necessary
3. **Reduce grid size** for optimization (30000-50000 points)
4. **Increase grid size** for final energies (75000-100000 points)

Here's an efficient workflow for expensive calculations:

.. code-block:: LAMMPS

   # Step 1: Rough optimization with minimal basis
   pair_style dft PBE sto-3g.json grid 30000 tol 1e-7
   pair_coeff * * 6.0 1.7
   minimize 1e-6 1e-8 100 500
   write_data rough.data
   
   # Step 2: Refine with better basis
   read_data rough.data
   pair_style dft PBE def2-svp.json grid 50000 tol 1e-8  
   minimize 1e-8 1e-10 500 2000
   write_data refined.data
   
   # Step 3: High-accuracy single point energy
   read_data refined.data
   pair_style dft PBE def2-tzvp.json grid 100000 tol 1e-10
   run 0  # Single point calculation

Grid integration and accuracy
------------------------------

Exchange-correlation energies for GGA and meta-GGA functionals require
numerical integration on a grid. The grid size dramatically affects
both accuracy and computational cost:

.. code-block:: LAMMPS

   # Coarse grid - fast but less accurate
   pair_style dft PBE basis.json grid 30000
   
   # Medium grid - default, good for most uses  
   pair_style dft PBE basis.json grid 50000
   
   # Fine grid - high accuracy
   pair_style dft PBE basis.json grid 100000
   
   # Ultra-fine - benchmark quality
   pair_style dft PBE basis.json grid 200000

The grid uses Lebedev-Laikov angular quadrature combined with 
Chebyshev-Gauss radial quadrature. Becke partitioning divides
space between atoms for multi-center integration.

SCF convergence control
------------------------

The self-consistent field (SCF) procedure iteratively solves for the
electron density. Convergence is controlled by tolerance and iteration
parameters:

.. code-block:: LAMMPS

   # Tight convergence for production
   pair_style dft B3LYP basis.json tol 1e-10 maxiter 200
   
   # Looser convergence for testing
   pair_style dft B3LYP basis.json tol 1e-7 maxiter 50

If SCF fails to converge:

* **Increase maxiter** - Some systems need more iterations
* **Start with a simpler functional** - Converge with PBE, then switch
* **Check molecular geometry** - Atoms too close can cause problems
* **Reduce mixing** - Modify mixing parameter in the source code

Including dispersion corrections
---------------------------------

DFT functionals typically underestimate van der Waals interactions.
Grimme's DFT-D methods add an empirical correction:

.. code-block:: LAMMPS

   # Add D3 dispersion with Becke-Johnson damping (recommended)
   pair_style dft PBE basis.json dispersion D3BJ
   
   # Original D3 with zero-damping
   pair_style dft PBE basis.json dispersion D3
   
   # Latest D4 with charge-dependent C6
   pair_style dft PBE basis.json dispersion D4

Dispersion corrections are essential for:

* Molecular crystals and polymers
* Host-guest complexes
* Protein-ligand binding
* Stacked aromatic systems
* Any system with significant vdW interactions

Working with water and hydrogen bonding
----------------------------------------

Water is challenging for DFT due to the importance of both covalent
and hydrogen bonding. Here's an optimized setup for water calculations:

.. code-block:: LAMMPS

   # Water dimer with hydrogen bonding
   units real
   atom_style atomic
   
   # Create two water molecules
   region box block -15 15 -15 15 -15 15  
   create_box 2 box
   
   # First water
   create_atoms 2 single 0.0 0.0 0.0      # O
   create_atoms 1 single 0.757 0.586 0.0  # H
   create_atoms 1 single -0.757 0.586 0.0 # H
   
   # Second water (H-bonded)
   create_atoms 2 single 2.8 0.0 0.0      # O
   create_atoms 1 single 3.557 0.586 0.0  # H  
   create_atoms 1 single 2.043 0.586 0.0  # H
   
   mass 1 1.00794   # H
   mass 2 15.9994   # O
   
   # BLYP-D3 works well for water
   pair_style dft BLYP def2-tzvp.json dispersion D3BJ grid 100000
   pair_coeff 1 1 1.0 0.5   # H
   pair_coeff 2 2 8.0 1.52  # O
   
   # Optimize with tight convergence
   minimize 1e-10 1e-12 2000 20000
   
   # Calculate interaction energy
   variable e_total equal pe
   print "Total energy: ${e_total} kcal/mol"

For liquid water simulations, consider:

* **Functionals**: BLYP-D3, PBE-D3, or SCAN work well
* **Basis sets**: At least def2-TZVP or aug-cc-pVDZ
* **Temperature effects**: Classical MD after DFT optimization

Studying organic molecules
--------------------------

Organic molecules benefit from specific functional choices. Here's
an example for benzene, showcasing aromatic systems:

.. code-block:: LAMMPS

   # Benzene with meta-GGA functional
   units real
   atom_style atomic
   boundary p p p
   
   region box block -15 15 -15 15 -15 15
   create_box 2 box
   
   # Create benzene ring (C6H6)
   # Carbon atoms in hexagon
   create_atoms 2 single  1.396  0.000 0.0
   create_atoms 2 single  0.698  1.209 0.0
   create_atoms 2 single -0.698  1.209 0.0  
   create_atoms 2 single -1.396  0.000 0.0
   create_atoms 2 single -0.698 -1.209 0.0
   create_atoms 2 single  0.698 -1.209 0.0
   
   # Hydrogen atoms
   create_atoms 1 single  2.476  0.000 0.0
   create_atoms 1 single  1.238  2.144 0.0
   create_atoms 1 single -1.238  2.144 0.0
   create_atoms 1 single -2.476  0.000 0.0
   create_atoms 1 single -1.238 -2.144 0.0
   create_atoms 1 single  1.238 -2.144 0.0
   
   mass 1 1.00794   # H
   mass 2 12.0107   # C
   
   # B3LYP is excellent for aromatics
   pair_style dft B3LYP def2-tzvp.json grid 100000 tol 1e-10
   pair_coeff 1 1 1.0 0.5   # H
   pair_coeff 2 2 6.0 1.7   # C
   
   # Tight optimization for planar structure
   min_style cg
   minimize 1e-11 1e-13 5000 50000
   
   # Verify planarity
   compute cz all property/atom z
   compute zmax all reduce max c_cz
   compute zmin all reduce min c_cz  
   variable flatness equal c_zmax-c_zmin
   print "Deviation from planarity: ${flatness} Angstroms"

Best practices for organics:

* **Aromatics**: B3LYP or wB97M-V with dispersion
* **Alkanes**: PBE-D3 or BLYP-D3
* **Biomolecules**: B3LYP-D3 or PBE0-D3
* **Reactions**: B3LYP or M06-2X for transition states

Understanding forces: Hellmann-Feynman and Pulay
-------------------------------------------------

A critical aspect of DFT with Gaussian basis sets is the calculation
of forces. Unlike plane-wave codes, atom-centered basis sets require
both Hellmann-Feynman and Pulay force contributions:

.. math::

   \mathbf{F}_{total} = \mathbf{F}_{HF} + \mathbf{F}_{Pulay}

The Hellmann-Feynman force arises from the derivative of the potential,
while the Pulay force accounts for the basis set moving with the atoms.
**The Pulay term is essential** - without it, forces can be wrong by
30-50%, leading to incorrect geometries and failed dynamics.

The DFT package automatically includes both contributions. This is why
DFT with Gaussian basis sets can give accurate forces despite using a
finite basis, unlike plane-wave methods that require very large cutoffs
for accurate forces.

Advanced example: Reaction pathway
-----------------------------------

Here's a complete example studying proton transfer in formamidine,
demonstrating how to use DFT for chemical reactions:

.. code-block:: LAMMPS

   # Formamidine tautomerization: HC(NH2)=NH <-> HC(NH)=NH2
   units real
   atom_style atomic
   boundary p p p
   
   # Read initial structure
   region box block -10 10 -10 10 -10 10
   create_box 3 box
   
   # Create formamidine HC(NH2)=NH
   # (coordinates would normally come from file)
   create_atoms 2 single 0.000  0.000  0.000  # C
   create_atoms 3 single 1.300  0.000  0.000  # N1  
   create_atoms 3 single -1.200 0.000  0.000  # N2
   create_atoms 1 single 0.000  0.900  0.000  # H-C
   create_atoms 1 single 1.800  0.800  0.000  # H-N1a
   create_atoms 1 single 1.800 -0.800  0.000  # H-N1b  
   create_atoms 1 single -2.000 0.000  0.000  # H-N2
   
   mass 1 1.00794   # H
   mass 2 12.0107   # C
   mass 3 14.0067   # N
   
   # Use high-level functional for reaction
   pair_style dft B3LYP def2-tzvpp.json grid 100000 tol 1e-10
   pair_coeff 1 1 1.0 0.5   # H
   pair_coeff 2 2 6.0 1.7   # C  
   pair_coeff 3 3 7.0 1.55  # N
   
   # Initial optimization
   minimize 1e-10 1e-12 1000 10000
   variable e_reactant equal pe
   write_data reactant.data
   
   # Set up reaction coordinate (H transfer from N1 to N2)
   group h_transfer id 5  # The transferring H
   group n_donor id 2     # N1 (donor)
   group n_accept id 3    # N2 (acceptor)
   
   # Constrained optimization along reaction path
   # Move H from N1 toward N2 in small steps
   variable step loop 20
   label reaction_loop
   
     # Calculate current position along path
     variable frac equal ${step}/20.0
     variable x_h equal 1.8-${frac}*3.0
     set atom 5 x ${x_h}  # Move H
     
     # Optimize with H position fixed
     fix freeze h_transfer setforce 0.0 0.0 0.0
     minimize 1e-8 1e-10 100 1000
     unfix freeze
     
     # Record energy
     variable e_path equal pe
     print "Step ${step}: E = ${e_path} kcal/mol"
     
   next step
   jump SELF reaction_loop
   
   # Final optimization of product
   minimize 1e-10 1e-12 1000 10000
   variable e_product equal pe
   write_data product.data
   
   # Calculate barrier
   variable barrier equal ${e_path_max}-${e_reactant}
   print "Reaction barrier: ${barrier} kcal/mol"

This example demonstrates:

* Setting up a chemical reaction coordinate
* Constrained optimization along a reaction path
* Calculating reaction barriers
* Using high-level functionals for accuracy

Performance optimization strategies
------------------------------------

DFT calculations are expensive, but several strategies can improve
performance:

**1. Hierarchical optimization**

Start with cheap methods and refine:

.. code-block:: LAMMPS

   # Stage 1: Minimal basis
   pair_style dft PBE sto-3g.json grid 20000 tol 1e-6
   minimize 1e-5 1e-7 50 200
   
   # Stage 2: Better basis
   pair_style dft PBE def2-svp.json grid 50000 tol 1e-8
   minimize 1e-8 1e-10 200 1000
   
   # Stage 3: Production quality
   pair_style dft PBE def2-tzvp.json grid 100000 tol 1e-10
   run 0  # Single point

**2. Mixed quantum/classical approaches**

For large systems, treat only the reactive region with DFT:

.. code-block:: LAMMPS

   # Define QM and MM regions
   group qm_atoms id 1:20     # Reactive site
   group mm_atoms id 21:1000  # Environment
   
   # Use hybrid potential (requires QM/MM implementation)
   # This is conceptual - full QM/MM requires additional setup

**3. Choosing efficient functionals**

For similar accuracy, computational cost varies significantly:

* LDA: Fastest (baseline)
* GGA (PBE, BLYP): 1.2-1.5x slower than LDA
* Meta-GGA (SCAN, TPSS): 2-3x slower than LDA  
* Hybrid (B3LYP, PBE0): 3-5x slower than LDA
* Double-hybrid: 10-20x slower than LDA

**4. Basis set considerations**

Computational scaling with basis size:

* STO-3G: ~10 functions per heavy atom
* 6-31G: ~15 functions per heavy atom
* def2-SVP: ~20 functions per heavy atom
* def2-TZVP: ~35 functions per heavy atom
* def2-QZVP: ~55 functions per heavy atom

Cost scales as N⁴ with number of basis functions for hybrids!

Troubleshooting common problems
--------------------------------

**SCF convergence failures**

If the self-consistent field fails to converge:

.. code-block:: LAMMPS

   # Increase iterations
   pair_style dft PBE basis.json maxiter 200
   
   # Loosen tolerance temporarily
   pair_style dft PBE basis.json tol 1e-7
   
   # Try a different starting guess (modify source code)
   # Or start from a converged similar structure

**"Electron density exceeds tabulated range"**

This error indicates atoms are too close:

.. code-block:: LAMMPS

   # Check for overlapping atoms
   delete_atoms overlap 0.5 all all
   
   # Increase box size
   change_box all x scale 1.1 y scale 1.1 z scale 1.1
   
   # Pre-optimize with classical potential
   pair_style lj/cut 10.0
   pair_coeff * * 0.1 3.0
   minimize 1e-4 1e-6 100 1000
   
   # Then switch to DFT
   pair_style dft PBE basis.json

**Memory issues**

DFT requires substantial memory for integral storage:

* Reduce basis set size
* Use fewer grid points
* Decrease system size
* Use pure functionals instead of hybrids

Validation and benchmarking
----------------------------

Always validate your DFT setup against known results:

.. code-block:: LAMMPS

   # Benchmark: H2 bond length should be ~0.74 Angstroms
   units real
   atom_style atomic
   
   region box block -10 10 -10 10 -10 10
   create_box 1 box
   create_atoms 1 single 0.0 0.0 -0.4
   create_atoms 1 single 0.0 0.0  0.4
   mass 1 1.00794
   
   # Test your functional/basis combination
   pair_style dft PBE cc-pvtz.json
   pair_coeff 1 1 1.0 0.5
   
   # Optimize
   minimize 1e-10 1e-12 1000 10000
   
   # Measure bond length
   compute bond all pair/local dist
   compute blen all reduce min c_bond
   variable blen equal c_blen
   print "H2 bond length: ${blen} Angstroms"
   print "Expected: ~0.74 Angstroms"
   
   # Calculate dissociation energy
   variable e_molecule equal pe
   delete_atoms group all
   create_atoms 1 single 0.0 0.0 0.0
   run 0
   variable e_atom equal pe
   variable de equal ${e_molecule}-2*${e_atom}
   print "Dissociation energy: ${de} kcal/mol"
   print "Expected: ~104 kcal/mol"

Compare with:

* Experimental values
* High-level quantum chemistry (CCSD(T))
* Published DFT benchmarks
* Basis set extrapolation

Summary and recommendations
----------------------------

DFT in LAMMPS provides powerful capabilities for quantum mechanical
simulations. Key points to remember:

1. **Start simple**: Use small basis sets and pure functionals for testing
2. **Choose functionals wisely**: B3LYP for organics, PBE for general use, SCAN for diverse systems
3. **Always include Pulay forces**: Automatically handled but essential for accuracy
4. **Use dispersion corrections**: Critical for non-covalent interactions
5. **Optimize workflow**: Hierarchical optimization saves computer time
6. **Validate results**: Benchmark against known systems
7. **Consider cost**: DFT is expensive - use judiciously

For production calculations:

* **Geometry optimization**: PBE/def2-SVP or B3LYP/6-31G*
* **Accurate energies**: B3LYP/def2-TZVP or SCAN/def2-TZVP
* **Weak interactions**: PBE-D3BJ/def2-TZVP or B3LYP-D3BJ/def2-TZVP
* **Best accuracy**: wB97M-V/def2-QZVPP (very expensive)

The combination of LAMMPS's MD capabilities with quantum mechanical
accuracy opens many possibilities for studying complex chemical systems,
reactions, and materials where electronic structure is crucial.

----------

References
""""""""""

For more information on DFT methods and their implementation:

* Kohn & Sham, Phys. Rev. 140, A1133 (1965) - Original DFT formulation
* Mardirossian & Head-Gordon, J. Chem. Phys. 144, 214110 (2016) - wB97M-V functional
* Grimme et al., J. Chem. Phys. 132, 154104 (2010) - DFT-D3 dispersion
* Pulay, Mol. Phys. 17, 197 (1969) - Pulay forces
* LibXC library: https://www.tddft.org/programs/libxc/
* Basis Set Exchange: https://www.basissetexchange.org/
