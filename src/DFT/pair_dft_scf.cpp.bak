/* ----------------------------------------------------------------------
   SCF calculation interface for LAMMPS DFT
------------------------------------------------------------------------- */

#include "pair_dft.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "update.h"
#include "force.h"

// NWChemEx includes
#include <pluginplay/pluginplay.hpp>
#include <simde/simde.hpp>
#include <chemist/chemist.hpp>
#include <scf/scf.hpp>

using namespace LAMMPS_NS;

static constexpr double ANGSTROM_TO_BOHR = 1.8897259886;
static constexpr double HARTREE_TO_EV = 27.211386245988;
static constexpr double HARTREE_TO_KCALMOL = 627.5094740631;
static constexpr double HARTREE_BOHR_TO_EV_ANGSTROM = HARTREE_TO_EV / ANGSTROM_TO_BOHR;

void PairDFT::run_scf()
{
  // Initialize module manager if needed
  if (!mm) {
    mm = std::make_shared<pluginplay::ModuleManager>();
    
    // Load SCF modules
    scf::load_modules(*mm);
    scf::xc::gauxc::load_modules(*mm);
    
    // Set defaults
    scf::set_defaults(*mm);
    scf::xc::gauxc::set_defaults(*mm);
  }
  
  // Create molecule from LAMMPS atoms
  chemist::Molecule mol;
  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  
  n_electrons = 0;
  for (int i = 0; i < nlocal; i++) {
    // Map LAMMPS type to atomic number (for now assume H=1)
    int Z = 1;  // TODO: proper type mapping
    double mass = 1.00794;
    
    // Convert position to Bohr
    double x_bohr = x[i][0] * ANGSTROM_TO_BOHR;
    double y_bohr = x[i][1] * ANGSTROM_TO_BOHR;
    double z_bohr = x[i][2] * ANGSTROM_TO_BOHR;
    
    mol.push_back(chemist::Atom(Z, mass, x_bohr, y_bohr, z_bohr));
    n_electrons += Z;
  }
  
  // Apply basis set
  chemist::AOBasisSet<double> aos;
  try {
    // Use basis_name directly - SCF knows how to handle standard basis names
    aos = chemist::apply_basis(basis_name, mol);
    n_basis = aos.size();
  } catch (const std::exception& e) {
    error->all(FLERR, fmt::format("Failed to apply basis {}: {}", basis_name, e.what()));
  }
  
  // Set up XC functional
  using xc_func = chemist::qm_operator::xc_functional;
  xc_func func_enum;
  
  // Map functional name to enum
  if (functional_name == "PBE") func_enum = xc_func::PBE;
  else if (functional_name == "BLYP") func_enum = xc_func::BLYP;
  else if (functional_name == "B3LYP") func_enum = xc_func::B3LYP;
  else if (functional_name == "PBE0") func_enum = xc_func::PBE0;
  else if (functional_name == "SVWN3") func_enum = xc_func::SVWN3;
  else if (functional_name == "SVWN5") func_enum = xc_func::SVWN5;
  else error->all(FLERR, fmt::format("Unknown functional: {}", functional_name));
  
  // Get the restricted SCF module
  auto& rscf_mod = mm->at("Restricted SCF");
  
  // Create input for SCF
  using rscf_pt = simde::RestrictedSCF;
  auto rscf_input = rscf_mod.make_input<rscf_pt>();
  
  // Create Hamiltonian
  simde::type::MolecularHamiltonian H(mol, aos);
  H.xc_functional(func_enum);
  
  // Create initial guess (core Hamiltonian)
  simde::type::rscf_wf guess(aos);
  
  // Set input values
  rscf_input.at("Hamiltonian") = H;
  rscf_input.at("Initial Guess") = guess;
  
  // Run SCF
  if (comm->me == 0) {
    utils::logmesg(lmp, "\n=== SCF Calculation ===\n");
    utils::logmesg(lmp, fmt::format("Functional: {}\n", functional_name));
    utils::logmesg(lmp, fmt::format("Basis: {} ({} functions)\n", basis_name, n_basis));
    utils::logmesg(lmp, fmt::format("Electrons: {}\n", n_electrons));
  }
  
  auto result = rscf_mod.run_as<simde::type::rscf_wf>(rscf_input);
  
  // Extract energy
  total_energy = result.energy();
  converged = true;  // TODO: get convergence status from result
  
  // Convert energy to LAMMPS units
  double energy_conversion = 1.0;
  if (strcmp(update->unit_style, "metal") == 0) {
    energy_conversion = HARTREE_TO_EV;
  } else if (strcmp(update->unit_style, "real") == 0) {
    energy_conversion = HARTREE_TO_KCALMOL;
  }
  
  eng_vdotr = total_energy * energy_conversion;
  
  if (comm->me == 0) {
    utils::logmesg(lmp, fmt::format("SCF Energy: {:.8f} Hartree = {:.8f} {}\n",
                                    total_energy, eng_vdotr,
                                    (strcmp(update->unit_style, "metal") == 0) ? "eV" : "kcal/mol"));
  }
  
  // Compute forces (numerical for now)
  double **f = atom->f;
  double delta = 1e-5;  // Displacement in Angstrom
  
  for (int i = 0; i < nlocal; i++) {
    for (int dim = 0; dim < 3; dim++) {
      // Save original position
      double x_orig = x[i][dim];
      
      // Forward displacement
      x[i][dim] = x_orig + delta;
      
      // Build molecule with displaced atom
      chemist::Molecule mol_plus;
      for (int j = 0; j < nlocal; j++) {
        int Z = 1;
        double mass = 1.00794;
        double xb = x[j][0] * ANGSTROM_TO_BOHR;
        double yb = x[j][1] * ANGSTROM_TO_BOHR;
        double zb = x[j][2] * ANGSTROM_TO_BOHR;
        mol_plus.push_back(chemist::Atom(Z, mass, xb, yb, zb));
      }
      
      // Run SCF for displaced geometry
      chemist::AOBasisSet<double> aos_plus = chemist::apply_basis(basis_name, mol_plus);
      simde::type::MolecularHamiltonian H_plus(mol_plus, aos_plus);
      H_plus.xc_functional(func_enum);
      simde::type::rscf_wf guess_plus(aos_plus);
      
      auto input_plus = rscf_mod.make_input<rscf_pt>();
      input_plus.at("Hamiltonian") = H_plus;
      input_plus.at("Initial Guess") = guess_plus;
      auto result_plus = rscf_mod.run_as<simde::type::rscf_wf>(input_plus);
      double e_plus = result_plus.energy();
      
      // Backward displacement
      x[i][dim] = x_orig - delta;
      
      chemist::Molecule mol_minus;
      for (int j = 0; j < nlocal; j++) {
        int Z = 1;
        double mass = 1.00794;
        double xb = x[j][0] * ANGSTROM_TO_BOHR;
        double yb = x[j][1] * ANGSTROM_TO_BOHR;
        double zb = x[j][2] * ANGSTROM_TO_BOHR;
        mol_minus.push_back(chemist::Atom(Z, mass, xb, yb, zb));
      }
      
      chemist::AOBasisSet<double> aos_minus = chemist::apply_basis(basis_name, mol_minus);
      simde::type::MolecularHamiltonian H_minus(mol_minus, aos_minus);
      H_minus.xc_functional(func_enum);
      simde::type::rscf_wf guess_minus(aos_minus);
      
      auto input_minus = rscf_mod.make_input<rscf_pt>();
      input_minus.at("Hamiltonian") = H_minus;
      input_minus.at("Initial Guess") = guess_minus;
      auto result_minus = rscf_mod.run_as<simde::type::rscf_wf>(input_minus);
      double e_minus = result_minus.energy();
      
      // Restore position
      x[i][dim] = x_orig;
      
      // Compute force (negative gradient)
      double force_hartree_bohr = -(e_plus - e_minus) / (2.0 * delta * ANGSTROM_TO_BOHR);
      
      // Convert to LAMMPS units
      double force_conversion = 1.0;
      if (strcmp(update->unit_style, "metal") == 0) {
        force_conversion = HARTREE_BOHR_TO_EV_ANGSTROM;
      } else if (strcmp(update->unit_style, "real") == 0) {
        force_conversion = HARTREE_TO_KCALMOL / ANGSTROM_TO_BOHR;
      }
      
      f[i][dim] = force_hartree_bohr * force_conversion;
    }
  }
}
