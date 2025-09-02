/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   DFT pair style using NWChemEx SCF module
------------------------------------------------------------------------- */

#include "pair_dft.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"
#include "utils.h"
#include "fmt/format.h"

// NWChemEx includes
#include <pluginplay/pluginplay.hpp>
#include <simde/simde.hpp>
#include <chemist/chemist.hpp>
#include <scf/scf.hpp>

#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairDFT::PairDFT(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;  // Single() not meaningful for DFT
  restartinfo = 0;    // No restart info for now
  one_coeff = 1;      // Only one coeff command
  manybody_flag = 1;  // DFT is many-body
  
  allocated = 0;
  
  // Default parameters
  functional_name = "PBE";
  basis_name = "def2-svp";
  max_scf_iterations = 100;
  energy_tolerance = 1e-8;
  density_tolerance = 1e-6;
  
  cut_global = 20.0;  // Default cutoff in Angstroms (not really used in DFT)
  
  // Initialize NWChemEx modules
  initialize_modules();
}

/* ---------------------------------------------------------------------- */

PairDFT::~PairDFT()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
  }
}

/* ----------------------------------------------------------------------
   Compute DFT energy and forces
------------------------------------------------------------------------- */

void PairDFT::compute(int eflag, int vflag)
{
  if (eflag || vflag) ev_setup(eflag, vflag);
  else evflag = vflag_fdotr = 0;

  // Convert LAMMPS atoms to NWChemEx molecule
  setup_molecule_from_lammps();
  
  // Setup basis set
  setup_basis_set();
  
  // Run SCF calculation
  run_scf_calculation();
  
  // Convert energy to LAMMPS units and store
  convert_energy_to_lammps_units();
  
  // Compute forces
  compute_forces();
  
  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   Global settings for DFT
------------------------------------------------------------------------- */

void PairDFT::settings(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR, "Illegal pair_style dft command");
  
  // Parse functional name
  functional_name = arg[0];
  
  // Parse basis set name
  basis_name = arg[1];
  
  // Parse optional arguments
  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "maxiter") == 0) {
      if (iarg + 1 >= narg) error->all(FLERR, "Illegal pair_style dft command");
      max_scf_iterations = atoi(arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "etol") == 0) {
      if (iarg + 1 >= narg) error->all(FLERR, "Illegal pair_style dft command");
      energy_tolerance = atof(arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "dtol") == 0) {
      if (iarg + 1 >= narg) error->all(FLERR, "Illegal pair_style dft command");
      density_tolerance = atof(arg[iarg + 1]);
      iarg += 2;
    } else {
      error->all(FLERR, fmt::format("Unknown pair_style dft option: {}", arg[iarg]));
    }
  }
  
  // Allocate arrays
  allocate();
  
  if (comm->me == 0) {
    utils::logmesg(lmp, "\n");
    utils::logmesg(lmp, "DFT Settings:\n");
    utils::logmesg(lmp, fmt::format("  Functional: {}\n", functional_name));
    utils::logmesg(lmp, fmt::format("  Basis set: {}\n", basis_name));
    utils::logmesg(lmp, fmt::format("  Max SCF iterations: {}\n", max_scf_iterations));
    utils::logmesg(lmp, fmt::format("  Energy tolerance: {:.2e}\n", energy_tolerance));
    utils::logmesg(lmp, fmt::format("  Density tolerance: {:.2e}\n", density_tolerance));
    utils::logmesg(lmp, "\n");
  }
}

/* ----------------------------------------------------------------------
   Set coefficients for one or more type pairs
------------------------------------------------------------------------- */

void PairDFT::coeff(int narg, char **arg)
{
  if (!allocated) allocate();
  
  if (narg != 3) error->all(FLERR, "Incorrect args for pair coefficients");
  
  // For DFT, we don't really use pair coefficients
  // Just set flags to indicate types are defined
  
  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);
  
  double cut_one = utils::numeric(FLERR, arg[2], false, lmp);
  
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }
  
  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   Initialize for a run
------------------------------------------------------------------------- */

void PairDFT::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR, "Pair style dft requires atom IDs");
  
  // Request a full neighbor list
  neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ----------------------------------------------------------------------
   Initialize one type pair
------------------------------------------------------------------------- */

double PairDFT::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    cut[i][j] = mix_distance(cut[i][i], cut[j][j]);
  }
  
  return cut[i][j];
}

/* ----------------------------------------------------------------------
   Allocate arrays
------------------------------------------------------------------------- */

void PairDFT::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  
  memory->create(setflag, n+1, n+1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;
  
  memory->create(cutsq, n+1, n+1, "pair:cutsq");
  memory->create(cut, n+1, n+1, "pair:cut");
}

/* ----------------------------------------------------------------------
   Write/read restart information
------------------------------------------------------------------------- */

void PairDFT::write_restart(FILE *fp)
{
  write_restart_settings(fp);
  
  int i, j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) fwrite(&cut[i][j], sizeof(double), 1, fp);
    }
}

void PairDFT::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();
  
  int i, j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j]) {
        if (me == 0) utils::sfread(FLERR, &cut[i][j], sizeof(double), 1, fp, nullptr, error);
        MPI_Bcast(&cut[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
}

void PairDFT::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global, sizeof(double), 1, fp);
  
  // Write functional and basis info
  int len = functional_name.length() + 1;
  fwrite(&len, sizeof(int), 1, fp);
  fwrite(functional_name.c_str(), sizeof(char), len, fp);
  
  len = basis_name.length() + 1;
  fwrite(&len, sizeof(int), 1, fp);
  fwrite(basis_name.c_str(), sizeof(char), len, fp);
  
  fwrite(&max_scf_iterations, sizeof(int), 1, fp);
  fwrite(&energy_tolerance, sizeof(double), 1, fp);
  fwrite(&density_tolerance, sizeof(double), 1, fp);
}

void PairDFT::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    utils::sfread(FLERR, &cut_global, sizeof(double), 1, fp, nullptr, error);
    
    int len;
    utils::sfread(FLERR, &len, sizeof(int), 1, fp, nullptr, error);
    char *buf = new char[len];
    utils::sfread(FLERR, buf, sizeof(char), len, fp, nullptr, error);
    functional_name = std::string(buf);
    delete[] buf;
    
    utils::sfread(FLERR, &len, sizeof(int), 1, fp, nullptr, error);
    buf = new char[len];
    utils::sfread(FLERR, buf, sizeof(char), len, fp, nullptr, error);
    basis_name = std::string(buf);
    delete[] buf;
    
    utils::sfread(FLERR, &max_scf_iterations, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &energy_tolerance, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &density_tolerance, sizeof(double), 1, fp, nullptr, error);
  }
  
  MPI_Bcast(&cut_global, 1, MPI_DOUBLE, 0, world);
  // Broadcast strings and other parameters...
}

/* ----------------------------------------------------------------------
   Single energy calculation (not implemented for DFT)
------------------------------------------------------------------------- */

double PairDFT::single(int i, int j, int itype, int jtype,
                       double rsq, double factor_coul, double factor_lj,
                       double &fforce)
{
  error->all(FLERR, "Single() not implemented for pair style dft");
  return 0.0;
}

/* ----------------------------------------------------------------------
   Extract data
------------------------------------------------------------------------- */

void *PairDFT::extract(const char *str, int &dim)
{
  dim = 0;
  
  if (strcmp(str, "dft_energy") == 0) return (void *) &total_energy;
  if (strcmp(str, "xc_energy") == 0) return (void *) &xc_energy;
  if (strcmp(str, "kinetic_energy") == 0) return (void *) &kinetic_energy;
  if (strcmp(str, "nuclear_repulsion") == 0) return (void *) &nuclear_repulsion;
  if (strcmp(str, "n_electrons") == 0) return (void *) &n_electrons;
  if (strcmp(str, "n_basis") == 0) return (void *) &n_basis_functions;
  if (strcmp(str, "converged") == 0) return (void *) &scf_converged;
  
  return nullptr;
}

/* ----------------------------------------------------------------------
   Convert energy to LAMMPS units
------------------------------------------------------------------------- */

void PairDFT::convert_energy_to_lammps_units()
{
  double conversion = 1.0;
  
  if (strcmp(update->unit_style, "metal") == 0) {
    conversion = HARTREE_TO_EV;
  } else if (strcmp(update->unit_style, "real") == 0) {
    conversion = HARTREE_TO_KCALMOL;
  }
  
  eng_vdotr = total_energy * conversion;
}

/* ----------------------------------------------------------------------
   Initialize NWChemEx modules
------------------------------------------------------------------------- */

void PairDFT::initialize_modules()
{
  // Initialize module manager
  module_manager = std::make_shared<pluginplay::ModuleManager>();
  
  // Load SCF modules
  scf::load_modules(*module_manager);
  scf::xc::gauxc::load_modules(*module_manager);
  
  // Set defaults
  scf::set_defaults(*module_manager);
  scf::xc::gauxc::set_defaults(*module_manager);
}

/* ----------------------------------------------------------------------
   Setup molecule from LAMMPS atoms
------------------------------------------------------------------------- */

void PairDFT::setup_molecule_from_lammps()
{
  molecule = std::make_unique<chemist::Molecule>();
  
  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  
  n_electrons = 0;
  
  for (int i = 0; i < nlocal; i++) {
    // Map LAMMPS type to atomic number
    // TODO: Implement proper type mapping
    int Z = type[i];  // For now, assume type = atomic number
    double mass = atom->mass[type[i]];
    
    // Convert position to Bohr
    double x_bohr = x[i][0] * ANGSTROM_TO_BOHR;
    double y_bohr = x[i][1] * ANGSTROM_TO_BOHR;
    double z_bohr = x[i][2] * ANGSTROM_TO_BOHR;
    
    molecule->push_back(chemist::Atom(Z, mass, x_bohr, y_bohr, z_bohr));
    n_electrons += Z;
  }
}

/* ----------------------------------------------------------------------
   Setup basis set
------------------------------------------------------------------------- */

void PairDFT::setup_basis_set()
{
  // Apply basis set to molecule
  try {
    auto aos = chemist::apply_basis(basis_name, *molecule);
    n_basis_functions = aos.size();
    basis_set = std::make_unique<simde::type::aos>(aos);
  } catch (const std::exception& e) {
    error->all(FLERR, fmt::format("Failed to apply basis {}: {}", basis_name, e.what()));
  }
}

/* ----------------------------------------------------------------------
   Run SCF calculation
------------------------------------------------------------------------- */

void PairDFT::run_scf_calculation()
{
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
  auto& rscf_mod = module_manager->at("Restricted SCF");
  
  // Create input for SCF
  using rscf_pt = simde::RestrictedSCF;
  auto rscf_input = rscf_mod.make_input<rscf_pt>();
  
  // Create Hamiltonian
  simde::type::MolecularHamiltonian H(*molecule, *basis_set);
  H.xc_functional(func_enum);
  
  // Create initial guess (core Hamiltonian)
  simde::type::rscf_wf guess(*basis_set);
  
  // Set input values
  rscf_input.at("Hamiltonian") = H;
  rscf_input.at("Initial Guess") = guess;
  
  // Run SCF
  if (comm->me == 0) {
    utils::logmesg(lmp, "\n=== SCF Calculation ===\n");
    utils::logmesg(lmp, fmt::format("Functional: {}\n", functional_name));
    utils::logmesg(lmp, fmt::format("Basis: {} ({} functions)\n", basis_name, n_basis_functions));
    utils::logmesg(lmp, fmt::format("Electrons: {}\n", n_electrons));
  }
  
  auto result = rscf_mod.run_as<simde::type::rscf_wf>(rscf_input);
  
  // Extract energy
  total_energy = result.energy();
  scf_converged = true;  // TODO: get convergence status from result
  
  if (comm->me == 0) {
    utils::logmesg(lmp, fmt::format("SCF Energy: {:.8f} Hartree\n", total_energy));
  }
}

/* ----------------------------------------------------------------------
   Compute forces from SCF
------------------------------------------------------------------------- */

void PairDFT::compute_forces()
{
  double **x = atom->x;
  double **f = atom->f;
  int nlocal = atom->nlocal;
  
  // For now, use numerical differentiation
  // TODO: Implement analytical gradients
  double delta = 1e-5;  // Displacement in Angstrom
  
  for (int i = 0; i < nlocal; i++) {
    for (int dim = 0; dim < 3; dim++) {
      // Save original position
      double x_orig = x[i][dim];
      
      // Forward displacement
      x[i][dim] = x_orig + delta;
      setup_molecule_from_lammps();
      setup_basis_set();
      
      // Run SCF for displaced geometry
      auto& rscf_mod = module_manager->at("Restricted SCF");
      using rscf_pt = simde::RestrictedSCF;
      auto input_plus = rscf_mod.make_input<rscf_pt>();
      
      using xc_func = chemist::qm_operator::xc_functional;
      xc_func func_enum;
      if (functional_name == "PBE") func_enum = xc_func::PBE;
      else if (functional_name == "BLYP") func_enum = xc_func::BLYP;
      else if (functional_name == "B3LYP") func_enum = xc_func::B3LYP;
      else if (functional_name == "PBE0") func_enum = xc_func::PBE0;
      else if (functional_name == "SVWN3") func_enum = xc_func::SVWN3;
      else if (functional_name == "SVWN5") func_enum = xc_func::SVWN5;
      
      simde::type::MolecularHamiltonian H_plus(*molecule, *basis_set);
      H_plus.xc_functional(func_enum);
      simde::type::rscf_wf guess_plus(*basis_set);
      
      input_plus.at("Hamiltonian") = H_plus;
      input_plus.at("Initial Guess") = guess_plus;
      auto result_plus = rscf_mod.run_as<simde::type::rscf_wf>(input_plus);
      double e_plus = result_plus.energy();
      
      // Backward displacement
      x[i][dim] = x_orig - delta;
      setup_molecule_from_lammps();
      setup_basis_set();
      
      auto input_minus = rscf_mod.make_input<rscf_pt>();
      simde::type::MolecularHamiltonian H_minus(*molecule, *basis_set);
      H_minus.xc_functional(func_enum);
      simde::type::rscf_wf guess_minus(*basis_set);
      
      input_minus.at("Hamiltonian") = H_minus;
      input_minus.at("Initial Guess") = guess_minus;
      auto result_minus = rscf_mod.run_as<simde::type::rscf_wf>(input_minus);
      double e_minus = result_minus.energy();
      
      // Restore position
      x[i][dim] = x_orig;
      
      // Compute force (negative gradient)
      double force_hartree_bohr = -(e_plus - e_minus) / (2.0 * delta * ANGSTROM_TO_BOHR);
      
      // Convert to LAMMPS units
      convert_forces_to_lammps_units();
      double force_conversion = 1.0;
      if (strcmp(update->unit_style, "metal") == 0) {
        force_conversion = HARTREE_BOHR_TO_EV_ANGSTROM;
      } else if (strcmp(update->unit_style, "real") == 0) {
        force_conversion = HARTREE_BOHR_TO_KCALMOL_ANGSTROM;
      }
      
      f[i][dim] = force_hartree_bohr * force_conversion;
    }
  }
}

/* ----------------------------------------------------------------------
   Convert forces to LAMMPS units
------------------------------------------------------------------------- */

void PairDFT::convert_forces_to_lammps_units()
{
  // This is handled inline in compute_forces() for now
  // Could be refactored to work on the entire force array
}
