/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   DFT pair style using NWChemEx SCF module
------------------------------------------------------------------------- */

#ifndef LMP_PAIR_DFT_H
#define LMP_PAIR_DFT_H

#include "pair.h"
#include <memory>
#include <string>

// Forward declarations
namespace pluginplay {
class ModuleManager;
}

namespace chemist {
class Molecule;
template<typename T> class AOBasisSet;
}

namespace simde {
namespace type {
class aos;
class chemical_system;
class hamiltonian;
class rscf_wf;
}
}

namespace LAMMPS_NS {

class PairDFT : public Pair {
public:
  PairDFT(class LAMMPS *);
  ~PairDFT() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  double single(int, int, int, int, double, double, double, double &) override;
  void *extract(const char *, int &) override;

protected:
  // NWChemEx module manager
  std::shared_ptr<pluginplay::ModuleManager> module_manager;
  
  // Current molecular system
  std::unique_ptr<chemist::Molecule> molecule;
  std::unique_ptr<simde::type::aos> basis_set;
  std::unique_ptr<simde::type::chemical_system> chem_system;
  
  // SCF parameters
  std::string functional_name;
  std::string basis_name;
  int max_scf_iterations;
  double energy_tolerance;
  double density_tolerance;
  
  // Results from last SCF calculation
  double total_energy;       // In Hartree
  double xc_energy;          // In Hartree
  double kinetic_energy;     // In Hartree
  double nuclear_repulsion;  // In Hartree
  int n_electrons;
  int n_basis_functions;
  bool scf_converged;
  
  // Internal methods
  void initialize_modules();
  void setup_molecule_from_lammps();
  void setup_basis_set();
  void run_scf_calculation();
  void compute_forces();
  void convert_energy_to_lammps_units();
  void convert_forces_to_lammps_units();
  
  // Utilities
  void allocate();
  int allocated;
  double cut_global;
  double **cut;
  
  // Unit conversion constants
  static constexpr double ANGSTROM_TO_BOHR = 1.8897259886;
  static constexpr double HARTREE_TO_EV = 27.211386245988;
  static constexpr double HARTREE_TO_KCALMOL = 627.5094740631;
  static constexpr double HARTREE_BOHR_TO_EV_ANGSTROM = HARTREE_TO_EV / ANGSTROM_TO_BOHR;
  static constexpr double HARTREE_BOHR_TO_KCALMOL_ANGSTROM = HARTREE_TO_KCALMOL / ANGSTROM_TO_BOHR;
};

}  // namespace LAMMPS_NS

#endif
