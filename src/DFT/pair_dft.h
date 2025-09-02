/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   DFT pair style using NWChemEx SCF
------------------------------------------------------------------------- */

#ifndef LMP_PAIR_DFT_H
#define LMP_PAIR_DFT_H

#include "pair.h"
#include <memory>
#include <string>

namespace pluginplay {
class ModuleManager;
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
  // SCF module manager
  std::shared_ptr<pluginplay::ModuleManager> mm;
  
  // Parameters
  std::string functional_name;
  std::string basis_name;
  int max_scf_iter;
  double e_tol;
  double d_tol;
  
  // Results
  double total_energy;
  double xc_energy;
  int n_electrons;
  int n_basis;
  bool converged;
  
  // Internal methods
  void run_scf();
  void allocate();
  int allocated;
  double cut_global;
  double **cut;
};

}

#endif
