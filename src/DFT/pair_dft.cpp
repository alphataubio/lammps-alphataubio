/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   DFT pair style using NWChemEx SCF
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

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

PairDFT::PairDFT(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  
  allocated = 0;
  
  // Defaults
  functional_name = "PBE";
  basis_name = "def2-svp";
  max_scf_iter = 100;
  e_tol = 1e-8;
  d_tol = 1e-6;
  cut_global = 20.0;
}

PairDFT::~PairDFT()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
  }
}

void PairDFT::compute(int eflag, int vflag)
{
  if (eflag || vflag) ev_setup(eflag, vflag);
  else evflag = vflag_fdotr = 0;

  // Run SCF calculation
  run_scf();
  
  // Energy is set in run_scf()
  
  if (vflag_fdotr) virial_fdotr_compute();
}

void PairDFT::settings(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR, "Illegal pair_style dft command");
  
  functional_name = arg[0];
  basis_name = arg[1];
  
  // Parse optional arguments
  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "maxiter") == 0) {
      if (iarg + 1 >= narg) error->all(FLERR, "Illegal pair_style dft command");
      max_scf_iter = atoi(arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "etol") == 0) {
      if (iarg + 1 >= narg) error->all(FLERR, "Illegal pair_style dft command");
      e_tol = atof(arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "dtol") == 0) {
      if (iarg + 1 >= narg) error->all(FLERR, "Illegal pair_style dft command");
      d_tol = atof(arg[iarg + 1]);
      iarg += 2;
    } else {
      error->all(FLERR, fmt::format("Unknown dft option: {}", arg[iarg]));
    }
  }
  
  allocate();
  
  if (comm->me == 0) {
    utils::logmesg(lmp, fmt::format("DFT: functional={} basis={}\n", 
                                    functional_name, basis_name));
  }
}

void PairDFT::coeff(int narg, char **arg)
{
  if (!allocated) allocate();
  if (narg != 3) error->all(FLERR, "Incorrect args for pair coefficients");
  
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

void PairDFT::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR, "Pair style dft requires atom IDs");
  
  neighbor->add_request(this, NeighConst::REQ_FULL);
}

double PairDFT::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    cut[i][j] = mix_distance(cut[i][i], cut[j][j]);
  }
  return cut[i][j];
}

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

void PairDFT::write_restart(FILE *fp)
{
  write_restart_settings(fp);
  
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) fwrite(&cut[i][j], sizeof(double), 1, fp);
    }
}

void PairDFT::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();
  
  int me = comm->me;
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++) {
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
}

void PairDFT::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    utils::sfread(FLERR, &cut_global, sizeof(double), 1, fp, nullptr, error);
  }
  MPI_Bcast(&cut_global, 1, MPI_DOUBLE, 0, world);
}

double PairDFT::single(int, int, int, int, double, double, double, double &fforce)
{
  error->all(FLERR, "Single() not implemented for pair style dft");
  return 0.0;
}

void *PairDFT::extract(const char *str, int &dim)
{
  dim = 0;
  
  if (strcmp(str, "dft_energy") == 0) return (void *) &total_energy;
  if (strcmp(str, "xc_energy") == 0) return (void *) &xc_energy;
  if (strcmp(str, "n_electrons") == 0) return (void *) &n_electrons;
  if (strcmp(str, "n_basis") == 0) return (void *) &n_basis;
  
  return nullptr;
}
