/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_libxc.h"
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
#include <xc.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLibXC::PairLibXC(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 1;
  restartinfo = 1;
  one_coeff = 0;
  
  xc_func_x = nullptr;
  xc_func_c = nullptr;
  xc_func_xc = nullptr;
  
  xc_functional_x = -1;
  xc_functional_c = -1;
  xc_functional_xc = -1;
  
  use_combined_xc = false;
  
  ngrid = 100;
  grid_spacing = 0.01;
  
  rho0 = nullptr;
  decay_length = nullptr;
  atomic_volume = nullptr;
}

/* ---------------------------------------------------------------------- */

PairLibXC::~PairLibXC()
{
  if (xc_func_x) {
    xc_func_end(xc_func_x);
    delete xc_func_x;
  }
  if (xc_func_c) {
    xc_func_end(xc_func_c);
    delete xc_func_c;
  }
  if (xc_func_xc) {
    xc_func_end(xc_func_xc);
    delete xc_func_xc;
  }
  
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(offset);
    memory->destroy(rho0);
    memory->destroy(decay_length);
    memory->destroy(atomic_volume);
  }
}

/* ---------------------------------------------------------------------- */

void PairLibXC::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;
  double rsq, r, rho_i, rho_j, rho_avg, drho_dr;
  double exc_energy, vxc, d_vxc_drho;
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  // Allocate arrays for density computation
  double *rho_total = new double[nlocal];
  double *vxc_local = new double[nlocal];
  
  // Initialize density arrays
  for (i = 0; i < nlocal; i++) {
    rho_total[i] = 0.0;
    vxc_local[i] = 0.0;
  }
  
  // First pass: compute electron density at each atom
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    
    // Add self-density
    rho_total[i] = rho0[itype][itype];
    
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];
      
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      
      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);
        // Add neighbor contribution to density
        rho_total[i] += electron_density(r, jtype);
      }
    }
  }
  
  // Compute XC energy and potential for each atom
  for (i = 0; i < nlocal; i++) {
    if (rho_total[i] > 0.0) {
      double rho = rho_total[i];
      double sigma = 0.0; // For spin-unpolarized calculation
      
      if (use_combined_xc) {
        // Use combined XC functional
        xc_lda_exc_vxc(xc_func_xc, 1, &rho, &exc_energy, &vxc_local[i]);
      } else {
        // Use separate X and C functionals
        double ex, ec, vx, vc;
        if (xc_func_x) {
          xc_lda_exc_vxc(xc_func_x, 1, &rho, &ex, &vx);
        } else {
          ex = vx = 0.0;
        }
        if (xc_func_c) {
          xc_lda_exc_vxc(xc_func_c, 1, &rho, &ec, &vc);
        } else {
          ec = vc = 0.0;
        }
        exc_energy = ex + ec;
        vxc_local[i] = vx + vc;
      }
    }
  }
  
  // Second pass: compute forces and energy
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);
        
        // Compute force from XC functional derivative
        drho_dr = electron_density_derivative(r, jtype);
        
        // Force contribution from XC potential
        fpair = -(vxc_local[i] + vxc_local[j]) * drho_dr / r;
        
        // Add pairwise repulsion or other terms if needed
        // This is a simplified implementation
        
        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx * fpair;
          f[j][1] -= dely * fpair;
          f[j][2] -= delz * fpair;
        }

        if (eflag) {
          // Simplified energy calculation
          evdwl = rho_total[i] * vxc_local[i] * atomic_volume[itype][jtype];
          evdwl -= offset[itype][jtype];
        }

        if (evflag) ev_tally(i, j, nlocal, newton_pair,
                            evdwl, 0.0, fpair, delx, dely, delz);
      }
    }
  }
  
  delete[] rho_total;
  delete[] vxc_local;

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   electron density function
------------------------------------------------------------------------- */

double PairLibXC::electron_density(double r, int itype)
{
  // Simple exponential decay model for electron density
  double decay = decay_length[itype][itype];
  double rho = rho0[itype][itype] * exp(-r / decay);
  return rho;
}

/* ----------------------------------------------------------------------
   derivative of electron density function
------------------------------------------------------------------------- */

double PairLibXC::electron_density_derivative(double r, int itype)
{
  // Derivative of exponential decay model
  double decay = decay_length[itype][itype];
  double drho = -rho0[itype][itype] * exp(-r / decay) / decay;
  return drho;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLibXC::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
  memory->create(cut, n + 1, n + 1, "pair:cut");
  memory->create(offset, n + 1, n + 1, "pair:offset");
  memory->create(rho0, n + 1, n + 1, "pair:rho0");
  memory->create(decay_length, n + 1, n + 1, "pair:decay_length");
  memory->create(atomic_volume, n + 1, n + 1, "pair:atomic_volume");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLibXC::settings(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR, "Illegal pair_style libxc command");
  
  // Parse functional specification
  // Format: pair_style libxc <functional> <cutoff> [options]
  // Example: pair_style libxc LDA_X+LDA_C_VWN 10.0
  // Example: pair_style libxc PBE 12.0
  
  char *functional_str = arg[0];
  cut_global = utils::numeric(FLERR, arg[1], false, lmp);
  
  // Parse functional string
  if (strchr(functional_str, '+') != nullptr) {
    // Separate X and C functionals
    use_combined_xc = false;
    char *x_func = strtok(functional_str, "+");
    char *c_func = strtok(nullptr, "+");
    
    if (x_func) {
      parse_functional_name(x_func, xc_functional_x);
      xc_func_x = new xc_func_type;
      if (xc_func_init(xc_func_x, xc_functional_x, XC_UNPOLARIZED) != 0) {
        error->all(FLERR, "Failed to initialize exchange functional");
      }
    }
    
    if (c_func) {
      parse_functional_name(c_func, xc_functional_c);
      xc_func_c = new xc_func_type;
      if (xc_func_init(xc_func_c, xc_functional_c, XC_UNPOLARIZED) != 0) {
        error->all(FLERR, "Failed to initialize correlation functional");
      }
    }
  } else {
    // Combined XC functional
    use_combined_xc = true;
    parse_functional_name(functional_str, xc_functional_xc);
    xc_func_xc = new xc_func_type;
    if (xc_func_init(xc_func_xc, xc_functional_xc, XC_UNPOLARIZED) != 0) {
      error->all(FLERR, "Failed to initialize XC functional");
    }
  }
  
  // Parse additional options if provided
  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "ngrid") == 0) {
      if (iarg + 1 >= narg) error->all(FLERR, "Missing ngrid value");
      ngrid = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "spacing") == 0) {
      if (iarg + 1 >= narg) error->all(FLERR, "Missing spacing value");
      grid_spacing = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else {
      error->all(FLERR, "Unknown pair_style libxc option");
    }
  }

  // reset cutoffs that have been explicitly set
  if (allocated) {
    int i, j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   parse functional name to LibXC ID
------------------------------------------------------------------------- */

void PairLibXC::parse_functional_name(const char *name, int &func_id)
{
  // Map common functional names to LibXC IDs
  // This is a simplified mapping - a full implementation would
  // use xc_functional_get_number() from LibXC
  
  if (strcmp(name, "LDA_X") == 0) func_id = XC_LDA_X;
  else if (strcmp(name, "LDA_C_VWN") == 0) func_id = XC_LDA_C_VWN;
  else if (strcmp(name, "LDA_C_PW") == 0) func_id = XC_LDA_C_PW;
  else if (strcmp(name, "GGA_X_PBE") == 0) func_id = XC_GGA_X_PBE;
  else if (strcmp(name, "GGA_C_PBE") == 0) func_id = XC_GGA_C_PBE;
  else if (strcmp(name, "PBE") == 0) func_id = XC_GGA_XC_PBE;
  else if (strcmp(name, "BLYP") == 0) func_id = XC_GGA_XC_B88_LYP;
  else if (strcmp(name, "B3LYP") == 0) func_id = XC_HYB_GGA_XC_B3LYP;
  else {
    // Try to get functional number from LibXC
    func_id = xc_functional_get_number(name);
    if (func_id == -1) {
      error->all(FLERR, "Unknown LibXC functional name");
    }
  }
}

/* ----------------------------------------------------------------------
   set coefficients for one or more type pairs
------------------------------------------------------------------------- */

void PairLibXC::coeff(int narg, char **arg)
{
  if (narg < 6 || narg > 7) 
    error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double rho0_one = utils::numeric(FLERR, arg[2], false, lmp);
  double decay_one = utils::numeric(FLERR, arg[3], false, lmp);
  double volume_one = utils::numeric(FLERR, arg[4], false, lmp);
  double cut_one = utils::numeric(FLERR, arg[5], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      rho0[i][j] = rho0_one;
      decay_length[i][j] = decay_one;
      atomic_volume[i][j] = volume_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLibXC::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR, "Pair style libxc requires atom IDs");
  
  // Request standard neighbor list
  neighbor->add_request(this, NeighConst::REQ_DEFAULT);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLibXC::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    rho0[i][j] = mix_energy(rho0[i][i], rho0[j][j],
                            decay_length[i][i], decay_length[j][j]);
    decay_length[i][j] = mix_distance(decay_length[i][i], decay_length[j][j]);
    atomic_volume[i][j] = mix_energy(atomic_volume[i][i], atomic_volume[j][j],
                                     decay_length[i][i], decay_length[j][j]);
    cut[i][j] = mix_distance(cut[i][i], cut[j][j]);
  }

  rho0[j][i] = rho0[i][j];
  decay_length[j][i] = decay_length[i][j];
  atomic_volume[j][i] = atomic_volume[i][j];

  // Compute offset energy at cutoff
  if (offset_flag && (cut[i][j] > 0.0)) {
    double rho_cut = electron_density(cut[i][j], j);
    double vxc_cut = 0.0;
    
    if (use_combined_xc) {
      double exc;
      xc_lda_exc_vxc(xc_func_xc, 1, &rho_cut, &exc, &vxc_cut);
    } else {
      double ex, ec, vx, vc;
      if (xc_func_x) {
        xc_lda_exc_vxc(xc_func_x, 1, &rho_cut, &ex, &vx);
      } else {
        vx = 0.0;
      }
      if (xc_func_c) {
        xc_lda_exc_vxc(xc_func_c, 1, &rho_cut, &ec, &vc);
      } else {
        vc = 0.0;
      }
      vxc_cut = vx + vc;
    }
    
    offset[i][j] = rho_cut * vxc_cut * atomic_volume[i][j];
  } else {
    offset[i][j] = 0.0;
  }
  
  offset[j][i] = offset[i][j];
  
  return cut[i][j];
}

/* ----------------------------------------------------------------------
   write pair coeffs to restart file
------------------------------------------------------------------------- */

void PairLibXC::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) {
        fwrite(&rho0[i][j], sizeof(double), 1, fp);
        fwrite(&decay_length[i][j], sizeof(double), 1, fp);
        fwrite(&atomic_volume[i][j], sizeof(double), 1, fp);
        fwrite(&cut[i][j], sizeof(double), 1, fp);
      }
    }
}

/* ----------------------------------------------------------------------
   read pair coeffs from restart file
------------------------------------------------------------------------- */

void PairLibXC::read_restart(FILE *fp)
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
        if (me == 0) {
          utils::sfread(FLERR, &rho0[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &decay_length[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &atomic_volume[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cut[i][j], sizeof(double), 1, fp, nullptr, error);
        }
        MPI_Bcast(&rho0[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&decay_length[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&atomic_volume[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cut[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
}

/* ----------------------------------------------------------------------
   write pair style settings to restart file
------------------------------------------------------------------------- */

void PairLibXC::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global, sizeof(double), 1, fp);
  fwrite(&xc_functional_x, sizeof(int), 1, fp);
  fwrite(&xc_functional_c, sizeof(int), 1, fp);
  fwrite(&xc_functional_xc, sizeof(int), 1, fp);
  fwrite(&use_combined_xc, sizeof(bool), 1, fp);
  fwrite(&ngrid, sizeof(int), 1, fp);
  fwrite(&grid_spacing, sizeof(double), 1, fp);
  fwrite(&offset_flag, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
   read pair style settings from restart file
------------------------------------------------------------------------- */

void PairLibXC::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    utils::sfread(FLERR, &cut_global, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &xc_functional_x, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &xc_functional_c, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &xc_functional_xc, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &use_combined_xc, sizeof(bool), 1, fp, nullptr, error);
    utils::sfread(FLERR, &ngrid, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &grid_spacing, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &offset_flag, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&cut_global, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&xc_functional_x, 1, MPI_INT, 0, world);
  MPI_Bcast(&xc_functional_c, 1, MPI_INT, 0, world);
  MPI_Bcast(&xc_functional_xc, 1, MPI_INT, 0, world);
  MPI_Bcast(&use_combined_xc, 1, MPI_C_BOOL, 0, world);
  MPI_Bcast(&ngrid, 1, MPI_INT, 0, world);
  MPI_Bcast(&grid_spacing, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
  
  // Re-initialize LibXC functionals
  if (xc_functional_x >= 0) {
    xc_func_x = new xc_func_type;
    xc_func_init(xc_func_x, xc_functional_x, XC_UNPOLARIZED);
  }
  if (xc_functional_c >= 0) {
    xc_func_c = new xc_func_type;
    xc_func_init(xc_func_c, xc_functional_c, XC_UNPOLARIZED);
  }
  if (xc_functional_xc >= 0) {
    xc_func_xc = new xc_func_type;
    xc_func_init(xc_func_xc, xc_functional_xc, XC_UNPOLARIZED);
  }
}

/* ----------------------------------------------------------------------
   write pair data to text format
------------------------------------------------------------------------- */

void PairLibXC::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp, "%d %g %g %g\n", i, rho0[i][i], decay_length[i][i], atomic_volume[i][i]);
}

/* ----------------------------------------------------------------------
   write all pair data to text format
------------------------------------------------------------------------- */

void PairLibXC::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp, "%d %d %g %g %g %g\n", i, j,
              rho0[i][j], decay_length[i][j], atomic_volume[i][j], cut[i][j]);
}

/* ----------------------------------------------------------------------
   single() function for energy and force calculations
------------------------------------------------------------------------- */

double PairLibXC::single(int i, int j, int itype, int jtype,
                         double rsq, double factor_coul, double factor_lj,
                         double &fforce)
{
  double r = sqrt(rsq);
  double rho_ij = electron_density(r, jtype);
  double drho_dr = electron_density_derivative(r, jtype);
  
  // Simplified single interaction calculation
  // In a full implementation, this would need to consider the
  // local density environment
  
  double vxc = 0.0;
  if (use_combined_xc) {
    double exc;
    xc_lda_exc_vxc(xc_func_xc, 1, &rho_ij, &exc, &vxc);
  } else {
    double ex, ec, vx, vc;
    if (xc_func_x) {
      xc_lda_exc_vxc(xc_func_x, 1, &rho_ij, &ex, &vx);
    } else {
      vx = 0.0;
    }
    if (xc_func_c) {
      xc_lda_exc_vxc(xc_func_c, 1, &rho_ij, &ec, &vc);
    } else {
      vc = 0.0;
    }
    vxc = vx + vc;
  }
  
  fforce = -vxc * drho_dr / r * factor_lj;
  
  double phi = rho_ij * vxc * atomic_volume[itype][jtype];
  return factor_lj * phi;
}
