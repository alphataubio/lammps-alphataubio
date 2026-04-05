// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org
   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_uf3.h"
#include "uf3_potential.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <cstdio>

using namespace LAMMPS_NS;

ComputeUF3::ComputeUF3(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), list(nullptr), array_local(nullptr),
    c_pe(nullptr), c_virial(nullptr)
{
  array_flag = 1;
  extarray = 0;
  virial_flag = 0;

  if (narg < 3) error->all(FLERR,"Illegal compute uf3 command");
  const int nbody = utils::inumeric(FLERR, arg[3], true, lmp);
  if (nbody == 2) pot_3b = false;
  else if (nbody == 3) pot_3b = true;
  else error->all(FLERR, "compute uf3 not (yet) implemented for {}-body terms", nbody);

  const int np1 = atom->ntypes + 1;
  memory->create(setflag, np1, np1, "uf3:setflag");
  memory->create(cutsq, np1, np1, "uf3:cutsq");
  if (pot_3b) memory->create(neighshort, maxshort, "uf3:neighshort");

  std::vector<std::string> elements_(1); // blank [0] to use ntypes+1;
  for(int i=5; i<narg ; i++) elements_.push_back(arg[i]);
  uf3_potential = new UF3Potential(lmp, arg[4], cutsq, setflag, elements_, pot_3b);

  if (virial_flag) size_array_rows = 1 + 3*(atom->natoms) + 6;
  else size_array_rows = 1 + 3*(atom->natoms);
  size_array_cols = ncoeff + 1;
  lastcol = size_array_cols-1;

}

/* ---------------------------------------------------------------------- */

ComputeUF3::~ComputeUF3()
{
  memory->destroy(setflag);
  memory->destroy(cutsq);
  if (pot_3b) memory->destroy(neighshort);
  memory->destroy(array_local);
  memory->destroy(array);
  if( virial_flag && modify->find_compute(id_virial) != -1 ) modify->delete_compute(id_virial);
}

/* ---------------------------------------------------------------------- */

void ComputeUF3::init()
{
  if (force->pair == nullptr)
    error->all(FLERR,"Compute uf3 requires a pair style be defined");

  // need an occasional full neighbor list
  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);

  if (modify->get_compute_by_style("uf3").size() > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute uf3");

  // allocate memory for global array
  memory->create(array_local,size_array_rows,size_array_cols, "uf3:array_local");
  memory->create(array,size_array_rows,size_array_cols, "uf3:array");

  // find compute for reference energy
  c_pe = modify->get_compute_by_id("thermo_pe");
  if (!c_pe) error->all(FLERR,"Compute thermo_pe does not exist.");

  // add compute for reference virial tensor
  if (virial_flag) {
    id_virial = id + std::string("_press");
    c_virial = modify->add_compute(id_virial + " all pressure NULL virial");
  }
}

/* ---------------------------------------------------------------------- */

void ComputeUF3::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeUF3::compute_array()
{
  int ntotal = atom->nlocal + atom->nghost;
  invoked_array = update->ntimestep;

  // clear global array
  for (int i = 0; i < size_array_rows; i++) {
    for (int j = 0; j < size_array_cols; j++) array_local[i][j] = 0.0;
  }

  // invoke full neighbor list (will copy or build if necessary)
  neighbor->build_one(list);

  const int inum = list->inum;
  const int* const ilist = list->ilist;
  const int* const numneigh = list->numneigh;
  int** const firstneigh = list->firstneigh;
  int * const type = atom->type;
  double **x = atom->x;

  // compute uf3 derivatives for each atom in group
  // use full neighbor list to count atoms less than cutoff

  const int* const mask = atom->mask;
  const int ntypes = atom->ntypes;

  for (int ii = 0; ii < inum; ii++) {
    int irow = 0;
    const int i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;

    const int itype = type[i];
    const int* const jlist = firstneigh[i];
    const int jnum = numneigh[i];
    const int row_offset_i = 1 + 3*(atom->tag[i]-1);
    const int type_offset = type_offsets.at(itype);



    for (int jj = 0; jj < jnum; jj++) {
      const int j = jlist[jj];
      const int jtype = type[j];
      const int row_offset_j = 1 + 3*(atom->tag[j]-1);

      const double delx = x[i][0] - x[j][0];
      const double dely = x[i][1] - x[j][1];
      const double delz = x[i][2] - x[j][2];
      const double rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutsq[itype][jtype]) continue;

      const double rij = sqrt(rsq);
      const int start_idx = uf3_potential->get_starting_index_2b(itype, jtype, rij);

      // ENERGY
      const double rth = rsq * rij;
      double **cached_constants_2b = uf3_potential->cached_constants_2b[itype][jtype];
      double evdwl =        cached_constants_2b[start_idx    ][0];
        evdwl += rij * cached_constants_2b[start_idx    ][1];
        evdwl += rsq * cached_constants_2b[start_idx    ][2];
        evdwl += rth * cached_constants_2b[start_idx    ][3];
        evdwl +=       cached_constants_2b[start_idx - 1][4];
        evdwl += rij * cached_constants_2b[start_idx - 1][5];
        evdwl += rsq * cached_constants_2b[start_idx - 1][6];
        evdwl += rth * cached_constants_2b[start_idx - 1][7];
        evdwl +=       cached_constants_2b[start_idx - 2][8];
        evdwl += rij * cached_constants_2b[start_idx - 2][9];
        evdwl += rsq * cached_constants_2b[start_idx - 2][10];
        evdwl += rth * cached_constants_2b[start_idx - 2][11];
        evdwl +=       cached_constants_2b[start_idx - 3][12];
        evdwl += rij * cached_constants_2b[start_idx - 3][13];
        evdwl += rsq * cached_constants_2b[start_idx - 3][14];
        evdwl += rth * cached_constants_2b[start_idx - 3][15];


      // FORCES
      double **cached_constants_2b_deri = uf3_potential->cached_constants_2b_deri[itype][jtype];
      double force_2b = cached_constants_2b_deri[start_idx - 1][0];
      force_2b += rij * cached_constants_2b_deri[start_idx - 1][1];
      force_2b += rsq * cached_constants_2b_deri[start_idx - 1][2];
      force_2b +=       cached_constants_2b_deri[start_idx - 2][3];
      force_2b += rij * cached_constants_2b_deri[start_idx - 2][4];
      force_2b += rsq * cached_constants_2b_deri[start_idx - 2][5];
      force_2b +=       cached_constants_2b_deri[start_idx - 3][6];
      force_2b += rij * cached_constants_2b_deri[start_idx - 3][7];
      force_2b += rsq * cached_constants_2b_deri[start_idx - 3][8];

        const double fpair = -1 * force_2b / rij;
        const double fx = delx * fpair;
        const double fy = dely * fpair;
        const double fz = delz * fpair;













/*


        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;
        f[j][0] -= fx;
        f[j][1] -= fy;
        f[j][2] -= fz;



      // Pseudocode for inside the i-j neighbor loop in compute_uf3.cpp
      double b_val[4]; // Store b_m(rij) here
      double b_der[4]; // Store db_m/drij here

      // ... [Evaluate the 4 active unweighted B-splines based on rij] ...

      for (int local_m = 0; local_m < 4; local_m++) {
        const int func_ind = start_idx - 3 + local_m; // Map to global basis index
        const int col = type_offset + func_ind;

        // 1. Energy Feature
        array_local[0][col] += 0.5 * b_val[local_m]; // 0.5 to avoid double counting if full neighbor list

        // 2. Force Features
        const double force_factor = -b_der[local_m] / rij;
        const double fx_feature = delx * force_factor;
        const double fy_feature = dely * force_factor;
        const double fz_feature = delz * force_factor;

        double **cached_constants_2b_deri = uf3_potential->cached_constants_2b_deri[itype][jtype];
        double force_2b = cached_constants_2b_deri[knot_start_index - 1][0];
        force_2b += rij * cached_constants_2b_deri[knot_start_index - 1][1];
        force_2b += rsq * cached_constants_2b_deri[knot_start_index - 1][2];
        force_2b +=       cached_constants_2b_deri[knot_start_index - 2][3];
        force_2b += rij * cached_constants_2b_deri[knot_start_index - 2][4];
        force_2b += rsq * cached_constants_2b_deri[knot_start_index - 2][5];
        force_2b +=       cached_constants_2b_deri[knot_start_index - 3][6];
        force_2b += rij * cached_constants_2b_deri[knot_start_index - 3][7];
        force_2b += rsq * cached_constants_2b_deri[knot_start_index - 3][8];

        const double fpair = -1 * force_2b / rij;
        const double fx = delx * fpair;
        const double fy = dely * fpair;
        const double fz = delz * fpair;
        
        array_local[row_offset_i    ][col] += fx_feature;
        array_local[row_offset_i + 1][col] += fy_feature;
        array_local[row_offset_i + 2][col] += fz_feature;
    
        array_local[row_offset_j    ][col] -= fx_feature;
        array_local[row_offset_j + 1][col] -= fy_feature;
        array_local[row_offset_j + 2][col] -= fz_feature;



        // 3. Virial Features (if requested)
        if (virial_flag) {
          array_local[size_array_rows-6][col] += delx * fx_feature; // W_xx
          array_local[size_array_rows-5][col] += dely * fy_feature; // W_yy
          array_local[size_array_rows-4][col] += delz * fz_feature; // W_zz
          array_local[size_array_rows-3][col] += delz * fy_feature; // W_zy
          array_local[size_array_rows-2][col] += delz * fx_feature; // W_zx
          array_local[size_array_rows-1][col] += dely * fx_feature; // W_yx
        }

        */

    } // loop over jj inside
  } // for ii loop

  // accumulate forces to global array
  for (int i = 0; i < atom->nlocal; i++) {
    int iglobal = atom->tag[i];
    int irow = 3*(iglobal-1)+1;
    array_local[irow++][lastcol] = atom->f[i][0];
    array_local[irow++][lastcol] = atom->f[i][1];
    array_local[irow][lastcol] = atom->f[i][2];
  }

  // sum up over all processes
  MPI_Allreduce(&array_local[0][0],&array[0][0],size_array_rows*size_array_cols,MPI_DOUBLE,MPI_SUM,world);

  // assign energy to last column
  array[0][lastcol] = c_pe->compute_scalar();

  // assign virial stress to last column
  // switch to Voigt notation
  if (virial_flag) {
    c_virial->compute_vector();
    int irow = 1 + 3*(atom->natoms);
    array[irow++][lastcol] = c_virial->vector[0];
    array[irow++][lastcol] = c_virial->vector[1];
    array[irow++][lastcol] = c_virial->vector[2];
    array[irow++][lastcol] = c_virial->vector[5];
    array[irow++][lastcol] = c_virial->vector[4];
    array[irow++][lastcol] = c_virial->vector[3];
  }
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double ComputeUF3::memory_usage()
{
  double bytes = (double)size_array_rows*size_array_cols*sizeof(double)*2; // uf3 and uf3_all
  return bytes;
}


