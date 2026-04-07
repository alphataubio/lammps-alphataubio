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

static constexpr int leading_trim = 3;
// allow column index (n_WW - 1 - trailing_trim), e.g. four active bases up to WW8 when size is 11
static constexpr int trailing_trim = 2;

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

  size_array_cols = 1; // last column is regression b vector
  const int ntypes = atom->ntypes;
  for (int i = 1; i <= ntypes; i++) {
    for (int j = 1; j <= ntypes; j++) {
      const double cut_2b_ij = uf3_potential->cut_2b[i][j];
      cutsq[i][j] = cut_2b_ij * cut_2b_ij;
      size_array_cols += uf3_potential->n2b_coeff_array_size[i][j];
      //fprintf(stderr, "*** n2b_coeff_array_size[%i][%i] %i\n", i, j, uf3_potential->n2b_coeff_array_size[i][j]);
      if (pot_3b) {
        for (int k = 1; k <= ntypes; k++) {
          const int map_to = uf3_potential->map_3b[i][j][k];
          size_array_cols += uf3_potential->n3b_coeff_array_size[map_to][0];
          size_array_cols += uf3_potential->n3b_coeff_array_size[map_to][1];
          size_array_cols += uf3_potential->n3b_coeff_array_size[map_to][2];

          //fprintf(stderr, "*** map_3b[%i][%i][%i] %i n3b_coeff_array_size[map_to] %i %i %i\n", i, j, k, map_to, uf3_potential->n3b_coeff_array_size[map_to][0], uf3_potential->n3b_coeff_array_size[map_to][1], uf3_potential->n3b_coeff_array_size[map_to][2]);
        }
      }
    }
  }
  //fprintf(stderr, "*** size_array_cols %i\n", size_array_cols);
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
  const int* const type = atom->type;
  double** const x = atom->x;
  const int* const mask = atom->mask;

  // compute uf3 derivatives for each atom in group
  // use full neighbor list to count atoms less than cutoff


  for (int ii = 0; ii < inum; ii++) {

    const int i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;

    const int itype = type[i];
    const int* const jlist = firstneigh[i];
    const int jnum = numneigh[i];
    const int row_offset_i = 1 + 3*(atom->tag[i]-1);

    for (int jj = 0; jj < jnum; jj++) {
      const int j = jlist[jj];
      const int jtype = type[j];
      const int row_offset_j = 1 + 3*(atom->tag[j]-1);

      const double delx = x[i][0] - x[j][0];
      const double dely = x[i][1] - x[j][1];
      const double delz = x[i][2] - x[j][2];
      const double rsq = delx*delx + dely*dely + delz*delz;
      if (rsq >= cutsq[itype][jtype]) continue;

      const double rij = sqrt(rsq);
      const double rth = rsq * rij;
      const int start_idx = uf3_potential->get_starting_index_2b(itype, jtype, rij);
      //fprintf(stderr, "*** rsq %f start_idx %i x[%i] %f %f %f x[%i] %f %f %f\n", rsq, start_idx, i, x[i][0], x[i][1], x[i][2], j, x[j][0], x[j][1], x[j][2]);

      // ENERGY & FORCE FEATURE EXTRACTION
      double **constants_2b = &(uf3_potential->cached_constants_2b[itype][jtype][start_idx-3]);
      double **constants_2b_deri = &(uf3_potential->cached_constants_2b_deri[itype][jtype][start_idx-3]);

      // Set your column offset per interaction type (e.g., A-A vs A-B)
      int type_offset = 0; 

      // Extract the 4 active basis functions directly from your pre-computed matrices
      for (int m = 0; m < 4; m++) {
        // map knot segment (start_idx-3)+m to WW column matching external tables (1-based segment index)
        const int col_offset = type_offset + m + start_idx - 3;
        if ( col_offset < leading_trim || col_offset > size_array_cols-trailing_trim-3 ) continue;
        // 1. Energy Feature (Cubic: 4 coefficients)
        const int n = (3 - m) * 4;
        array_local[0][col_offset] +=         constants_2b[m][n]
                                      + rij * constants_2b[m][n+1]
                                      + rsq * constants_2b[m][n+2]
                                      + rth * constants_2b[m][n+3];
        // 2. Force Feature (Derivative of the cubic: C1 + 2*C2*r + 3*C3*r^2)
        const double b_der = constants_2b[m][n+1] + 2.0 * rij * constants_2b[m][n+2] + 3.0 * rsq * constants_2b[m][n+3];
        const double force_factor = -b_der / rij;
        const double fx_feature = delx * force_factor;
        const double fy_feature = dely * force_factor;
        const double fz_feature = delz * force_factor;
        array_local[row_offset_i  ][col_offset] += fx_feature;
        array_local[row_offset_i+1][col_offset] += fy_feature;
        array_local[row_offset_i+2][col_offset] += fz_feature;
        array_local[row_offset_j  ][col_offset] -= fx_feature;
        array_local[row_offset_j+1][col_offset] -= fy_feature;
        array_local[row_offset_j+2][col_offset] -= fz_feature;

        // Virials (If requested)
        if (virial_flag) {
          array_local[size_array_rows-6][col_offset] += delx * fx_feature; // W_xx
          array_local[size_array_rows-5][col_offset] += dely * fy_feature; // W_yy
          array_local[size_array_rows-4][col_offset] += delz * fz_feature; // W_zz
          array_local[size_array_rows-3][col_offset] += delz * fy_feature; // W_zy
          array_local[size_array_rows-2][col_offset] += delz * fx_feature; // W_zx
          array_local[size_array_rows-1][col_offset] += dely * fx_feature; // W_yx
        }
      }
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


