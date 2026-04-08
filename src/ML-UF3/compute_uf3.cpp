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

#include <algorithm>
#include <cstdio>

using namespace LAMMPS_NS;

static constexpr int leading_trim = 3;
static constexpr int trailing_trim = 2;

ComputeUF3::ComputeUF3(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), cutsq(nullptr), setflag(nullptr),
    neighshort(nullptr), type_offset_2b(nullptr), type_offset_3b(nullptr),
    list(nullptr), array_local(nullptr), c_pe(nullptr), c_virial(nullptr)
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
  if (pot_3b) {
    maxshort = 20;
    memory->create(neighshort, maxshort, "uf3:neighshort");
    memory->create(cutsq, np1, np1, "uf3:cutsq");
  }

  std::vector<std::string> elements_(1); // blank [0] to use ntypes+1;
  for(int i=5; i<narg ; i++) elements_.push_back(arg[i]);
  uf3_potential = new UF3Potential(lmp, arg[4], cutsq, setflag, elements_, pot_3b);
  memory->create(type_offset_2b, np1, np1, "uf3:type_offset_2b");
  memory->create(type_offset_3b, uf3_potential->tot_interaction_count_3b, "uf3:type_offset_3b");

  if (virial_flag) size_array_rows = 1 + 3*(atom->natoms) + 6;
  else size_array_rows = 1 + 3*(atom->natoms);

  size_array_cols = 0;
  const int ntypes = atom->ntypes;
  for (int i = 1; i <= ntypes; i++) {
    for (int j = i; j <= ntypes; j++) {
      const double cut_2b_ij = uf3_potential->cut_2b[i][j];
      cutsq[i][j] = cut_2b_ij * cut_2b_ij;
      type_offset_2b[i][j] = type_offset_2b[j][i] = size_array_cols;
      size_array_cols += uf3_potential->n2b_coeff_array_size[i][j];
      //fprintf(stderr, "*** n2b_coeff_array_size[%i][%i] %i\n", i, j, uf3_potential->n2b_coeff_array_size[i][j]);
    }
  }

  if (pot_3b) {
    for (int i = 1; i <= ntypes; i++) {
      for (int j = i; j <= ntypes; j++) {
        for (int k = j; k <= ntypes; k++) {
          const int map_to = uf3_potential->map_3b[i][j][k];
          type_offset_3b[map_to] = size_array_cols;
          auto n3b_coeff_size = uf3_potential->n3b_coeff_array_size[map_to];
          size_array_cols += n3b_coeff_size[0] * n3b_coeff_size[1] * n3b_coeff_size[2];

          //fprintf(stderr, "*** map_3b[%i][%i][%i] %i n3b_coeff_array_size[map_to] %i %i %i\n", i, j, k, map_to, uf3_potential->n3b_coeff_array_size[map_to][0], uf3_potential->n3b_coeff_array_size[map_to][1], uf3_potential->n3b_coeff_array_size[map_to][2]);
        }
      }
    }
  }

  size_array_cols++; // last column is regression b vector
  fprintf(stderr, "*** size_array_cols %i\n", size_array_cols);
  lastcol = size_array_cols-1;

}

/* ---------------------------------------------------------------------- */

ComputeUF3::~ComputeUF3()
{
  memory->destroy(setflag);
  memory->destroy(cutsq);
  memory->destroy(type_offset_2b);
  if (pot_3b) {
    memory->destroy(neighshort);
    memory->destroy(type_offset_3b);
  }
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

    const double xi0 = x[i][0];
    const double xi1 = x[i][1];
    const double xi2 = x[i][2];
    const int itype = type[i];
    const int* const jlist = firstneigh[i];
    int numshort = 0;
    const int row_offset_i = 1 + 3*(atom->tag[i]-1);

    for (int jj = 0; jj < numneigh[i]; jj++) {
      const int j = jlist[jj];
      const int jtype = type[j];
      const int row_offset_j = 1 + 3*(atom->tag[j]-1);

      const double delx = xi0 - x[j][0];
      const double dely = xi1 - x[j][1];
      const double delz = xi2 - x[j][2];
      const double rsq = delx*delx + dely*dely + delz*delz;
      if (rsq >= cutsq[itype][jtype]) continue;

      const double rij = sqrt(rsq);
      const double rth = rsq * rij;
      const int start_idx = uf3_potential->get_starting_index_2b(itype, jtype, rij);
      //fprintf(stderr, "*** rsq %f start_idx %i x[%i] %f %f %f x[%i] %f %f %f\n", rsq, start_idx, i, x[i][0], x[i][1], x[i][2], j, x[j][0], x[j][1], x[j][2]);

      if (pot_3b) {
        if (rij <= uf3_potential->cut_3b_list[itype][jtype]) {
          neighshort[numshort] = j;
          if (numshort >= maxshort - 1) {
            maxshort += maxshort / 2;
            memory->grow(neighshort, maxshort, "pair:neighshort");
          }
          numshort = numshort + 1;
        }
      }

      // ENERGY & FORCE FEATURE EXTRACTION
      double **cc_2b = &(uf3_potential->cached_constants_2b[itype][jtype][start_idx-3]);
      double **cc_2b_deri = &(uf3_potential->cached_constants_2b_deri[itype][jtype][start_idx-3]);

      // Extract the 4 active basis functions directly from your pre-computed matrices
      for (int m = 0; m < 4; m++) {
        const int col_offset = type_offset_2b[itype][jtype] + m + start_idx - 3;
        if ( col_offset < leading_trim || col_offset > size_array_cols-trailing_trim-3 ) continue; //FIXME multielement
        // 1. Energy Feature (Cubic: 4 coefficients)
        const int n = (3 - m) * 4;
        array_local[0][col_offset] += cc_2b[m][n] + rij * cc_2b[m][n+1] + rsq * cc_2b[m][n+2] + rth * cc_2b[m][n+3];
        // 2. Force Feature (Derivative of the cubic: C1 + 2*C2*r + 3*C3*r^2)
        const double b_der = cc_2b[m][n+1] + 2.0 * rij * cc_2b[m][n+2] + 3.0 * rsq * cc_2b[m][n+3];
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

    // 3-body interaction
    // jth atom

    for (int jj = 0; jj < numshort - 1; jj++) {

      double del_rji[3], del_rki[3], del_rkj[3];

      const int j = neighshort[jj];
      const int jtype = type[j];
      const int row_offset_j = 1 + 3*(atom->tag[j]-1);

      del_rji[0] = x[j][0] - xi0;
      del_rji[1] = x[j][1] - xi1;
      del_rji[2] = x[j][2] - xi2;
      const double rij_sq = (del_rji[0] * del_rji[0]) + (del_rji[1] * del_rji[1]) + (del_rji[2] * del_rji[2]);
      const double rij = sqrt(rij_sq);

      // kth atom
      for (int kk = jj + 1; kk < numshort; kk++) {

        const int k = neighshort[kk];
        const int ktype = type[k];
        const int row_offset_k = 1 + 3*(atom->tag[k]-1);

        del_rki[0] = x[k][0] - xi0;
        del_rki[1] = x[k][1] - xi1;
        del_rki[2] = x[k][2] - xi2;
        const double rik_sq = (del_rki[0] * del_rki[0]) + (del_rki[1] * del_rki[1]) + (del_rki[2] * del_rki[2]);
        const double rik = sqrt(rik_sq);
        auto cut_3b_i = uf3_potential->cut_3b[itype];
        auto min_cut_3b_ijk = uf3_potential->min_cut_3b[itype][jtype][ktype];
        if ( rij > cut_3b_i[jtype][ktype] || rij < min_cut_3b_ijk[2] ) continue;
        if ( rik > cut_3b_i[ktype][jtype] || rik < min_cut_3b_ijk[1] ) continue;

        del_rkj[0] = x[k][0] - x[j][0];
        del_rkj[1] = x[k][1] - x[j][1];
        del_rkj[2] = x[k][2] - x[j][2];
        const double rjk_sq =(del_rkj[0] * del_rkj[0]) + (del_rkj[1] * del_rkj[1]) + (del_rkj[2] * del_rkj[2]);
        const double rjk = sqrt(rjk_sq);
        if (rjk < min_cut_3b_ijk[0]) continue;

        const double rij_th = rij * rij_sq;
        const double rik_th = rik * rik_sq;
        const double rjk_th = rjk * rjk_sq;

        const int map_to = uf3_potential->map_3b[itype][jtype][ktype];
        double ***cached_constants_3b = uf3_potential->cached_constants_3b[map_to];
        double ***cached_constants_3b_deri = uf3_potential->cached_constants_3b_deri[map_to];
        const int iknot_ij = uf3_potential->get_starting_index_3b(itype, jtype, ktype, rij, 2) - 3;
        const int iknot_ik = uf3_potential->get_starting_index_3b(itype, jtype, ktype, rik, 1) - 3;
        const int iknot_jk = uf3_potential->get_starting_index_3b(itype, jtype, ktype, rjk, 0) - 3;
        double basis_ij[4], basis_ik[4], basis_jk[4], basis_ij_der[3], basis_ik_der[3], basis_jk_der[3];

        fprintf(stderr, "*** type_offset_3b[%i] %i iknot_ij %i iknot_ik %i iknot_jk %i\n", map_to, type_offset_3b[map_to], iknot_ij, iknot_ik, iknot_jk);



        // -------- basis_ij_der --------
        auto cc_3b_deri_ij = &(cached_constants_3b_deri[0][iknot_ij]);
        basis_ij_der[0] = cc_3b_deri_ij[0][6] + rij * cc_3b_deri_ij[0][7] + rij_sq * cc_3b_deri_ij[0][8];
        basis_ij_der[1] = cc_3b_deri_ij[1][3] + rij * cc_3b_deri_ij[1][4] + rij_sq * cc_3b_deri_ij[1][5];
        basis_ij_der[2] = cc_3b_deri_ij[2][0] + rij * cc_3b_deri_ij[2][1] + rij_sq * cc_3b_deri_ij[2][2];

        // -------- basis_ik_der --------
        auto cc_3b_deri_ik = &(cached_constants_3b_deri[1][iknot_ik]);
        basis_ik_der[0] = cc_3b_deri_ik[0][6] + rik * cc_3b_deri_ik[0][7] + rik_sq * cc_3b_deri_ik[0][8];
        basis_ik_der[1] = cc_3b_deri_ik[1][3] + rik * cc_3b_deri_ik[1][4] + rik_sq * cc_3b_deri_ik[1][5];
        basis_ik_der[2] = cc_3b_deri_ik[2][0] + rik * cc_3b_deri_ik[2][1] + rik_sq * cc_3b_deri_ik[2][2];

        // -------- basis_jk_der --------
        auto cc_3b_deri_jk = &(cached_constants_3b_deri[2][iknot_jk]);
        basis_jk_der[0] = cc_3b_deri_jk[0][6] + rjk * cc_3b_deri_jk[0][7] + rjk_sq * cc_3b_deri_jk[0][8];
        basis_jk_der[1] = cc_3b_deri_jk[1][3] + rjk * cc_3b_deri_jk[1][4] + rjk_sq * cc_3b_deri_jk[1][5];
        basis_jk_der[2] = cc_3b_deri_jk[2][0] + rjk * cc_3b_deri_jk[2][1] + rjk_sq * cc_3b_deri_jk[2][2];

        // Ensure 4th element padding for derivative arrays to unify the loops
        const double d_bij[4] = {basis_ij_der[0], basis_ij_der[1], basis_ij_der[2], 0.0};
        const double d_bik[4] = {basis_ik_der[0], basis_ik_der[1], basis_ik_der[2], 0.0};
        const double d_bjk[4] = {basis_jk_der[0], basis_jk_der[1], basis_jk_der[2], 0.0};

        // Precompute spatial unit vectors
        const double dx_ij = del_rji[0] / rij;
        const double dy_ij = del_rji[1] / rij;
        const double dz_ij = del_rji[2] / rij;

        const double dx_ik = del_rki[0] / rik;
        const double dy_ik = del_rki[1] / rik;
        const double dz_ik = del_rki[2] / rik;

        const double dx_jk = del_rkj[0] / rjk;
        const double dy_jk = del_rkj[1] / rjk;
        const double dz_jk = del_rkj[2] / rjk;

        const int base_offset = type_offset_3b[map_to];
        const int K_m = uf3_potential->n3b_coeff_array_size[map_to][1];
        const int K_n = uf3_potential->n3b_coeff_array_size[map_to][0];

        // Unified Descriptor Accumulation
        for (int l = 0; l < 4; l++) {
          const double b_ij  = basis_ij[l];
          const double db_ij = d_bij[l];
          
          for (int m = 0; m < 4; m++) {
            const double b_ik  = basis_ik[m];
            const double db_ik = d_bik[m];
            
            // Flattened 1D matrix column for this (l, m) tensor slice
            const int col_offset = base_offset + (iknot_ij + l) * K_m * K_n + (iknot_ik + m) * K_n + iknot_jk;

            // Precompute constant base products for the final n-loop
            const double term1_base = db_ij * b_ik; 
            const double term2_base = b_ij * db_ik; 
            const double term3_base = b_ij * b_ik;  

            for (int n = 0; n < 4; n++) {
              const int current_col = col_offset + n;
              
              const double b_jk  = basis_jk[n];
              const double db_jk = d_bjk[n];

              // --- 1. ENERGY DESCRIPTOR ---
              array_local[0][current_col] += term3_base * b_jk;

              // --- 2. FORCE DESCRIPTORS ---
              const double d_ij_part = term1_base * b_jk;
              const double d_ik_part = term2_base * b_jk;
              const double d_jk_part = term3_base * db_jk;

              // Project onto Cartesian axes
              const double fij_x = d_ij_part * dx_ij;
              const double fik_x = d_ik_part * dx_ik;
              const double fjk_x = d_jk_part * dx_jk;

              const double fij_y = d_ij_part * dy_ij;
              const double fik_y = d_ik_part * dy_ik;
              const double fjk_y = d_jk_part * dy_jk;

              const double fij_z = d_ij_part * dz_ij;
              const double fik_z = d_ik_part * dz_ik;
              const double fjk_z = d_jk_part * dz_jk;

              // Force on atom i
              array_local[row_offset_i    ][current_col] += (fij_x + fik_x);
              array_local[row_offset_i + 1][current_col] += (fij_y + fik_y);
              array_local[row_offset_i + 2][current_col] += (fij_z + fik_z);

              // Force on atom j
              array_local[row_offset_j    ][current_col] += (-fij_x + fjk_x);
              array_local[row_offset_j + 1][current_col] += (-fij_y + fjk_y);
              array_local[row_offset_j + 2][current_col] += (-fij_z + fjk_z);

              // Force on atom k
              array_local[row_offset_k    ][current_col] -= (fik_x + fjk_x);
              array_local[row_offset_k + 1][current_col] -= (fik_y + fjk_y);
              array_local[row_offset_k + 2][current_col] -= (fik_z + fjk_z);
            }
          }
        }

















        // -------- basis_ij --------
        auto cc_3b_ij = &(cached_constants_3b[0][iknot_ij]);
        basis_ij[0] = cc_3b_ij[0][12] + rij * cc_3b_ij[0][13] + rij_sq * cc_3b_ij[0][14] + rij_th * cc_3b_ij[0][15];
        basis_ij[1] = cc_3b_ij[1][8]  + rij * cc_3b_ij[1][9]  + rij_sq * cc_3b_ij[1][10] + rij_th * cc_3b_ij[1][11];
        basis_ij[2] = cc_3b_ij[2][4]  + rij * cc_3b_ij[2][5]  + rij_sq * cc_3b_ij[2][6]  + rij_th * cc_3b_ij[2][7];
        basis_ij[3] = cc_3b_ij[3][0]  + rij * cc_3b_ij[3][1]  + rij_sq * cc_3b_ij[3][2]  + rij_th * cc_3b_ij[3][3];

        // -------- basis_ik --------
        auto cc_3b_ik = &(cached_constants_3b[1][iknot_ik]);
        basis_ik[0] = cc_3b_ik[0][12] + rik * cc_3b_ik[0][13] + rik_sq * cc_3b_ik[0][14] + rik_th * cc_3b_ik[0][15];
        basis_ik[1] = cc_3b_ik[1][8]  + rik * cc_3b_ik[1][9]  + rik_sq * cc_3b_ik[1][10] + rik_th * cc_3b_ik[1][11];
        basis_ik[2] = cc_3b_ik[2][4]  + rik * cc_3b_ik[2][5]  + rik_sq * cc_3b_ik[2][6]  + rik_th * cc_3b_ik[2][7];
        basis_ik[3] = cc_3b_ik[3][0]  + rik * cc_3b_ik[3][1]  + rik_sq * cc_3b_ik[3][2]  + rik_th * cc_3b_ik[3][3];

        // -------- basis_jk --------
        auto cc_3b_jk = &(cached_constants_3b[2][iknot_jk]);
        basis_jk[0] = cc_3b_jk[0][12] + rjk * cc_3b_jk[0][13] + rjk_sq * cc_3b_jk[0][14] + rjk_th * cc_3b_jk[0][15];
        basis_jk[1] = cc_3b_jk[1][8]  + rjk * cc_3b_jk[1][9]  + rjk_sq * cc_3b_jk[1][10] + rjk_th * cc_3b_jk[1][11];
        basis_jk[2] = cc_3b_jk[2][4]  + rjk * cc_3b_jk[2][5]  + rjk_sq * cc_3b_jk[2][6]  + rjk_th * cc_3b_jk[2][7];
        basis_jk[3] = cc_3b_jk[3][0]  + rjk * cc_3b_jk[3][1]  + rjk_sq * cc_3b_jk[3][2]  + rjk_th * cc_3b_jk[3][3];


        const int base_offset = type_offset_3b[map_to];
        const int K_m = uf3_potential->n3b_coeff_array_size[map_to][1]; // Total basis functions for r_ik
        const int K_n = uf3_potential->n3b_coeff_array_size[map_to][0]; // Total basis functions for r_jk
        for (int l = 0; l < 4; l++) {
          const double basis_ij_i = basis_ij[l];
          for (int m = 0; m < 4; m++) {
            const double factor = basis_ij_i * basis_ik[m];
            const int col_offset = base_offset + (iknot_ij + l) * K_m * K_n + (iknot_ik + m) * K_n + iknot_jk;
            //if ( col_offset < leading_trim || col_offset > size_array_cols-trailing_trim-3 ) continue;
            array_local[0][col_offset]   += factor * basis_jk[0];
            array_local[0][col_offset+1] += factor * basis_jk[1];
            array_local[0][col_offset+2] += factor * basis_jk[2];
            array_local[0][col_offset+3] += factor * basis_jk[3];
          }
        }

        // -------- basis_ij_der --------
        auto cc_3b_deri_ij = &(cached_constants_3b_deri[0][iknot_ij]);
        basis_ij_der[0] = cc_3b_deri_ij[0][6] + rij * cc_3b_deri_ij[0][7] + rij_sq * cc_3b_deri_ij[0][8];
        basis_ij_der[1] = cc_3b_deri_ij[1][3] + rij * cc_3b_deri_ij[1][4] + rij_sq * cc_3b_deri_ij[1][5];
        basis_ij_der[2] = cc_3b_deri_ij[2][0] + rij * cc_3b_deri_ij[2][1] + rij_sq * cc_3b_deri_ij[2][2];

        // -------- basis_ik_der --------
        auto cc_3b_deri_ik = &(cached_constants_3b_deri[1][iknot_ik]);
        basis_ik_der[0] = cc_3b_deri_ik[0][6] + rik * cc_3b_deri_ik[0][7] + rik_sq * cc_3b_deri_ik[0][8];
        basis_ik_der[1] = cc_3b_deri_ik[1][3] + rik * cc_3b_deri_ik[1][4] + rik_sq * cc_3b_deri_ik[1][5];
        basis_ik_der[2] = cc_3b_deri_ik[2][0] + rik * cc_3b_deri_ik[2][1] + rik_sq * cc_3b_deri_ik[2][2];

        // -------- basis_jk_der --------
        auto cc_3b_deri_jk = &(cached_constants_3b_deri[2][iknot_jk]);
        basis_jk_der[0] = cc_3b_deri_jk[0][6] + rjk * cc_3b_deri_jk[0][7] + rjk_sq * cc_3b_deri_jk[0][8];
        basis_jk_der[1] = cc_3b_deri_jk[1][3] + rjk * cc_3b_deri_jk[1][4] + rjk_sq * cc_3b_deri_jk[1][5];
        basis_jk_der[2] = cc_3b_deri_jk[2][0] + rjk * cc_3b_deri_jk[2][1] + rjk_sq * cc_3b_deri_jk[2][2];

        double triangle_eval1 = 0.0;
        for (int l = 0; l < 3; l++) {
          const double basis_ij_der_i = basis_ij_der[l];
          for (int m = 0; m < 4; m++) {
            const double factor = basis_ij_der_i * basis_ik[m];
            triangle_eval1 += factor * (basis_jk[0] + basis_jk[1] + basis_jk[2] + basis_jk[3]);
          }
        }

        double triangle_eval2 = 0.0;
        for (int l = 0; l < 4; l++) {
          const double basis_ij_i = basis_ij[l];
          for (int m = 0; m < 3; m++) {
            const double factor = basis_ij_i * basis_ik_der[m];
            triangle_eval2 += factor * (basis_jk[0] + basis_jk[1] + basis_jk[2] + basis_jk[3]);
          }
        }

        double triangle_eval3 = 0.0;
        for (int l = 0; l < 4; l++) {
          const double basis_ij_i = basis_ij[l];
          for (int m = 0; m < 4; m++) {
            const double factor = basis_ij_i * basis_ik[m];
            triangle_eval3 += factor * (basis_jk_der[0] + basis_jk_der[1] + basis_jk_der[2]);
          }
        }

        const double fij0 = triangle_eval1 * del_rji[0] / rij;
        const double fik0 = triangle_eval2 * del_rki[0] / rik;
        const double fjk0 = triangle_eval3 * del_rkj[0] / rjk;

        const double fij1 = triangle_eval1 * del_rji[1] / rij;
        const double fik1 = triangle_eval2 * del_rki[1] / rik;
        const double fjk1 = triangle_eval3 * del_rkj[1] / rjk;

        const double fij2 = triangle_eval1 * del_rji[2] / rij;
        const double fik2 = triangle_eval2 * del_rki[2] / rik;
        const double fjk2 = triangle_eval3 * del_rkj[2] / rjk;

        double Fi[3], Fj[3], Fk[3];

        Fi[0] = fij0 + fik0;
        Fi[1] = fij1 + fik1;
        Fi[2] = fij2 + fik2;

        Fj[0] = -fij0 + fjk0;
        Fj[1] = -fij1 + fjk1;
        Fj[2] = -fij2 + fjk2;

        Fk[0] = -(fik0 + fjk0);
        Fk[1] = -(fik1 + fjk1);
        Fk[2] = -(fik2 + fjk2);

/*
        array_local[row_offset_i  ][col_offset] += Fi[0];
        array_local[row_offset_i+1][col_offset] += Fi[1];
        array_local[row_offset_i+2][col_offset] += Fi[2];

        array_local[row_offset_j  ][col_offset] += Fj[0];
        array_local[row_offset_j+1][col_offset] += Fj[1];
        array_local[row_offset_j+2][col_offset] += Fj[2];

        array_local[row_offset_k  ][col_offset] += Fk[0];
        array_local[row_offset_k+1][col_offset] += Fk[1];
        array_local[row_offset_k+2][col_offset] += Fk[2];
*/



      }
    }

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


