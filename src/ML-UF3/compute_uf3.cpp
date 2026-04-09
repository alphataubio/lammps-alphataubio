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
#include <cmath>

using namespace LAMMPS_NS;

ComputeUF3::ComputeUF3(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), cutsq(nullptr), setflag(nullptr),
    neighshort(nullptr), list(nullptr), array_local(nullptr), 
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

  virial_flag = utils::logical(FLERR, arg[4], false, lmp);
  if (virial_flag) size_array_rows = 1 + 3*(atom->natoms) + 6;
  else size_array_rows = 1 + 3*(atom->natoms);

  const int np1 = atom->ntypes + 1;
  memory->create(setflag, np1, np1, "uf3:setflag");
  memory->create(cutsq, np1, np1, "uf3:cutsq");
  if (pot_3b) {
    maxshort = 20;
    memory->create(neighshort, maxshort, "uf3:neighshort");
  }

  std::vector<std::string> elements_(1); 
  for(int i=6; i<narg ; i++) elements_.push_back(arg[i]);
  uf3_potential = new UF3Potential(lmp, arg[5], cutsq, setflag, elements_, pot_3b);

  lastcol = 0;
  const int ntypes = atom->ntypes;
  for (int i = 1; i <= ntypes; i++) {
    for (int j = i; j <= ntypes; j++) {
      const double cut_2b_ij = uf3_potential->cut_2b[i][j];
      cutsq[i][j] = cutsq[j][i] = cut_2b_ij * cut_2b_ij;
      for (int l = 0; l < uf3_potential->n2b_coeff_array_size[i][j]; l++) {
        int idx = static_cast<int>(std::round(uf3_potential->n2b_coeff_array[i][j][l]));
        if (idx > lastcol) lastcol = idx;
      }
    }
  }

  if (pot_3b) {
    for (int i = 1; i <= ntypes; i++) {
      for (int j = i; j <= ntypes; j++) {
        for (int k = j; k <= ntypes; k++) {
          int map_to = uf3_potential->map_3b[i][j][k];
          auto n3b_size = uf3_potential->n3b_coeff_array_size[map_to];
          for (int l = 0; l < n3b_size[0]; l++) {
            for (int m = 0; m < n3b_size[1]; m++) {
              for (int n = 0; n < n3b_size[2]; n++) {
                int idx = static_cast<int>(std::round(uf3_potential->n3b_coeff_array[map_to][l][m][n]));
                if (idx > lastcol) lastcol = idx;
              }
            }
          }
        }
      }
    }
  }
  lastcol++; // reference energy/forces last column
  size_array_cols = lastcol + 1;
  //fprintf(stderr, "*** Automatically sized descriptor matrix to %i columns\n", size_array_cols);
  memory->create(array_local, size_array_rows, size_array_cols, "uf3:array_local");
  memory->create(array, size_array_rows, size_array_cols, "uf3:array");

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
  delete uf3_potential;
}

/* ---------------------------------------------------------------------- */

void ComputeUF3::init()
{
  if (force->pair == nullptr)
    error->all(FLERR,"Compute uf3 requires a pair style be defined");

  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);

  if (modify->get_compute_by_style("uf3").size() > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute uf3");


  c_pe = modify->get_compute_by_id("thermo_pe");
  if (!c_pe) error->all(FLERR,"Compute thermo_pe does not exist.");

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
  invoked_array = update->ntimestep;

  for (int i = 0; i < size_array_rows; i++) {
    for (int j = 0; j < size_array_cols; j++) array_local[i][j] = 0.0;
  }

  neighbor->build_one(list);

  const int inum = list->inum;
  const int* const ilist = list->ilist;
  const int* const numneigh = list->numneigh;
  int** const firstneigh = list->firstneigh;
  const int* const type = atom->type;
  double** const x = atom->x;
  const int* const mask = atom->mask;

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

      // --- 2-BODY DESCRIPTORS ---
      double **cc_2b = &(uf3_potential->cached_constants_2b[itype][jtype][start_idx-3]);

      for (int m = 0; m < 4; m++) {
        // Direct sparse mapping from the .pot file
        double map_val = uf3_potential->n2b_coeff_array[itype][jtype][start_idx - 3 + m];
        if (map_val < -0.5) continue; // Python dropped this column
        const int sparse_col = static_cast<int>(std::round(map_val));

        const int n = (3 - m) * 4;
        const double basis_val = cc_2b[m][n] + rij * cc_2b[m][n+1] + rsq * cc_2b[m][n+2] + rth * cc_2b[m][n+3];
        const double b_der = cc_2b[m][n+1] + 2.0 * rij * cc_2b[m][n+2] + 3.0 * rsq * cc_2b[m][n+3];
        
        const double force_factor = -b_der / rij;
        const double fx_feature = delx * force_factor;
        const double fy_feature = dely * force_factor;
        const double fz_feature = delz * force_factor;

        // Energy Accumulation
        array_local[0][sparse_col] += basis_val;
        
        // Force Accumulation
        array_local[row_offset_i  ][sparse_col] += fx_feature;
        array_local[row_offset_i+1][sparse_col] += fy_feature;
        array_local[row_offset_i+2][sparse_col] += fz_feature;
        
        array_local[row_offset_j  ][sparse_col] -= fx_feature;
        array_local[row_offset_j+1][sparse_col] -= fy_feature;
        array_local[row_offset_j+2][sparse_col] -= fz_feature;

        // Virials
        if (virial_flag) {
          array_local[size_array_rows-6][sparse_col] += delx * fx_feature; 
          array_local[size_array_rows-5][sparse_col] += dely * fy_feature; 
          array_local[size_array_rows-4][sparse_col] += delz * fz_feature; 
          array_local[size_array_rows-3][sparse_col] += delz * fy_feature; 
          array_local[size_array_rows-2][sparse_col] += delz * fx_feature; 
          array_local[size_array_rows-1][sparse_col] += dely * fx_feature; 
        }
      }
    } 

    // --- 3-BODY DESCRIPTORS ---
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

        // Evaluate pure basis functions using the cache
        auto cc_3b_ij = &(cached_constants_3b[0][iknot_ij]);
        basis_ij[0] = cc_3b_ij[0][12] + rij * cc_3b_ij[0][13] + rij_sq * cc_3b_ij[0][14] + rij_th * cc_3b_ij[0][15];
        basis_ij[1] = cc_3b_ij[1][8]  + rij * cc_3b_ij[1][9]  + rij_sq * cc_3b_ij[1][10] + rij_th * cc_3b_ij[1][11];
        basis_ij[2] = cc_3b_ij[2][4]  + rij * cc_3b_ij[2][5]  + rij_sq * cc_3b_ij[2][6]  + rij_th * cc_3b_ij[2][7];
        basis_ij[3] = cc_3b_ij[3][0]  + rij * cc_3b_ij[3][1]  + rij_sq * cc_3b_ij[3][2]  + rij_th * cc_3b_ij[3][3];

        auto cc_3b_ik = &(cached_constants_3b[1][iknot_ik]);
        basis_ik[0] = cc_3b_ik[0][12] + rik * cc_3b_ik[0][13] + rik_sq * cc_3b_ik[0][14] + rik_th * cc_3b_ik[0][15];
        basis_ik[1] = cc_3b_ik[1][8]  + rik * cc_3b_ik[1][9]  + rik_sq * cc_3b_ik[1][10] + rik_th * cc_3b_ik[1][11];
        basis_ik[2] = cc_3b_ik[2][4]  + rik * cc_3b_ik[2][5]  + rik_sq * cc_3b_ik[2][6]  + rik_th * cc_3b_ik[2][7];
        basis_ik[3] = cc_3b_ik[3][0]  + rik * cc_3b_ik[3][1]  + rik_sq * cc_3b_ik[3][2]  + rik_th * cc_3b_ik[3][3];

        auto cc_3b_jk = &(cached_constants_3b[2][iknot_jk]);
        basis_jk[0] = cc_3b_jk[0][12] + rjk * cc_3b_jk[0][13] + rjk_sq * cc_3b_jk[0][14] + rjk_th * cc_3b_jk[0][15];
        basis_jk[1] = cc_3b_jk[1][8]  + rjk * cc_3b_jk[1][9]  + rjk_sq * cc_3b_jk[1][10] + rjk_th * cc_3b_jk[1][11];
        basis_jk[2] = cc_3b_jk[2][4]  + rjk * cc_3b_jk[2][5]  + rjk_sq * cc_3b_jk[2][6]  + rjk_th * cc_3b_jk[2][7];
        basis_jk[3] = cc_3b_jk[3][0]  + rjk * cc_3b_jk[3][1]  + rjk_sq * cc_3b_jk[3][2]  + rjk_th * cc_3b_jk[3][3];

        auto cc_3b_deri_ij = &(cached_constants_3b_deri[0][iknot_ij]);
        basis_ij_der[0] = cc_3b_deri_ij[0][6] + rij * cc_3b_deri_ij[0][7] + rij_sq * cc_3b_deri_ij[0][8];
        basis_ij_der[1] = cc_3b_deri_ij[1][3] + rij * cc_3b_deri_ij[1][4] + rij_sq * cc_3b_deri_ij[1][5];
        basis_ij_der[2] = cc_3b_deri_ij[2][0] + rij * cc_3b_deri_ij[2][1] + rij_sq * cc_3b_deri_ij[2][2];

        auto cc_3b_deri_ik = &(cached_constants_3b_deri[1][iknot_ik]);
        basis_ik_der[0] = cc_3b_deri_ik[0][6] + rik * cc_3b_deri_ik[0][7] + rik_sq * cc_3b_deri_ik[0][8];
        basis_ik_der[1] = cc_3b_deri_ik[1][3] + rik * cc_3b_deri_ik[1][4] + rik_sq * cc_3b_deri_ik[1][5];
        basis_ik_der[2] = cc_3b_deri_ik[2][0] + rik * cc_3b_deri_ik[2][1] + rik_sq * cc_3b_deri_ik[2][2];

        auto cc_3b_deri_jk = &(cached_constants_3b_deri[2][iknot_jk]);
        basis_jk_der[0] = cc_3b_deri_jk[0][6] + rjk * cc_3b_deri_jk[0][7] + rjk_sq * cc_3b_deri_jk[0][8];
        basis_jk_der[1] = cc_3b_deri_jk[1][3] + rjk * cc_3b_deri_jk[1][4] + rjk_sq * cc_3b_deri_jk[1][5];
        basis_jk_der[2] = cc_3b_deri_jk[2][0] + rjk * cc_3b_deri_jk[2][1] + rjk_sq * cc_3b_deri_jk[2][2];

        const double d_bij[4] = {basis_ij_der[0], basis_ij_der[1], basis_ij_der[2], 0.0};
        const double d_bik[4] = {basis_ik_der[0], basis_ik_der[1], basis_ik_der[2], 0.0};
        const double d_bjk[4] = {basis_jk_der[0], basis_jk_der[1], basis_jk_der[2], 0.0};

        const double dx_ij = del_rji[0] / rij;
        const double dy_ij = del_rji[1] / rij;
        const double dz_ij = del_rji[2] / rij;

        const double dx_ik = del_rki[0] / rik;
        const double dy_ik = del_rki[1] / rik;
        const double dz_ik = del_rki[2] / rik;

        const double dx_jk = del_rkj[0] / rjk;
        const double dy_jk = del_rkj[1] / rjk;
        const double dz_jk = del_rkj[2] / rjk;

        for (int l = 0; l < 4; l++) {
          const double b_ij  = basis_ij[l];
          const double db_ij = d_bij[l];
          
          for (int m = 0; m < 4; m++) {
            const double b_ik  = basis_ik[m];
            const double db_ik = d_bik[m];
            
            const double term1_base = db_ij * b_ik; 
            const double term2_base = b_ij * db_ik; 
            const double term3_base = b_ij * b_ik;  

            for (int n = 0; n < 4; n++) {
              
              // Direct sparse mapping from the .pot file
              const double map_val = uf3_potential->n3b_coeff_array[map_to][iknot_ij + l][iknot_ik + m][iknot_jk + n];
              if (map_val < -0.5) continue; // Python dropped this column

              const int sparse_col = static_cast<int>(std::round(map_val));
              
              const double b_jk  = basis_jk[n];
              const double db_jk = d_bjk[n];

              // --- Energy Descriptor ---
              array_local[0][sparse_col] += term3_base * b_jk;

              // --- Force Descriptors ---
              const double d_ij_part = term1_base * b_jk;
              const double d_ik_part = term2_base * b_jk;
              const double d_jk_part = term3_base * db_jk;

              const double fij_x = d_ij_part * dx_ij;
              const double fik_x = d_ik_part * dx_ik;
              const double fjk_x = d_jk_part * dx_jk;

              const double fij_y = d_ij_part * dy_ij;
              const double fik_y = d_ik_part * dy_ik;
              const double fjk_y = d_jk_part * dy_jk;

              const double fij_z = d_ij_part * dz_ij;
              const double fik_z = d_ik_part * dz_ik;
              const double fjk_z = d_jk_part * dz_jk;

              // Pre-calculate specific atom force components
              const double Fi_x = fij_x + fik_x;
              const double Fi_y = fij_y + fik_y;
              const double Fi_z = fij_z + fik_z;

              const double Fj_x = -fij_x + fjk_x;
              const double Fj_y = -fij_y + fjk_y;
              const double Fj_z = -fij_z + fjk_z;

              const double Fk_x = -(fik_x + fjk_x);
              const double Fk_y = -(fik_y + fjk_y);
              const double Fk_z = -(fik_z + fjk_z);

              // Force Accumulation
              array_local[row_offset_i    ][sparse_col] += Fi_x;
              array_local[row_offset_i + 1][sparse_col] += Fi_y;
              array_local[row_offset_i + 2][sparse_col] += Fi_z;

              array_local[row_offset_j    ][sparse_col] += Fj_x;
              array_local[row_offset_j + 1][sparse_col] += Fj_y;
              array_local[row_offset_j + 2][sparse_col] += Fj_z;

              array_local[row_offset_k    ][sparse_col] += Fk_x;
              array_local[row_offset_k + 1][sparse_col] += Fk_y;
              array_local[row_offset_k + 2][sparse_col] += Fk_z;

              // 3-Body Virial Accumulation
              if (virial_flag) {
                array_local[size_array_rows-6][sparse_col] += (del_rji[0] * Fj_x) + (del_rki[0] * Fk_x); // W_xx
                array_local[size_array_rows-5][sparse_col] += (del_rji[1] * Fj_y) + (del_rki[1] * Fk_y); // W_yy
                array_local[size_array_rows-4][sparse_col] += (del_rji[2] * Fj_z) + (del_rki[2] * Fk_z); // W_zz
                array_local[size_array_rows-3][sparse_col] += (del_rji[2] * Fj_y) + (del_rki[2] * Fk_y); // W_zy
                array_local[size_array_rows-2][sparse_col] += (del_rji[2] * Fj_x) + (del_rki[2] * Fk_x); // W_zx
                array_local[size_array_rows-1][sparse_col] += (del_rji[1] * Fj_x) + (del_rki[1] * Fk_x); // W_yx
              }

            }
          }
        }
      }
    }
  }

  // Assign reference force targets
  for (int i = 0; i < atom->nlocal; i++) {
    int iglobal = atom->tag[i];
    int irow = 3*(iglobal-1)+1;
    array_local[irow++][lastcol] = atom->f[i][0];
    array_local[irow++][lastcol] = atom->f[i][1];
    array_local[irow][lastcol] = atom->f[i][2];
  }

  MPI_Allreduce(&array_local[0][0],&array[0][0],size_array_rows*size_array_cols,MPI_DOUBLE,MPI_SUM,world);

  // Assign reference energy target
  array[0][lastcol] = c_pe->compute_scalar();

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

double ComputeUF3::memory_usage()
{
  return (double)size_array_rows*size_array_cols*sizeof(double)*2;
}
