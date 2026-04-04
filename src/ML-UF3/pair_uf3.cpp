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

/* ----------------------------------------------------------------------
   Contributing author: Ajinkya Hire (Univ. of Florida),
                        Hendrik Krass (Univ. of Constance),
                        Matthias Rupp (Luxembourg Institute of Science and Technology),
                        Richard Hennig (Univ of Florida)
---------------------------------------------------------------------- */

#include "pair_uf3.h"
#include "uf3_potential.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <algorithm>
#include <cmath>

using namespace LAMMPS_NS;
using MathConst::THIRD;

/* ---------------------------------------------------------------------- */

PairUF3::PairUF3(LAMMPS *lmp) : Pair(lmp), uf3_potential(nullptr), neighshort(nullptr)
{
  single_enable = 1;    // 1 if single() routine exists
  one_coeff = 1;        // 1 if allows only one coeff * * call
  restartinfo = 0;      // 1 if pair style writes restart info
  maxshort = 20;
  centroidstressflag = CENTROID_AVAIL;
  manybody_flag = 1;
  pot_3b = false;
}

/* ---------------------------------------------------------------------- */

PairUF3::~PairUF3()
{
  if (copymode) return;
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    if (pot_3b) memory->destroy(neighshort);
  }
  delete uf3_potential;
}

/* ----------------------------------------------------------------------
 *     global settings
 * ---------------------------------------------------------------------- */

void PairUF3::settings(int narg, char **arg)
{

  if (narg != 1)
    error->all(FLERR,
               "Invalid number of arguments for pair_style uf3. "
               "Are you using a 2-body or 2 & 3-body UF potential?");
  const int nbody = utils::inumeric(FLERR, arg[0], true, lmp);
  if (nbody == 2) {
    pot_3b = false;
    manybody_flag = 0;
  } else if (nbody == 3) {
    pot_3b = true;
    single_enable = 0;
  } else
    error->all(FLERR, "Pair style uf3 not (yet) implemented for {}-body terms", nbody);
}

/* ----------------------------------------------------------------------
 *    set coeffs for one or more type pairs
 * ---------------------------------------------------------------------- */
void PairUF3::coeff(int narg, char **arg)
{
  if (narg != 3 + atom->ntypes)
    error->all(FLERR, "Invalid number of arguments uf3 in pair coeffs.");
  if (!allocated) allocate();
  std::vector<std::string> elements_(narg - 3);
  for(int i=3; i<narg ; i++) elements_.push_back(arg[i]);
  uf3_potential = new UF3Potential(lmp, arg[2], cutsq, setflag, elements_, pot_3b);
}

/* ---------------------------------------------------------------------- */

void PairUF3::allocate()
{
  allocated = 1;
  const int np1 = atom->ntypes + 1;
  map = new int[np1];    //No need to delete map as ~Pair deletes map
  memory->create(setflag, np1, np1, "pair:setflag");
  memory->create(cutsq, np1, np1, "pair:cutsq");
  if (pot_3b) memory->create(neighshort, maxshort, "pair:neighshort");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */
void PairUF3::init_style()
{
  if (force->newton_pair == 0) error->all(FLERR, "UF3: Pair style requires newton pair on");
  // request a default neighbor list
  neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ----------------------------------------------------------------------
   init list sets the pointer to full neighbour list requested in previous function
------------------------------------------------------------------------- */

void PairUF3::init_list(int /*id*/, class NeighList *ptr)
{
  list = ptr;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */
double PairUF3::init_one(int i, int j)
{
  //init_one is called by pair.cpp at line 267 where it is squred
  //at line 268
  return uf3_potential->cut_2b[i][j];
}

/* ---------------------------------------------------------------------- */

void PairUF3::compute(int eflag, int vflag)
{

  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  const int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  // loop over neighbors of my atoms
  for (int ii = 0; ii < inum; ii++) {
    double evdwl = 0;
    const int i = ilist[ii];
    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];
    const int itype = type[i];
    int *jlist = firstneigh[i];
    int numshort = 0;
    for (int jj = 0; jj < numneigh[i]; jj++) {
      const int j = jlist[jj] & NEIGHMASK;
      const double delx = xtmp - x[j][0];
      const double dely = ytmp - x[j][1];
      const double delz = ztmp - x[j][2];
      const double rsq = delx * delx + dely * dely + delz * delz;
      const int jtype = type[j];
      if (rsq < cutsq[itype][jtype]) {
        const double rij = sqrt(rsq);

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

        const int knot_start_index = uf3_potential->get_starting_index_2b(itype, jtype, rij);
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

        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;
        f[j][0] -= fx;
        f[j][1] -= fy;
        f[j][2] -= fz;

        if (eflag) {
          const double rth = rsq * rij;
          double **cached_constants_2b = uf3_potential->cached_constants_2b[itype][jtype];
          evdwl =        cached_constants_2b[knot_start_index    ][0];
          evdwl += rij * cached_constants_2b[knot_start_index    ][1];
          evdwl += rsq * cached_constants_2b[knot_start_index    ][2];
          evdwl += rth * cached_constants_2b[knot_start_index    ][3];
          evdwl +=       cached_constants_2b[knot_start_index - 1][4];
          evdwl += rij * cached_constants_2b[knot_start_index - 1][5];
          evdwl += rsq * cached_constants_2b[knot_start_index - 1][6];
          evdwl += rth * cached_constants_2b[knot_start_index - 1][7];
          evdwl +=       cached_constants_2b[knot_start_index - 2][8];
          evdwl += rij * cached_constants_2b[knot_start_index - 2][9];
          evdwl += rsq * cached_constants_2b[knot_start_index - 2][10];
          evdwl += rth * cached_constants_2b[knot_start_index - 2][11];
          evdwl +=       cached_constants_2b[knot_start_index - 3][12];
          evdwl += rij * cached_constants_2b[knot_start_index - 3][13];
          evdwl += rsq * cached_constants_2b[knot_start_index - 3][14];
          evdwl += rth * cached_constants_2b[knot_start_index - 3][15];
        }

        if (evflag) {
          ev_tally_xyz(i, j, nlocal, newton_pair, evdwl, 0.0, fx, fy, fz, delx, dely, delz);

          // Centroid Stress
          if (vflag_either && cvflag_atom) {
            double v[6];

            v[0] = delx * fx;
            v[1] = dely * fy;
            v[2] = delz * fz;
            v[3] = delx * fy;
            v[4] = delx * fz;
            v[5] = dely * fz;

            cvatom[i][0] += 0.5 * v[0];
            cvatom[i][1] += 0.5 * v[1];
            cvatom[i][2] += 0.5 * v[2];
            cvatom[i][3] += 0.5 * v[3];
            cvatom[i][4] += 0.5 * v[4];
            cvatom[i][5] += 0.5 * v[5];
            cvatom[i][6] += 0.5 * v[3];
            cvatom[i][7] += 0.5 * v[4];
            cvatom[i][8] += 0.5 * v[5];

            cvatom[j][0] += 0.5 * v[0];
            cvatom[j][1] += 0.5 * v[1];
            cvatom[j][2] += 0.5 * v[2];
            cvatom[j][3] += 0.5 * v[3];
            cvatom[j][4] += 0.5 * v[4];
            cvatom[j][5] += 0.5 * v[5];
            cvatom[j][6] += 0.5 * v[3];
            cvatom[j][7] += 0.5 * v[4];
            cvatom[j][8] += 0.5 * v[5];
          }
        }
      }
    }

    // 3-body interaction
    // jth atom

    for (int jj = 0; jj < numshort - 1; jj++) {

      double del_rji[3], del_rki[3], del_rkj[3];

      const int j = neighshort[jj];
      const int jtype = type[j];
      del_rji[0] = x[j][0] - xtmp;
      del_rji[1] = x[j][1] - ytmp;
      del_rji[2] = x[j][2] - ztmp;
      const double rij_sq = (del_rji[0] * del_rji[0]) + (del_rji[1] * del_rji[1]) + (del_rji[2] * del_rji[2]);
      const double rij = sqrt(rij_sq);

      // kth atom
      for (int kk = jj + 1; kk < numshort; kk++) {

        const int k = neighshort[kk];
        const int ktype = type[k];

        del_rki[0] = x[k][0] - xtmp;
        del_rki[1] = x[k][1] - ytmp;
        del_rki[2] = x[k][2] - ztmp;
        const double rik_sq = (del_rki[0] * del_rki[0]) + (del_rki[1] * del_rki[1]) + (del_rki[2] * del_rki[2]);
        const double rik = sqrt(rik_sq);

        if ((rij <= uf3_potential->cut_3b[itype][jtype][ktype]) &&
            (rik <= uf3_potential->cut_3b[itype][ktype][jtype]) &&
            (rij >= uf3_potential->min_cut_3b[itype][jtype][ktype][2]) &&
            (rik >= uf3_potential->min_cut_3b[itype][jtype][ktype][1])) {

          del_rkj[0] = x[k][0] - x[j][0];
          del_rkj[1] = x[k][1] - x[j][1];
          del_rkj[2] = x[k][2] - x[j][2];

          const double rjk_sq =(del_rkj[0] * del_rkj[0]) + (del_rkj[1] * del_rkj[1]) + (del_rkj[2] * del_rkj[2]);
          const double rjk = sqrt(rjk_sq);

          if (rjk >= uf3_potential->min_cut_3b[itype][jtype][ktype][0]) {
            const double rij_th = rij * rij_sq;
            const double rik_th = rik * rik_sq;
            const double rjk_th = rjk * rjk_sq;

            const int map_to = uf3_potential->map_3b[itype][jtype][ktype];
            double ***cached_constants_3b = uf3_potential->cached_constants_3b[map_to];
            double ***cached_constants_3b_deri = uf3_potential->cached_constants_3b_deri[map_to];

            const int knot_start_index_ij = uf3_potential->get_starting_index_3b(itype, jtype, ktype, rij, 2);
            const int knot_start_index_ik = uf3_potential->get_starting_index_3b(itype, jtype, ktype, rik, 1);
            const int knot_start_index_jk = uf3_potential->get_starting_index_3b(itype, jtype, ktype, rjk, 0);
            double basis_ij[4], basis_ik[4], basis_jk[4], basis_ij_der[3], basis_ik_der[3], basis_jk_der[3];

            //--------------basis_ij
            basis_ij[0] =           cached_constants_3b[0][knot_start_index_ij - 3][12];
            basis_ij[0] += rij    * cached_constants_3b[0][knot_start_index_ij - 3][13];
            basis_ij[0] += rij_sq * cached_constants_3b[0][knot_start_index_ij - 3][14];
            basis_ij[0] += rij_th * cached_constants_3b[0][knot_start_index_ij - 3][15];

            basis_ij[1] =           cached_constants_3b[0][knot_start_index_ij - 2][8];
            basis_ij[1] += rij    * cached_constants_3b[0][knot_start_index_ij - 2][9];
            basis_ij[1] += rij_sq * cached_constants_3b[0][knot_start_index_ij - 2][10];
            basis_ij[1] += rij_th * cached_constants_3b[0][knot_start_index_ij - 2][11];

            basis_ij[2] =           cached_constants_3b[0][knot_start_index_ij - 1][4];
            basis_ij[2] += rij    * cached_constants_3b[0][knot_start_index_ij - 1][5];
            basis_ij[2] += rij_sq * cached_constants_3b[0][knot_start_index_ij - 1][6];
            basis_ij[2] += rij_th * cached_constants_3b[0][knot_start_index_ij - 1][7];

            basis_ij[3] =           cached_constants_3b[0][knot_start_index_ij    ][0];
            basis_ij[3] += rij    * cached_constants_3b[0][knot_start_index_ij    ][1];
            basis_ij[3] += rij_sq * cached_constants_3b[0][knot_start_index_ij    ][2];
            basis_ij[3] += rij_th * cached_constants_3b[0][knot_start_index_ij    ][3];

            //--------------basis_ik
            basis_ik[0] =           cached_constants_3b[1][knot_start_index_ik - 3][12];
            basis_ik[0] += rik    * cached_constants_3b[1][knot_start_index_ik - 3][13];
            basis_ik[0] += rik_sq * cached_constants_3b[1][knot_start_index_ik - 3][14];
            basis_ik[0] += rik_th * cached_constants_3b[1][knot_start_index_ik - 3][15];

            basis_ik[1] =           cached_constants_3b[1][knot_start_index_ik - 2][8];
            basis_ik[1] += rik    * cached_constants_3b[1][knot_start_index_ik - 2][9];
            basis_ik[1] += rik_sq * cached_constants_3b[1][knot_start_index_ik - 2][10];
            basis_ik[1] += rik_th * cached_constants_3b[1][knot_start_index_ik - 2][11];

            basis_ik[2] =           cached_constants_3b[1][knot_start_index_ik - 1][4];
            basis_ik[2] += rik    * cached_constants_3b[1][knot_start_index_ik - 1][5];
            basis_ik[2] += rik_sq * cached_constants_3b[1][knot_start_index_ik - 1][6];
            basis_ik[2] += rik_th * cached_constants_3b[1][knot_start_index_ik - 1][7];

            basis_ik[3] =           cached_constants_3b[1][knot_start_index_ik    ][0];
            basis_ik[3] += rik    * cached_constants_3b[1][knot_start_index_ik    ][1];
            basis_ik[3] += rik_sq * cached_constants_3b[1][knot_start_index_ik    ][2];
            basis_ik[3] += rik_th * cached_constants_3b[1][knot_start_index_ik    ][3];

            //--------------basis_jk
            basis_jk[0] =           cached_constants_3b[2][knot_start_index_jk - 3][12];
            basis_jk[0] += rjk    * cached_constants_3b[2][knot_start_index_jk - 3][13];
            basis_jk[0] += rjk_sq * cached_constants_3b[2][knot_start_index_jk - 3][14];
            basis_jk[0] += rjk_th * cached_constants_3b[2][knot_start_index_jk - 3][15];

            basis_jk[1] =           cached_constants_3b[2][knot_start_index_jk - 2][8];
            basis_jk[1] += rjk    * cached_constants_3b[2][knot_start_index_jk - 2][9];
            basis_jk[1] += rjk_sq * cached_constants_3b[2][knot_start_index_jk - 2][10];
            basis_jk[1] += rjk_th * cached_constants_3b[2][knot_start_index_jk - 2][11];

            basis_jk[2] =           cached_constants_3b[2][knot_start_index_jk - 1][4];
            basis_jk[2] += rjk    * cached_constants_3b[2][knot_start_index_jk - 1][5];
            basis_jk[2] += rjk_sq * cached_constants_3b[2][knot_start_index_jk - 1][6];
            basis_jk[2] += rjk_th * cached_constants_3b[2][knot_start_index_jk - 1][7];

            basis_jk[3] =           cached_constants_3b[2][knot_start_index_jk    ][0];
            basis_jk[3] += rjk    * cached_constants_3b[2][knot_start_index_jk    ][1];
            basis_jk[3] += rjk_sq * cached_constants_3b[2][knot_start_index_jk    ][2];
            basis_jk[3] += rjk_th * cached_constants_3b[2][knot_start_index_jk    ][3];

            //----------------basis_ij_der
            basis_ij_der[0] =           cached_constants_3b_deri[0][knot_start_index_ij - 3][6];
            basis_ij_der[0] += rij    * cached_constants_3b_deri[0][knot_start_index_ij - 3][7];
            basis_ij_der[0] += rij_sq * cached_constants_3b_deri[0][knot_start_index_ij - 3][8];

            basis_ij_der[1] =           cached_constants_3b_deri[0][knot_start_index_ij - 2][3];
            basis_ij_der[1] += rij    * cached_constants_3b_deri[0][knot_start_index_ij - 2][4];
            basis_ij_der[1] += rij_sq * cached_constants_3b_deri[0][knot_start_index_ij - 2][5];

            basis_ij_der[2] =           cached_constants_3b_deri[0][knot_start_index_ij - 1][0];
            basis_ij_der[2] += rij    * cached_constants_3b_deri[0][knot_start_index_ij - 1][1];
            basis_ij_der[2] += rij_sq * cached_constants_3b_deri[0][knot_start_index_ij - 1][2];

            //----------------basis_ik_der
            basis_ik_der[0] =           cached_constants_3b_deri[1][knot_start_index_ik - 3][6];
            basis_ik_der[0] += rik    * cached_constants_3b_deri[1][knot_start_index_ik - 3][7];
            basis_ik_der[0] += rik_sq * cached_constants_3b_deri[1][knot_start_index_ik - 3][8];

            basis_ik_der[1] =           cached_constants_3b_deri[1][knot_start_index_ik - 2][3];
            basis_ik_der[1] += rik    * cached_constants_3b_deri[1][knot_start_index_ik - 2][4];
            basis_ik_der[1] += rik_sq * cached_constants_3b_deri[1][knot_start_index_ik - 2][5];

            basis_ik_der[2] =           cached_constants_3b_deri[1][knot_start_index_ik - 1][0];
            basis_ik_der[2] += rik    * cached_constants_3b_deri[1][knot_start_index_ik - 1][1];
            basis_ik_der[2] += rik_sq * cached_constants_3b_deri[1][knot_start_index_ik - 1][2];

            //----------------basis_jk_der
            basis_jk_der[0] =           cached_constants_3b_deri[2][knot_start_index_jk - 3][6];
            basis_jk_der[0] += rjk    * cached_constants_3b_deri[2][knot_start_index_jk - 3][7];
            basis_jk_der[0] += rjk_sq * cached_constants_3b_deri[2][knot_start_index_jk - 3][8];

            basis_jk_der[1] =           cached_constants_3b_deri[2][knot_start_index_jk - 2][3];
            basis_jk_der[1] += rjk    * cached_constants_3b_deri[2][knot_start_index_jk - 2][4];
            basis_jk_der[1] += rjk_sq * cached_constants_3b_deri[2][knot_start_index_jk - 2][5];

            basis_jk_der[2] =           cached_constants_3b_deri[2][knot_start_index_jk - 1][0];
            basis_jk_der[2] += rjk    * cached_constants_3b_deri[2][knot_start_index_jk - 1][1];
            basis_jk_der[2] += rjk_sq * cached_constants_3b_deri[2][knot_start_index_jk - 1][2];

            double triangle_eval0 = 0.0;
            double triangle_eval1 = 0.0;
            double triangle_eval2 = 0.0;
            double triangle_eval3 = 0.0;
            const int iknot_ij = knot_start_index_ij - 3;
            const int iknot_ik = knot_start_index_ik - 3;
            const int iknot_jk = knot_start_index_jk - 3;

            for (int l = 0; l < 3; l++) {
              const double basis_ij_der_i = basis_ij_der[l];
              for (int m = 0; m < 4; m++) {
                const double factor = basis_ij_der_i * basis_ik[m];
                double *slice = &(uf3_potential->coeff_for_der_ij[map_to][iknot_ij + l][iknot_ik + m][iknot_jk]);
                const double tmp0 = slice[0] * basis_jk[0];
                const double tmp1 = slice[1] * basis_jk[1];
                const double tmp2 = slice[2] * basis_jk[2];
                const double tmp3 = slice[3] * basis_jk[3];
                triangle_eval1 += factor * (tmp0 + tmp1 + tmp2 + tmp3);
              }
            }

            for (int l = 0; l < 4; l++) {
              const double basis_ij_i = basis_ij[l];
              for (int m = 0; m < 3; m++) {
                const double factor = basis_ij_i * basis_ik_der[m];
                double *slice = &(uf3_potential->coeff_for_der_ik[map_to][iknot_ij + l][iknot_ik + m][iknot_jk]);
                const double tmp0 = slice[0] * basis_jk[0];
                const double tmp1 = slice[1] * basis_jk[1];
                const double tmp2 = slice[2] * basis_jk[2];
                const double tmp3 = slice[3] * basis_jk[3];
                triangle_eval2 += factor * (tmp0 + tmp1 + tmp2 + tmp3);
              }
            }

            for (int l = 0; l < 4; l++) {
              const double basis_ij_i = basis_ij[l];
              for (int m = 0; m < 4; m++) {
                const double factor = basis_ij_i * basis_ik[m];
                double *slice = &(uf3_potential->coeff_for_der_jk[map_to][iknot_ij + l][iknot_ik + m][iknot_jk]);
                const double tmp0 = slice[0] * basis_jk_der[0];
                const double tmp1 = slice[1] * basis_jk_der[1];
                const double tmp2 = slice[2] * basis_jk_der[2];
                triangle_eval3 += factor * (tmp0 + tmp1 + tmp2);
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
            f[i][0] += Fi[0];
            f[i][1] += Fi[1];
            f[i][2] += Fi[2];

            Fj[0] = -fij0 + fjk0;
            Fj[1] = -fij1 + fjk1;
            Fj[2] = -fij2 + fjk2;
            f[j][0] += Fj[0];
            f[j][1] += Fj[1];
            f[j][2] += Fj[2];

            Fk[0] = -(fik0 + fjk0);
            Fk[1] = -(fik1 + fjk1);
            Fk[2] = -(fik2 + fjk2);
            f[k][0] += Fk[0];
            f[k][1] += Fk[1];
            f[k][2] += Fk[2];

            if (eflag) {
              for (int l = 0; l < 4; l++) {
                const double basis_ij_i = basis_ij[l];
                for (int m = 0; m < 4; m++) {
                  const double factor = basis_ij_i * basis_ik[m];
                  const double *slice =
                      &(uf3_potential->n3b_coeff_array[map_to][iknot_ij + l][iknot_ik + m][iknot_jk]);
                  const double tmp0 = slice[0] * basis_jk[0];
                  const double tmp1 = slice[1] * basis_jk[1];
                  const double tmp2 = slice[2] * basis_jk[2];
                  const double tmp3 = slice[3] * basis_jk[3];
                  triangle_eval0 += factor * (tmp0 + tmp1 + tmp2 + tmp3);
                }
              }
              evdwl = triangle_eval0;
            }

            if (evflag) {
              ev_tally3(i, j, k, evdwl, 0, Fj, Fk, del_rji, del_rki);
              // Centroid stress 3-body term
              if (vflag_either && cvflag_atom) {

                const double ric0 = THIRD * (-del_rji[0] - del_rki[0]);
                const double ric1 = THIRD * (-del_rji[1] - del_rki[1]);
                const double ric2 = THIRD * (-del_rji[2] - del_rki[2]);
                cvatom[i][0] += ric0 * Fi[0];
                cvatom[i][1] += ric1 * Fi[1];
                cvatom[i][2] += ric2 * Fi[2];
                cvatom[i][3] += ric0 * Fi[1];
                cvatom[i][4] += ric0 * Fi[2];
                cvatom[i][5] += ric1 * Fi[2];
                cvatom[i][6] += ric1 * Fi[0];
                cvatom[i][7] += ric2 * Fi[0];
                cvatom[i][8] += ric2 * Fi[1];

                const double rjc0 = THIRD * (del_rji[0] - del_rkj[0]);
                const double rjc1 = THIRD * (del_rji[1] - del_rkj[1]);
                const double rjc2 = THIRD * (del_rji[2] - del_rkj[2]);
                cvatom[j][0] += rjc0 * Fj[0];
                cvatom[j][1] += rjc1 * Fj[1];
                cvatom[j][2] += rjc2 * Fj[2];
                cvatom[j][3] += rjc0 * Fj[1];
                cvatom[j][4] += rjc0 * Fj[2];
                cvatom[j][5] += rjc1 * Fj[2];
                cvatom[j][6] += rjc1 * Fj[0];
                cvatom[j][7] += rjc2 * Fj[0];
                cvatom[j][8] += rjc2 * Fj[1];

                const double rkc0 = THIRD * (del_rki[0] + del_rkj[0]);
                const double rkc1 = THIRD * (del_rki[1] + del_rkj[1]);
                const double rkc2 = THIRD * (del_rki[2] + del_rkj[2]);
                cvatom[k][0] += rkc0 * Fk[0];
                cvatom[k][1] += rkc1 * Fk[1];
                cvatom[k][2] += rkc2 * Fk[2];
                cvatom[k][3] += rkc0 * Fk[1];
                cvatom[k][4] += rkc0 * Fk[2];
                cvatom[k][5] += rkc1 * Fk[2];
                cvatom[k][6] += rkc1 * Fk[0];
                cvatom[k][7] += rkc2 * Fk[0];
                cvatom[k][8] += rkc2 * Fk[1];
              }
            }
          }
        }
      }
    }
  }
  if (vflag_fdotr) virial_fdotr_compute();
}


double PairUF3::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
                       double /*factor_coul*/, double factor_lj, double &fforce)
{
  double value = 0.0;
  const double r = sqrt(rsq);

  if (rsq < cutsq[itype][jtype]) {
    const int knot_start_index = uf3_potential->get_starting_index_2b(itype, jtype, r);
    double **cached_constants_2b = uf3_potential->cached_constants_2b[itype][jtype];
    double **cached_constants_2b_deri = uf3_potential->cached_constants_2b_deri[itype][jtype];

    double force_2b = cached_constants_2b_deri[knot_start_index - 1][0];
    force_2b += r   * cached_constants_2b_deri[knot_start_index - 1][1];
    force_2b += rsq * cached_constants_2b_deri[knot_start_index - 1][2];
    force_2b +=       cached_constants_2b_deri[knot_start_index - 2][3];
    force_2b += r   * cached_constants_2b_deri[knot_start_index - 2][4];
    force_2b += rsq * cached_constants_2b_deri[knot_start_index - 2][5];
    force_2b +=       cached_constants_2b_deri[knot_start_index - 3][6];
    force_2b += r   * cached_constants_2b_deri[knot_start_index - 3][7];
    force_2b += rsq * cached_constants_2b_deri[knot_start_index - 3][8];
    fforce = factor_lj * force_2b;

    const double rth = rsq * r;
    value =        cached_constants_2b[knot_start_index    ][0];
    value += r   * cached_constants_2b[knot_start_index    ][1];
    value += rsq * cached_constants_2b[knot_start_index    ][2];
    value += rth * cached_constants_2b[knot_start_index    ][3];
    value +=       cached_constants_2b[knot_start_index - 1][4];
    value += r   * cached_constants_2b[knot_start_index - 1][5];
    value += rsq * cached_constants_2b[knot_start_index - 1][6];
    value += rth * cached_constants_2b[knot_start_index - 1][7];
    value +=       cached_constants_2b[knot_start_index - 2][8];
    value += r   * cached_constants_2b[knot_start_index - 2][9];
    value += rsq * cached_constants_2b[knot_start_index - 2][10];
    value += rth * cached_constants_2b[knot_start_index - 2][11];
    value +=       cached_constants_2b[knot_start_index - 3][12];
    value += r   * cached_constants_2b[knot_start_index - 3][13];
    value += rsq * cached_constants_2b[knot_start_index - 3][14];
    value += rth * cached_constants_2b[knot_start_index - 3][15];
  }

  return factor_lj * value;
}

double PairUF3::memory_usage()
{
  double bytes = Pair::memory_usage();
  if (uf3_potential) bytes += uf3_potential->memory_usage();
  bytes += (double) maxshort * sizeof(int);    //neighshort
  bytes += (double) 1 * sizeof(int);     //maxshort
  bytes += (double) 1 * sizeof(bool);    //pot_3b
  return bytes;
}
