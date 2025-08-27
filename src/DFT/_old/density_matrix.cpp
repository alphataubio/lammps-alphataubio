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

#include "pair_dft.h"
#include <Eigen/Dense>
#include <cmath>
#include <algorithm>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

DensityMatrix::DensityMatrix(int nbasis) : n_basis(nbasis)
{
  current.setZero(n_basis, n_basis);
  previous.setZero(n_basis, n_basis);
  difference.setZero(n_basis, n_basis);
  
  use_diis = true;
  diis_size = 6;
  
  diis_fock.clear();
  diis_error.clear();
}

/* ---------------------------------------------------------------------- */

DensityMatrix::~DensityMatrix()
{
}

/* ---------------------------------------------------------------------- */

void DensityMatrix::initialize_guess(const Eigen::MatrixXd &S)
{
  // Simple guess: set diagonal elements proportional to overlap
  current.setZero();
  
  for (int i = 0; i < n_basis; i++) {
    current(i, i) = 1.0 / sqrt(S(i, i));
  }
  
  // Normalize
  double trace = (current * S).trace();
  if (trace > 0) {
    current /= trace;
  }
  
  previous = current;
}

/* ---------------------------------------------------------------------- */

void DensityMatrix::update_from_mo(const Eigen::MatrixXd &C, int nocc)
{
  // Save previous density
  previous = current;
  
  // Build new density matrix from occupied orbitals
  // D = 2 * C_occ * C_occ^T for closed shell
  current.setZero();
  
  for (int i = 0; i < nocc; i++) {
    current += 2.0 * C.col(i) * C.col(i).transpose();
  }
  
  // Compute difference for convergence check
  difference = current - previous;
}

/* ---------------------------------------------------------------------- */

void DensityMatrix::mix_with_previous(double mixing_param)
{
  // Simple linear mixing for stability
  // D_new = (1-alpha)*D_old + alpha*D_current
  current = mixing_param * current + (1.0 - mixing_param) * previous;
}

/* ---------------------------------------------------------------------- */

double DensityMatrix::get_change() const
{
  // RMS change in density matrix
  double sum = 0.0;
  int count = 0;
  
  for (int i = 0; i < n_basis; i++) {
    for (int j = 0; j < n_basis; j++) {
      sum += difference(i, j) * difference(i, j);
      count++;
    }
  }
  
  return sqrt(sum / count);
}

/* ---------------------------------------------------------------------- */

std::vector<double> DensityMatrix::compute_mulliken_charges(const Eigen::MatrixXd &S)
{
  // Mulliken population analysis
  // q_A = Z_A - sum_mu(P_mu,mu * S_mu,mu) for mu on atom A
  
  // For now, return empty vector (need atom mapping)
  std::vector<double> charges;
  
  // TODO: Implement with proper basis-to-atom mapping
  
  return charges;
}

/* ---------------------------------------------------------------------- */

std::vector<double> DensityMatrix::compute_lowdin_charges(const Eigen::MatrixXd &S)
{
  // Löwdin population analysis
  // Uses S^(1/2) transformation
  
  // Compute S^(1/2)
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S);
  Eigen::MatrixXd S_sqrt = es.operatorSqrt();
  
  // Transform density matrix
  Eigen::MatrixXd P_lowdin = S_sqrt * current * S_sqrt;
  
  // For now, return empty vector (need atom mapping)
  std::vector<double> charges;
  
  // TODO: Implement with proper basis-to-atom mapping
  
  return charges;
}

/* ---------------------------------------------------------------------- */

void DensityMatrix::apply_diis(Eigen::MatrixXd &F)
{
  if (!use_diis) return;
  
  // DIIS (Direct Inversion of Iterative Subspace)
  // Accelerates SCF convergence
  
  // Store current Fock matrix and error
  diis_fock.push_back(F);
  
  // Compute error matrix: e = FDS - SDF
  // For simplicity, using gradient of energy wrt density
  Eigen::MatrixXd error = F * current - current * F;
  diis_error.push_back(error);
  
  // Keep only last diis_size iterations
  if (diis_fock.size() > static_cast<size_t>(diis_size)) {
    diis_fock.erase(diis_fock.begin());
    diis_error.erase(diis_error.begin());
  }
  
  int n = diis_fock.size();
  if (n < 2) return;  // Need at least 2 iterations
  
  // Build B matrix
  Eigen::MatrixXd B(n + 1, n + 1);
  B.setZero();
  
  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= i; j++) {
      double val = (diis_error[i].cwiseProduct(diis_error[j])).sum();
      B(i, j) = val;
      B(j, i) = val;
    }
    B(i, n) = -1.0;
    B(n, i) = -1.0;
  }
  B(n, n) = 0.0;
  
  // Solve for coefficients
  Eigen::VectorXd rhs(n + 1);
  rhs.setZero();
  rhs(n) = -1.0;
  
  Eigen::VectorXd c = B.colPivHouseholderQr().solve(rhs);
  
  // Build extrapolated Fock matrix
  F.setZero();
  for (int i = 0; i < n; i++) {
    F += c(i) * diis_fock[i];
  }
}
