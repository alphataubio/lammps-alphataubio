/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(write_psf,WritePsf);
// clang-format on
#else

#ifndef LMP_WRITE_PSF_H
#define LMP_WRITE_PSF_H

#include "command.h"

namespace LAMMPS_NS {

class WritePsf : public Command {
 public:
  WritePsf(class LAMMPS *);
  void command(int, char **) override;
  void write(const std::string &);

 private:
  int me, nprocs;
  int igroup, groupbit;    // group that WritePsf is performed on
  FILE *fp;
  bigint natoms_local, natoms;
  bigint nbonds_local, nbonds;
  bigint nangles_local, nangles;
  bigint ndihedrals_local, ndihedrals;
  bigint nimpropers_local, nimpropers;
  
  struct psf_atom {
    tagint tag, molecule;
    int type;
    double q;
    char segment[9], residue[9], name[9];
  };
  std::vector<psf_atom> psf_atoms;

  void header();
  void atoms();
  void bonds();
  void angles();
  void dihedrals();
  void impropers();
  int count();
  int pack_bond(tagint **);
  int pack_angle(tagint **);
  //int pack_dihedral(tagint **);
  //int pack_improper(tagint **);

};


}    // namespace LAMMPS_NS

#endif
#endif
