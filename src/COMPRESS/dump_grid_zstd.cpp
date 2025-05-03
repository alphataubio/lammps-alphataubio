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
   Contributing author: Mitch Murphy (alphataubio at gmail)
------------------------------------------------------------------------- */

#ifdef LAMMPS_ZSTD

#include "dump_grid_zstd.h"

#include "domain.h"
#include "error.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

DumpGridZStd::DumpGridZStd(LAMMPS *lmp, int narg, char **arg) : DumpGrid(lmp, narg, arg), writer(lmp)
{
  if (!compressed) error->all(FLERR,"Dump grid/zstd only writes compressed files");
}

/* ---------------------------------------------------------------------- */

void DumpGridZStd::openfile()
{
  if (singlefile_opened) return;
  if (multifile == 0) singlefile_opened = 1;

  // if one file per timestep, replace '*' with current timestep

  char *filecurrent = writer.get_file_name(filename, update->ntimestep, padflag);

  if (append_flag && access(filecurrent, F_OK) == 0) {
    error->all(FLERR, "Dump grid/zstd file {} already exists and cannot be appended to", filecurrent);
  }

  writer.open(filecurrent, append_flag);

  delete[] filecurrent;
}

/* ---------------------------------------------------------------------- */

void DumpGridZStd::write_header(bigint ndump)
{
  std::string header = DumpGrid::get_header(ndump);
  writer.write_header(header);
}

/* ---------------------------------------------------------------------- */

void DumpGridZStd::write_data(int n, double *mybuf)
{
  if (buffer_flag) {
    writer.write(n, mybuf, format_line_user,
                 format_int_user, format_bigint_user);
  } else if (format_column_user) {
    writer.write(n, mybuf, format_column_user,
                 format_int_user, format_bigint_user);
  } else {
    writer.write(n, mybuf, format_column_default);
  }
}

/* ---------------------------------------------------------------------- */

void DumpGridZStd::write()
{
  DumpGrid::write();
  if (filewriter) {
    if (fp) {
      writer.close();
      fp = nullptr;
    }
  }
}

/* ---------------------------------------------------------------------- */

int DumpGridZStd::modify_param(int narg, char **arg)
{
  int consumed = DumpGrid::modify_param(narg, arg);
  if (consumed == 0) {
    if (strcmp(arg[0], "compression_level") == 0) {
      if (narg < 2) error->all(FLERR, "Illegal dump_modify command");
      int min_level = -7;
      int max_level = 22;
      int compression_level = utils::inumeric(FLERR, arg[1], false, lmp);
      if ((compression_level < min_level) || (compression_level > max_level))
        error->all(FLERR, "Illegal dump_modify command");
      writer.compression_level = compression_level;
      consumed = 2;
    } else if (strcmp(arg[0], "checksum") == 0) {
      if (narg < 2) error->all(FLERR, "Illegal dump_modify command");
      writer.checksum = utils::logical(FLERR, arg[1], false, lmp);
      consumed = 2;
    }
  }
  return consumed;
}

#endif
