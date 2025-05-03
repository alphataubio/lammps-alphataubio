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

#include "lammps.h"
#include "utils.h"
#include "fmt/format.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <string>
#include <mpi.h>

using namespace LAMMPS_NS;

class DumpGridTest : public ::testing::Test {
protected:
    void SetUp() override
    {
        testbinary = "DumpGridTest_tmp.bin";
        char *args[] = {(char *)"DumpGridTest", (char *)"-log", (char *)"none",
                        (char *)"-echo", (char *)"screen", (char *)"-nocite"};

        int argc = sizeof(args) / sizeof(char *);
        lmp = new LAMMPS(argc, args, MPI_COMM_WORLD);
    }

    void TearDown() override
    {
        platform::unlink(testbinary);
        delete lmp;
    }

    void command(const std::string &cmd) { lmp->input->one(cmd); }

    std::string testbinary;
    LAMMPS *lmp;
};

TEST_F(DumpGridTest, dump_grid_run0)
{
    command("dimension 2");
    command("create_box 1 NULL 0 5 0 5 -0.5 0.5");
    command("create_atoms 1 single 0 0 0");
    command("create_atoms 1 single 1 0 0");
    command("create_atoms 1 single 0 1 0");
    command("create_atoms 1 single 2 1 0");
    command("create_atoms 1 single 1 2 0");
    command("create_atoms 1 single 2 2 0");
    command("create_atoms 1 single 3 2 0");
    command("create_atoms 1 single 2 3 0");
    command("create_atoms 1 single 3 3 0");
    command("create_atoms 1 single 4 3 0");
    command("create_atoms 1 single 3 4 0");
    command("create_atoms 1 single 4 4 0");
    command("create_atoms 1 single 5 4 0");
    command("pair_style zero 3.5");
    command("pair_coeff * *");
    command("pair_modify shift yes");
    command("mass * 1.0");
    command("fix ave_grid all ave/grid 1 1 1 10 10 1 vx vy vz");
    command("run 0 post no");

    int natoms = lmp->atom->natoms;
    ASSERT_EQ(natoms, 13);

    command("write_dump all grid test_run0.grid.text fix ave_grid:grid:data");

    TearDown();
    platform::unlink("test_run0.grid.text");
    SetUp();

    command("boundary f f f");
    command("dimension 2");
    command("read_dump test_run0.grid.text 0 x y z box no format grid fix ave_grid:grid:data");
    command("fix ave_grid all ave/grid 1 1 1 10 10 1 vx vy vz");

    natoms = lmp->atom->natoms;
    ASSERT_EQ(natoms, 13);

    platform::unlink("test_run0.grid.text");
}

TEST_F(DumpGridTest, run1plus1)
{
    command("dimension 3");
    command("create_box 1 NULL 0 2 0 2 0 2");
    command("create_atoms 1 random 100 27632 NULL");
    command("mass * 1.0");
    command("pair_style lj/cut 2.5");
    command("pair_coeff * * 1.0 1.0");
    command("velocity all create 1.0 76287");
    command("compute 1 all stress/atom NULL");
    command("fix ave_grid all ave/grid 1 1 1 10 10 10 vx vy vz");
    command("run 1 post no");

    command("write_dump all grid test_run1.grid.text fix ave_grid:grid:data");

    int natoms = lmp->atom->natoms;

    TearDown();
    platform::unlink("test_run1.grid.text");
    SetUp();

    command("dimension 3");
    command("read_dump test_run1.grid.text 1 grid gdim box no format grid fix ave_grid:grid:data");
    command("fix ave_grid all ave/grid 1 1 1 10 10 10 vx vy vz");

    int natoms2 = lmp->atom->natoms;
    ASSERT_EQ(natoms, natoms2);

    // Continue from last timestep
    command("reset_timestep 1");
    command("mass * 1.0");
    command("pair_style lj/cut 2.5");
    command("pair_coeff * * 1.0 1.0");
    command("run 1 post no");

    command("write_dump all grid test_run2.grid.text fix ave_grid:grid:data");

    platform::unlink("test_run1.grid.text");
    platform::unlink("test_run2.grid.text");
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
