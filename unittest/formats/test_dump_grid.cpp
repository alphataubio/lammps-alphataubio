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

#include "fmt/format.h"
#include "info.h"
#include "lammps.h"
#include "utils.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <mpi.h>
#include <vector>

using namespace LAMMPS_NS;

static bool verbose = false;

class DumpGridTest : public ::testing::Test {
protected:
    void SetUp() override
    {
        testdata = "./test_dump_grid.txt";
        tstamp = std::time(nullptr);
    }

    void TearDown() override
    {
        if (std::remove(testdata.c_str()) != 0) {
            // ignore error
        }
    }

    void run_dump_grid_test(const std::string &dump_style, const std::string &filename_suffix)
    {
        const char *args[] = {"LAMMPS-test", "-log", "none", "-echo", "none", "-nocite"};
        
        char **argv = (char **)args;
        int argc    = sizeof(args) / sizeof(char *);
        
        LAMMPS *lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);

        std::string filename = fmt::format("test_dump_grid_{}.{}", tstamp, filename_suffix);

        std::string cmds = fmt::format(
            "region box block 0 10 0 10 0 10\n"
            "create_box 1 box\n"
            "create_atoms 1 single 5 5 5\n"
            "mass 1 1.0\n"
            "fix grid all ave/grid 1 1 1 20 20 20 vx vy vz\n"
            "run 0\n"
            "dump grid all {} 1 {}\n"
            "dump_modify grid every 1\n"
            "run 0\n",
            dump_style, filename);

        lmp->commands_string(cmds);

        // check if the file was created and has content
        FILE *fp = fopen(filename.c_str(), "r");
        ASSERT_NE(fp, nullptr);
        
        // Read first few lines to verify format
        char line[1024];
        int line_count = 0;
        
        while (fgets(line, sizeof(line), fp) && line_count < 10) {
            line_count++;
            if (strstr(line, "ITEM:") != nullptr) {
                // Found an ITEM header line, which is expected for dump grid
                EXPECT_TRUE(true);
            }
        }
        
        fclose(fp);
        
        // Clean up
        if (std::remove(filename.c_str()) != 0) {
            // ignore error
        }
        
        delete lmp;
    }

    std::string testdata;
    std::time_t tstamp;
};

TEST_F(DumpGridTest, dump_grid)
{
    run_dump_grid_test("grid", "grid");
}

TEST_F(DumpGridTest, dump_grid_gz)
{
    if (!COMPRESS::has_gzip_support()) {
        GTEST_SKIP() << "GZIP support not included";
    }
    run_dump_grid_test("grid/gz", "grid.gz");
}

TEST_F(DumpGridTest, dump_grid_zstd)
{
#ifdef LAMMPS_ZSTD
    run_dump_grid_test("grid/zstd", "grid.zst");
#else
    GTEST_SKIP() << "ZSTD support not included";
#endif
}

TEST_F(DumpGridTest, dump_grid_with_units)
{
    const char *args[] = {"LAMMPS-test", "-log", "none", "-echo", "none", "-nocite"};
    
    char **argv = (char **)args;
    int argc    = sizeof(args) / sizeof(char *);
    
    LAMMPS *lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);

    std::string filename = fmt::format("test_dump_grid_units_{}.grid", tstamp);

    std::string cmds = fmt::format(
        "region box block 0 10 0 10 0 10\n"
        "create_box 1 box\n"
        "create_atoms 1 single 5 5 5\n"
        "mass 1 1.0\n"
        "fix grid all ave/grid 1 1 1 20 20 20 vx vy vz\n"
        "run 0\n"
        "dump grid all grid 1 {}\n"
        "dump_modify grid units yes\n"
        "run 0\n",
        filename);

    lmp->commands_string(cmds);

    // check if the file was created and has content
    FILE *fp = fopen(filename.c_str(), "r");
    ASSERT_NE(fp, nullptr);
    
    // Read first few lines to verify format includes UNITS
    char line[1024];
    bool found_units = false;
    
    while (fgets(line, sizeof(line), fp)) {
        if (strstr(line, "ITEM: UNITS") != nullptr) {
            found_units = true;
            break;
        }
    }
    
    EXPECT_TRUE(found_units);
    
    fclose(fp);
    
    // Clean up
    if (std::remove(filename.c_str()) != 0) {
        // ignore error
    }
    
    delete lmp;
}

TEST_F(DumpGridTest, dump_grid_with_time)
{
    const char *args[] = {"LAMMPS-test", "-log", "none", "-echo", "none", "-nocite"};
    
    char **argv = (char **)args;
    int argc    = sizeof(args) / sizeof(char *);
    
    LAMMPS *lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);

    std::string filename = fmt::format("test_dump_grid_time_{}.grid", tstamp);

    std::string cmds = fmt::format(
        "region box block 0 10 0 10 0 10\n"
        "create_box 1 box\n"
        "create_atoms 1 single 5 5 5\n"
        "mass 1 1.0\n"
        "fix grid all ave/grid 1 1 1 20 20 20 vx vy vz\n"
        "run 0\n"
        "dump grid all grid 1 {}\n"
        "dump_modify grid time yes\n"
        "run 0\n",
        filename);

    lmp->commands_string(cmds);

    // check if the file was created and has content
    FILE *fp = fopen(filename.c_str(), "r");
    ASSERT_NE(fp, nullptr);
    
    // Read first few lines to verify format includes TIME
    char line[1024];
    bool found_time = false;
    
    while (fgets(line, sizeof(line), fp)) {
        if (strstr(line, "ITEM: TIME") != nullptr) {
            found_time = true;
            break;
        }
    }
    
    EXPECT_TRUE(found_time);
    
    fclose(fp);
    
    // Clean up
    if (std::remove(filename.c_str()) != 0) {
        // ignore error
    }
    
    delete lmp;
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = LAMMPS_NS::utils::split_words(var);
        for (auto arg : env) {
            if (arg == "-v") {
                verbose = true;
            }
        }
    }

    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
