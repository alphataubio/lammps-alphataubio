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

#include "lammps.h"
#include "utils.h"
#include "fmt/format.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <cstdlib>
#include <string>

#include "compressed_dump_test.h"

using namespace LAMMPS_NS;

class DumpGridCompressTest : public ::testing::Test {
protected:
    void SetUp() override
    {
        compression_style  = compression_style_global;
        compression_suffix = fmt::format("{}", compression_style);
        testbinary         = "DumpGridCompressTest.bin";
        char *args[]       = {(char *)"DumpGridCompressTest", (char *)"-log", (char *)"none",
                        (char *)"-echo", (char *)"screen", (char *)"-nocite"};

        int argc = sizeof(args) / sizeof(char *);
        lmp      = new LAMMPS(argc, args, MPI_COMM_WORLD);
    }

    void TearDown() override
    {
        platform::unlink(testbinary);
        delete lmp;
    }

    void generate_text_and_compressed_dump(const std::string &text_file,
                                           const std::string &compressed_file,
                                           const std::string &dump_modify_options,
                                           const std::string &fields)
    {
        command(fmt::format("dump id1 all grid 1 {} {}", text_file, fields));
        command("dump_modify id1 units yes");
        command("run 0");
        command("undump id1");

        command(fmt::format("dump id2 all grid/{} 1 {} {}", 
                compression_style, compressed_file, fields));
        command("dump_modify id2 units yes");
        if (!dump_modify_options.empty()) {
            command(fmt::format("dump_modify id2 {}", dump_modify_options));
        }
        command("run 0");
        command("undump id2");
    }

    std::string get_dump_string(const std::string &dump_style, const std::string &file,
                                const std::string &fields)
    {
        command(fmt::format("dump id all {} 1 {} {}", dump_style, file, fields));
        command("run 0");
        command("undump id");
        std::string content = read_text_file(file);
        platform::unlink(file);
        return content;
    }

    void generate_grid_dump(const std::string &dump_style, const std::string &dump_file,
                            const std::string &fields, int ntimesteps)
    {
        command(fmt::format("dump id all {} 1 {} {}", dump_style, dump_file, fields));
        command(fmt::format("run {}", ntimesteps));
    }

    void continue_grid_dump(int ntimesteps)
    {
        command(fmt::format("run {}", ntimesteps));
    }

    void generate_compressed_dump(const std::string &dump_file, const std::string &fields,
                                  const std::string &dump_modify_options, int ntimesteps)
    {
        command(fmt::format("dump id all grid/{} 1 {} {}", compression_style, 
                dump_file, fields));
        if (!dump_modify_options.empty()) {
            command(fmt::format("dump_modify id {}", dump_modify_options));
        }
        command(fmt::format("run {}", ntimesteps));
    }

    void continue_compressed_dump(int ntimesteps)
    {
        command(fmt::format("run {}", ntimesteps));
    }

    void command(const std::string &cmd) { lmp->input->one(cmd); }

    std::string compression_style;
    std::string compression_suffix;
    std::string testbinary;
    LAMMPS *lmp;
};

TEST_F(DumpGridCompressTest, compressed_run0)
{
    if (!COMPRESS_EXECUTABLE) GTEST_SKIP();

    command("dimension 2");
    command("create_box 1 NULL 0 5 0 5 -0.5 0.5");
    command("create_atoms 1 single 0 0 0");
    command("create_atoms 1 single 1 0 0");
    command("create_atoms 1 single 0 1 0");
    command("create_atoms 1 single 2 1 0");
    command("mass * 1.0");
    command("fix ave_grid all ave/grid 1 1 1 10 10 1 vx vy vz");

    std::string text_file = "dump_grid_text_run0.grid";
    std::string compressed_file = fmt::format("dump_grid_compressed_run0.{}", 
                                            compression_suffix);
    std::string decompressed_file = compressed_file + ".decompressed";

    generate_text_and_compressed_dump(text_file, compressed_file, "", 
                                      "fix ave_grid:grid:data");

    compress_file(text_file, compressed_file + ".text");
    convert_compressed_binary_to_text(compressed_file, decompressed_file);

    ASSERT_FILE_EXISTS(text_file);
    ASSERT_FILE_EXISTS(compressed_file);
    ASSERT_FILE_EXISTS(compressed_file + ".text");
    ASSERT_FILE_EXISTS(decompressed_file);

    auto converted_content = read_text_file(decompressed_file);
    auto text_content      = read_text_file(text_file);

    ASSERT_THAT(converted_content, StrEq(text_content));

    platform::unlink(text_file);
    platform::unlink(compressed_file);
    platform::unlink(compressed_file + ".text");
    platform::unlink(decompressed_file);
}

TEST_F(DumpGridCompressTest, compressed_run1)
{
    if (!COMPRESS_EXECUTABLE) GTEST_SKIP();

    command("dimension 3");
    command("create_box 1 NULL 0 2 0 2 0 2");
    command("create_atoms 1 random 100 27632 NULL");
    command("mass * 1.0");
    command("pair_style lj/cut 2.5");
    command("pair_coeff * * 1.0 1.0");
    command("velocity all create 1.0 76287");
    command("fix ave_grid all ave/grid 1 1 1 10 10 10 vx vy vz");

    std::string text_file = "dump_grid_text_run1.grid";
    std::string compressed_file = fmt::format("dump_grid_compressed_run1.{}", 
                                            compression_suffix);
    std::string decompressed_file = compressed_file + ".decompressed";

    generate_text_and_compressed_dump(text_file, compressed_file, "", 
                                      "fix ave_grid:grid:data");

    compress_file(text_file, compressed_file + ".text");
    convert_compressed_binary_to_text(compressed_file, decompressed_file);

    ASSERT_FILE_EXISTS(text_file);
    ASSERT_FILE_EXISTS(compressed_file);
    ASSERT_FILE_EXISTS(compressed_file + ".text");
    ASSERT_FILE_EXISTS(decompressed_file);

    auto converted_content = read_text_file(decompressed_file);
    auto text_content      = read_text_file(text_file);

    ASSERT_THAT(converted_content, StrEq(text_content));

    platform::unlink(text_file);
    platform::unlink(compressed_file);
    platform::unlink(compressed_file + ".text");
    platform::unlink(decompressed_file);
}

TEST_F(DumpGridCompressTest, compressed_multi_file_run1)
{
    if (!COMPRESS_EXECUTABLE) GTEST_SKIP();

    command("dimension 3");
    command("create_box 1 NULL 0 2 0 2 0 2");
    command("create_atoms 1 random 100 27632 NULL");
    command("mass * 1.0");
    command("pair_style lj/cut 2.5");
    command("pair_coeff * * 1.0 1.0");
    command("velocity all create 1.0 76287");
    command("fix ave_grid all ave/grid 1 1 1 10 10 10 vx vy vz");

    std::string text_file = "dump_grid_text_multi_run1_*.grid";
    std::string compressed_file = fmt::format("dump_grid_compressed_multi_run1_*.{}", 
                                            compression_suffix);

    generate_text_and_compressed_dump(text_file, compressed_file, "", 
                                      "fix ave_grid:grid:data");

    std::string text_file_0 = "dump_grid_text_multi_run1_0.grid";
    std::string compressed_file_0 = fmt::format("dump_grid_compressed_multi_run1_0.{}", 
                                              compression_suffix);
    std::string decompressed_file_0 = compressed_file_0 + ".decompressed";

    ASSERT_FILE_EXISTS(text_file_0);
    ASSERT_FILE_EXISTS(compressed_file_0);

    convert_compressed_binary_to_text(compressed_file_0, decompressed_file_0);

    ASSERT_FILE_EXISTS(decompressed_file_0);

    auto converted_content = read_text_file(decompressed_file_0);
    auto text_content      = read_text_file(text_file_0);

    ASSERT_THAT(converted_content, StrEq(text_content));

    platform::unlink(text_file_0);
    platform::unlink(compressed_file_0);
    platform::unlink(decompressed_file_0);
}
