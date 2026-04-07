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
---------------------------------------------------------------------- */

#include "uf3_potential.h"

#include "uf3_bspline_basis2.h"
#include "uf3_bspline_basis3.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "platform.h"
#include "text_file_reader.h"

#include <algorithm>
#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

UF3Potential::UF3Potential(LAMMPS *lmp, const std::string &potf_name,
                           double **cutsq_, int **setflag_, std::vector<std::string> &elements_, bool pot_3b_) :
    Pointers(lmp),
    cutsq(cutsq_), setflag(setflag_), elements(elements_), element_to_type(), pot_3b(pot_3b_),
    setflag_3b(nullptr), knot_spacing_type_2b(nullptr), knot_spacing_type_3b(nullptr),
    cut_2b(nullptr), cut_3b(nullptr), cut_3b_list(nullptr), min_cut_3b(nullptr),
    knot_spacing_2b(nullptr), knot_spacing_3b(nullptr), n2b_knots_array(nullptr),
    n2b_coeff_array(nullptr), n2b_knots_array_size(nullptr), n2b_coeff_array_size(nullptr),
    cached_constants_2b(nullptr), cached_constants_2b_deri(nullptr), map_3b(nullptr),
    n3b_knots_array(nullptr), n3b_coeff_array(nullptr), n3b_knots_array_size(nullptr),
    n3b_coeff_array_size(nullptr), coeff_for_der_jk(nullptr), coeff_for_der_ik(nullptr),
    coeff_for_der_ij(nullptr), cached_constants_3b(nullptr), cached_constants_3b_deri(nullptr),
    get_starting_index_2b_ptr(nullptr), get_starting_index_3b_ptr(nullptr)
{

  allocated = 0;
  max_num_knots_2b = 0;
  max_num_coeff_2b = 0;
  max_num_knots_3b = 0;
  max_num_coeff_3b = 0;
  tot_interaction_count_3b = 0;

  for (int i=1; i < elements.size(); i++) element_to_type[elements[i]] = i;

  if (comm->me == 0) utils::logmesg(lmp, "Reading UF3 potential {}... ", potf_name);
  double time1 = platform::walltime();
  allocate();
  uf3_read_unified_pot_file(potf_name);
  communicate();
  create_bsplines();
  if (comm->me == 0) utils::logmesg(lmp, "done ({:.2f} seconds).", platform::walltime() - time1);

}

/* ---------------------------------------------------------------------- */

UF3Potential::~UF3Potential()
{
  if (allocated) {
    memory->destroy(cut_2b);
    memory->destroy(knot_spacing_type_2b);
    memory->destroy(knot_spacing_2b);
    memory->destroy(n2b_knots_array_size);
    memory->destroy(n2b_coeff_array_size);
    memory->destroy(n2b_knots_array);
    memory->destroy(n2b_coeff_array);
    memory->destroy(cached_constants_2b);
    memory->destroy(cached_constants_2b_deri);

    if (pot_3b) {
      memory->destroy(setflag_3b);
      memory->destroy(cut_3b);
      memory->destroy(cut_3b_list);
      memory->destroy(min_cut_3b);
      memory->destroy(knot_spacing_type_3b);
      memory->destroy(knot_spacing_3b);
      memory->destroy(map_3b);
      memory->destroy(n3b_knots_array_size);
      memory->destroy(n3b_coeff_array_size);
      memory->destroy(n3b_knots_array);
      memory->destroy(n3b_coeff_array);
      memory->destroy(coeff_for_der_jk);
      memory->destroy(coeff_for_der_ik);
      memory->destroy(coeff_for_der_ij);
      memory->destroy(cached_constants_3b);
      memory->destroy(cached_constants_3b_deri);
    }
  }
}



void UF3Potential::allocate()
{
  allocated = 1;
  const int np1 = atom->ntypes + 1;
  // cut_2b is specific to this pair style. We will set the values in cut_2b
  memory->create(cut_2b, np1, np1, "uf3:cut_2b");
  // Contains info about type of knot_spacing--> 0 = uniform knot spacing (default)
  // 1 = non-uniform knot spacing
  memory->create(knot_spacing_type_2b, np1, np1, "uf3:knot_spacing_type_2b");
  memory->create(knot_spacing_2b, np1, np1, "uf3:knot_spacing_2b");

  // Contains size of 2b knots vectors and 2b coeff matrices
  memory->create(n2b_knots_array_size, np1, np1, "uf3:n2b_knots_array_size");
  memory->create(n2b_coeff_array_size, np1, np1, "uf3:n2b_coeff_array_size");

  if (pot_3b) {
    // Contains info about wether UF potential were found for type i, j and k
    memory->create(setflag_3b, np1, np1, np1, "uf3:setflag_3b");
    // Contains info about 3-body cutoff distance for type i, j and k
    memory->create(cut_3b, np1, np1, np1, "uf3:cut_3b");
    // Contains info about 3-body cutoff distance for type i, j and k
    // for constructing 3-body list
    memory->create(cut_3b_list, np1, np1, "uf3:cut_3b_list");
    // Contains info about minimum 3-body cutoff distance for type i, j and k
    memory->create(min_cut_3b, np1, np1, np1, 3, "uf3:min_cut_3b");
    // Contains info about type of knot_spacing--> 0 = uniform knot spacing (default)
    // 1 = non-uniform knot spacing
    memory->create(knot_spacing_type_3b, np1, np1, np1, "uf3:knot_spacing_type_3b");
    memory->create(knot_spacing_3b, np1, np1, np1, 3, "uf3:knot_spacing_3b");

    tot_interaction_count_3b = 0;
    // conatins map of I-J-K interaction
    memory->create(map_3b, np1, np1, np1, "uf3:map_3b");

    // setting cut_3b, setflag = 0 and map_3b
    for (int i = 1; i < np1; i++) {
      for (int j = 1; j < np1; j++) {
        cut_3b_list[i][j] = setflag[i][j] = 0;
        n2b_coeff_array_size[i][j] = n2b_knots_array_size[i][j] = 0;
        for (int k = 1; k < np1; k++) {
          cut_3b[i][j][k] = setflag_3b[i][j][k] = 0;
          min_cut_3b[i][j][k][0] = min_cut_3b[i][j][k][1] = min_cut_3b[i][j][k][2] = 0;
          map_3b[i][j][k] = tot_interaction_count_3b++;
        }
      }
    }

    // contains sizes of 3b knots vectors and 3b coeff matrices
    memory->create(n3b_knots_array_size, tot_interaction_count_3b, 3, "uf3:n3b_knots_array_size");
    memory->create(n3b_coeff_array_size, tot_interaction_count_3b, 3, "uf3:n3b_coeff_array_size");
    for (int i = 0; i < tot_interaction_count_3b; i++) {
      n3b_coeff_array_size[i][0] = n3b_coeff_array_size[i][1] = n3b_coeff_array_size[i][2] = 0;
      n3b_knots_array_size[i][0] = n3b_knots_array_size[i][1] = n3b_knots_array_size[i][2] = 0;
    }

  }
}

void UF3Potential::uf3_read_unified_pot_file(const std::string &potf_name)
{
  //Go through the entire file and get the sizes of knot vectors and
  //coeff vectors/matrices
  //
  //Create arrays
  //
  //Go through the file again and read the knots and coefficients
  //

  const int np1 = atom->ntypes + 1;

  FILE *fp = utils::open_potential(potf_name, lmp, nullptr);
  if (!fp)
    error->all(FLERR, "Cannot open UF3 potential file {}: {}", potf_name, utils::getsyserror());

  TextFileReader txtfilereader(fp, "UF3:POTFP");
  txtfilereader.ignore_comments = false;

  //while loop over the entire file, find blocks starting with #UF3 POT
  //if block found read the very next line to determine 2B or 3B block
  //if 2B read the knot vector and coeff vector size
  //if 3B read the knot vectors and coeff matrix size
  int line_counter = 1;
  char *line;
  while ((line = txtfilereader.next_line(1))) {
    Tokenizer line_token(line);

    // Detect start of a block
    if (line_token.contains("#UF3 POT")) {
      // Block start detected
      if (line_token.contains("UNITS:") == 0)
        error->all(FLERR, "UF3: {} does not contain UNITS metadata in header", potf_name);

      // Read the 2nd line of the block
      std::string temp_line = txtfilereader.next_line(1);
      line_counter++;
      ValueTokenizer fp2nd_line(temp_line);

      std::string nbody_on_file = fp2nd_line.next_string();
      if (nbody_on_file == "2B") {
        //2B block
        if (fp2nd_line.count() != 6)
          error->all(FLERR, "UF3: Expected 6 words on line {} of {} file but found {} word/s",
                     line_counter, potf_name, fp2nd_line.count());

        //get the elements
        const std::string element1 = fp2nd_line.next_string();
        const std::string element2 = fp2nd_line.next_string();
        const int itype = element_to_type[element1];
        const int jtype = element_to_type[element2];

        if ((itype != 0) && (jtype != 0)) {
          //Trailing and leading trim check
          int leading_trim = fp2nd_line.next_int();
          int trailing_trim = fp2nd_line.next_int();
          if (leading_trim != 0) error->all(FLERR, "UF3 implemented only for leading_trim=0");
          if (trailing_trim != 3) error->all(FLERR, "UF3 implemented only for trailing_trim=3");

          //read next line, should contain cutoff and size of knot vector
          temp_line = txtfilereader.next_line(1);
          line_counter++;
          ValueTokenizer fp3rd_line(temp_line);
          if (fp3rd_line.count() != 2)
            error->all(FLERR,
                       "UF3: Expected only 2 words on 3rd line => "
                       "Rij_CUTOFF NUM_OF_KNOTS. Found {} word/s",
                       fp3rd_line.count());

          //cut is used in init_one which is called by pair.cpp at line 267
          //where the return of init_one is squared
          cut_2b[itype][jtype] = fp3rd_line.next_double();
          cut_2b[jtype][itype] = cut_2b[itype][jtype];

          int num_knots_2b = fp3rd_line.next_int();
          n2b_knots_array_size[itype][jtype] = num_knots_2b;
          n2b_knots_array_size[jtype][itype] = num_knots_2b;
          max_num_knots_2b = std::max(max_num_knots_2b, num_knots_2b);

          //skip next line
          txtfilereader.skip_line();
          line_counter++;

          //read number of coeff
          temp_line = txtfilereader.next_line(1);
          line_counter++;
          ValueTokenizer fp5th_line(temp_line);

          int num_coeff_2b = fp5th_line.next_int();
          if (num_coeff_2b <= 0)
            error->all(FLERR,
                       "UF3: 0 or negative number found for num_coeff_2b"
                       " on line {} of the potential file",
                       line_counter);
          n2b_coeff_array_size[itype][jtype] = num_coeff_2b;
          n2b_coeff_array_size[jtype][itype] = num_coeff_2b;
          max_num_coeff_2b = std::max(max_num_coeff_2b, num_coeff_2b);
        } //FIXME no else clause in case element not found
      } else if ((nbody_on_file == "3B") && (pot_3b)) {
        //3B block
        if (fp2nd_line.count() != 7)
          error->all(FLERR,
                     "UF3: Expected 7 words on line {} of {} file"
                     "but found {} word/s",
                     line_counter, potf_name, fp2nd_line.count());

        if (nbody_on_file == "3B") {
          //get the elements
          const std::string element1 = fp2nd_line.next_string();
          const std::string element2 = fp2nd_line.next_string();
          const std::string element3 = fp2nd_line.next_string();
          const int itype = element_to_type[element1];
          const int jtype = element_to_type[element2];
          const int ktype = element_to_type[element3];


          if ((itype != 0) && (jtype != 0) && (ktype != 0)) {
            //Trailing and leading trim check
            int leading_trim = fp2nd_line.next_int();
            int trailing_trim = fp2nd_line.next_int();
            if (leading_trim != 0) error->all(FLERR, "UF3 implemented only for leading_trim=0");
            if (trailing_trim != 3) error->all(FLERR, "UF3 implemented only for trailing_trim=3");

            //read next line, should contain cutoffs and size of knot vectors
            temp_line = txtfilereader.next_line(6);
            line_counter++;
            ValueTokenizer fp3rd_line(temp_line);

            if (fp3rd_line.count() != 6)
              error->all(FLERR,
                         "UF3: Expected only 6 numbers on 3rd line => "
                         "Rjk_CUTOFF Rik_CUTOFF Rij_CUTOFF NUM_OF_KNOTS_JK "
                         "NUM_OF_KNOTS_IK NUM_OF_KNOTS_IJ Found {} number/s",
                         fp3rd_line.count());

            double cut3b_rjk = fp3rd_line.next_double();
            double cut3b_rij = fp3rd_line.next_double();
            double cut3b_rik = fp3rd_line.next_double();

            if (cut3b_rij != cut3b_rik)
              error->all(FLERR,
                         "UF3: rij!=rik for {}-{}-{}. "
                         "Current implementation only works for rij=rik",
                         element1, element2, element3);

            if (2 * cut3b_rik != cut3b_rjk)
              error->all(FLERR,
                         "UF3: 2rij=2rik!=rik for {}-{}-{}. "
                         "Current implementation only works for 2rij=2rik!=rik",
                         element1, element2, element3);

            cut_3b_list[itype][jtype] = std::max(cut3b_rij, cut_3b_list[itype][jtype]);
            cut_3b_list[itype][ktype] = std::max(cut_3b_list[itype][ktype], cut3b_rik);

            cut_3b[itype][jtype][ktype] = cut3b_rij;
            cut_3b[itype][ktype][jtype] = cut3b_rik;

            int num_knots_3b_jk = fp3rd_line.next_int();
            int num_knots_3b_ik = fp3rd_line.next_int();
            int num_knots_3b_ij = fp3rd_line.next_int();

            n3b_knots_array_size[map_3b[itype][jtype][ktype]][0] = num_knots_3b_jk;
            n3b_knots_array_size[map_3b[itype][jtype][ktype]][1] = num_knots_3b_ik;
            n3b_knots_array_size[map_3b[itype][jtype][ktype]][2] = num_knots_3b_ij;

            n3b_knots_array_size[map_3b[itype][ktype][jtype]][0] = num_knots_3b_jk;
            n3b_knots_array_size[map_3b[itype][ktype][jtype]][1] = num_knots_3b_ij;
            n3b_knots_array_size[map_3b[itype][ktype][jtype]][2] = num_knots_3b_ik;

            max_num_knots_3b = std::max(max_num_knots_3b, num_knots_3b_jk);
            max_num_knots_3b = std::max(max_num_knots_3b, num_knots_3b_ik);
            max_num_knots_3b = std::max(max_num_knots_3b, num_knots_3b_ij);

            //skip next 3 line
            txtfilereader.skip_line();
            txtfilereader.skip_line();
            txtfilereader.skip_line();
            line_counter += 3;

            //read number of coeff
            temp_line = txtfilereader.next_line(3);
            line_counter++;
            ValueTokenizer fp7th_line(temp_line);

            if (fp7th_line.count() != 3)
              error->all(FLERR,
                         "UF3: Expected 3 numbers on 7th line => "
                         "SHAPE_OF_COEFF_MATRIX[I][J][K] found {} numbers",
                         fp7th_line.count());

            int coeff_matrix_dim1 = fp7th_line.next_int();
            int coeff_matrix_dim2 = fp7th_line.next_int();
            int coeff_matrix_dim3 = fp7th_line.next_int();

            n3b_coeff_array_size[map_3b[itype][jtype][ktype]][0] = coeff_matrix_dim1;
            n3b_coeff_array_size[map_3b[itype][jtype][ktype]][1] = coeff_matrix_dim2;
            n3b_coeff_array_size[map_3b[itype][jtype][ktype]][2] = coeff_matrix_dim3;

            n3b_coeff_array_size[map_3b[itype][ktype][jtype]][0] = coeff_matrix_dim2;
            n3b_coeff_array_size[map_3b[itype][ktype][jtype]][1] = coeff_matrix_dim1;
            n3b_coeff_array_size[map_3b[itype][ktype][jtype]][2] = coeff_matrix_dim3;

            max_num_coeff_3b = std::max(max_num_coeff_3b, coeff_matrix_dim1);
            max_num_coeff_3b = std::max(max_num_coeff_3b, coeff_matrix_dim2);
            max_num_coeff_3b = std::max(max_num_coeff_3b, coeff_matrix_dim3);
          }
        }
      } else {
        if (!((nbody_on_file == "3B") && (!pot_3b)))
          error->all(FLERR, "UF3: Expected '2B' or '3B' on line {} of {}", line_counter, potf_name);
      }
    }    //if of #UF3 POT
    line_counter++;
  }    // while

  //Create knot and coeff arrays
  if (max_num_knots_2b <= 0)
    error->all(FLERR,
               "UF3: Error reading the size of 2B knot vector\n"
               "Possibly no 2B UF3 potential block detected in {} file",
               potf_name);
  memory->destroy(n2b_knots_array);
  memory->create(n2b_knots_array, np1, np1, max_num_knots_2b, "uf3:n2b_knots_array");

  if (max_num_coeff_2b <= 0)
    error->all(FLERR,
               "UF3: Error reading the size of 2B coeff vector\n"
               "Possibly no 2B UF3 potential block detected in {} file",
               potf_name);

  memory->destroy(n2b_coeff_array);
  memory->create(n2b_coeff_array, np1, np1, max_num_coeff_2b, "uf3:n2b_coeff_array");

  if (pot_3b) {
    if (max_num_knots_3b <= 0)
      error->all(FLERR,
                 "UF3: Error reading the size of 3B knot vector\n"
                 "Possibly no 3B UF3 potential block detected in {} file",
                 potf_name);
    memory->destroy(n3b_knots_array);
    memory->create(n3b_knots_array, tot_interaction_count_3b, 3, max_num_knots_3b,
                   "pair:n3b_knots_array");

    if (max_num_coeff_3b <= 0)
      error->all(FLERR,
                 "UF3: Error reading the size of 3B coeff matrices\n"
                 "Possibly no 3B UF3 potential block detected in {} file",
                 potf_name);
    memory->destroy(n3b_coeff_array);
    memory->create(n3b_coeff_array, tot_interaction_count_3b, max_num_coeff_3b, max_num_coeff_3b,
                   max_num_coeff_3b, "pair:n3b_coeff_array");
  }

  //Go back to the begning of the file
  txtfilereader.rewind();

  //Go through the file again and fill knot and coeff arrays
  //while loop to read the data
  while ((line = txtfilereader.next_line(1))) {
    Tokenizer line_token(line);

    //Detect start of a block
    if (line_token.contains("#UF3 POT")) {
      //Block start detected
      //Read the 2nd line of the block
      std::string temp_line = txtfilereader.next_line(1);
      ValueTokenizer fp2nd_line(temp_line);
      std::string nbody_on_file = fp2nd_line.next_string();

      if (nbody_on_file == "2B") {
        //get the elements
        const std::string element1 = fp2nd_line.next_string();
        const std::string element2 = fp2nd_line.next_string();
        const int itype = element_to_type[element1];
        const int jtype = element_to_type[element2];

        //skip the next two tokens
        fp2nd_line.skip(2);

        //uk or nk?
        std::string knot_type = fp2nd_line.next_string();
        if (knot_type == "uk") {
          knot_spacing_type_2b[itype][jtype] = knot_spacing_type_2b[jtype][itype] = 0;
        } else if (knot_type == "nk") {
          knot_spacing_type_2b[itype][jtype] = knot_spacing_type_2b[jtype][itype] = 1;
        } else
          error->all(FLERR,
                     "UF3: Expected uniform 'uk' or non-uniform 'nk' knots on line 2 of {}-{} block but found {}",
                     element1, element2, knot_type);

        if ((itype != 0) && (jtype != 0)) {
          //skip line containing info of cutoff and knot vect size
          txtfilereader.skip_line();

          int num_knots_2b = n2b_knots_array_size[itype][jtype];

          temp_line = txtfilereader.next_line(num_knots_2b);
          ValueTokenizer fp4th_line(temp_line);

          if ((int) fp4th_line.count() != num_knots_2b)
            error->all(FLERR,
                       "UF3: Error reading the 2B potential block for {}-{}\n"
                       "Expected {} numbers on 4th line of the block but found {} numbers",
                       element1, element2, num_knots_2b, fp4th_line.count());

          for (int k = 0; k < num_knots_2b; k++)
            n2b_knots_array[itype][jtype][k] = n2b_knots_array[jtype][itype][k] = fp4th_line.next_double();

          knot_spacing_2b[itype][jtype] =
              n2b_knots_array[itype][jtype][4] - n2b_knots_array[itype][jtype][3];
          knot_spacing_2b[jtype][itype] = knot_spacing_2b[itype][jtype];

          //skip next line
          txtfilereader.skip_line();

          const int num_of_coeff_2b = n2b_coeff_array_size[itype][jtype];

          temp_line = txtfilereader.next_line(num_of_coeff_2b);
          ValueTokenizer fp6th_line(temp_line);

          if ((int) fp6th_line.count() != num_of_coeff_2b)
            error->all(FLERR, "UF3: Error reading 2B block for {}-{}, expected {} numbers on line 6 but found {}",
                       element1, element2, num_of_coeff_2b, fp6th_line.count());

          for (int k = 0; k < num_of_coeff_2b; k++)
            n2b_coeff_array[itype][jtype][k] = n2b_coeff_array[jtype][itype][k] = fp6th_line.next_double();

          if (num_knots_2b != num_of_coeff_2b + 4)
            error->all(FLERR,
                       "UF3: {}-{} interaction block has incorrect knot and "
                       "coeff data nknots (={}) != ncoeffs (={}) + 3 + 1",
                       element1, element2, num_knots_2b, num_of_coeff_2b);

          setflag[itype][jtype] = setflag[jtype][itype] = 1;
        }
      }

      if ((nbody_on_file == "3B") && (pot_3b)) {
        //get the elements
        const std::string element1 = fp2nd_line.next_string();
        const std::string element2 = fp2nd_line.next_string();
        const std::string element3 = fp2nd_line.next_string();
        const int itype = element_to_type[element1];
        const int jtype = element_to_type[element2];
        const int ktype = element_to_type[element3];

        //skip the next two tokens
        fp2nd_line.skip(2);

        //uk or nk?
        std::string knot_type = fp2nd_line.next_string();
        if (knot_type == "uk") {
          knot_spacing_type_3b[itype][jtype][ktype] = knot_spacing_type_3b[itype][ktype][jtype] = 0;
        } else if (knot_type == "nk") {
          knot_spacing_type_3b[itype][jtype][ktype] = knot_spacing_type_3b[itype][ktype][jtype] = 1;
        } else
          error->all(FLERR,
                     "UF3: Expected either 'uk'(uniform-knots) or 'nk'(non-uniform knots) "
                     "Found {} on the 2nd line of {}-{}-{} interaction block",
                     knot_type, element1, element2, element3);

        if ((itype != 0) && (jtype != 0) && (ktype != 0)) {
          //skip line containing info of cutoffs and knot vector sizes
          txtfilereader.skip_line();

          const int num_knots_3b_jk = n3b_knots_array_size[map_3b[itype][jtype][ktype]][0];
          const int num_knots_3b_ik = n3b_knots_array_size[map_3b[itype][jtype][ktype]][1];
          const int num_knots_3b_ij = n3b_knots_array_size[map_3b[itype][jtype][ktype]][2];

          temp_line = txtfilereader.next_line(num_knots_3b_jk);
          ValueTokenizer fp4th_line(temp_line);
          if ((int) fp4th_line.count() != num_knots_3b_jk)
            error->all(FLERR, "UF3: Error reading 3B block {}-{}-{}, expected {} numbers on line 4 but found {} ",
                       element1, element2, element3, num_knots_3b_jk, fp4th_line.count());

          for (int i = 0; i < num_knots_3b_jk; i++) {
            n3b_knots_array[map_3b[itype][jtype][ktype]][0][i] = fp4th_line.next_double();
            n3b_knots_array[map_3b[itype][ktype][jtype]][0][i] =
                n3b_knots_array[map_3b[itype][jtype][ktype]][0][i];
          }

          min_cut_3b[itype][jtype][ktype][0] = n3b_knots_array[map_3b[itype][jtype][ktype]][0][0];
          min_cut_3b[itype][ktype][jtype][0] = n3b_knots_array[map_3b[itype][ktype][jtype]][0][0];

          knot_spacing_3b[itype][jtype][ktype][0] =
              n3b_knots_array[map_3b[itype][jtype][ktype]][0][4] -
              n3b_knots_array[map_3b[itype][jtype][ktype]][0][3];
          knot_spacing_3b[itype][ktype][jtype][0] = knot_spacing_3b[itype][jtype][ktype][0];

          temp_line = txtfilereader.next_line(num_knots_3b_ik);
          ValueTokenizer fp5th_line(temp_line);
          if ((int) fp5th_line.count() != num_knots_3b_ik)
            error->all(FLERR,
                       "UF3: Error reading the 3B potential block for {}-{}-{}\n"
                       "Expected {} numbers on 5th line of the block but found {} "
                       "numbers",
                       element1, element2, element3, num_knots_3b_ik, fp5th_line.count());

          for (int i = 0; i < num_knots_3b_ik; i++) {
            n3b_knots_array[map_3b[itype][jtype][ktype]][1][i] = fp5th_line.next_double();
            n3b_knots_array[map_3b[itype][ktype][jtype]][2][i] =
                n3b_knots_array[map_3b[itype][jtype][ktype]][1][i];
          }

          min_cut_3b[itype][jtype][ktype][1] = n3b_knots_array[map_3b[itype][jtype][ktype]][1][0];
          min_cut_3b[itype][ktype][jtype][2] = n3b_knots_array[map_3b[itype][ktype][jtype]][2][0];

          knot_spacing_3b[itype][jtype][ktype][1] =
              n3b_knots_array[map_3b[itype][jtype][ktype]][1][4] -
              n3b_knots_array[map_3b[itype][jtype][ktype]][1][3];
          knot_spacing_3b[itype][ktype][jtype][2] = knot_spacing_3b[itype][jtype][ktype][1];

          temp_line = txtfilereader.next_line(num_knots_3b_ij);
          ValueTokenizer fp6th_line(temp_line);
          if ((int) fp6th_line.count() != num_knots_3b_ij)
            error->all(FLERR,
                       "UF3: Error reading the 3B potential block for {}-{}-{}\n"
                       "Expected {} numbers on 6th line of the block but found {} "
                       "numbers",
                       element1, element2, element3, num_knots_3b_ij, fp6th_line.count());

          for (int i = 0; i < num_knots_3b_ij; i++) {
            n3b_knots_array[map_3b[itype][jtype][ktype]][2][i] = fp6th_line.next_double();
            n3b_knots_array[map_3b[itype][ktype][jtype]][1][i] =
                n3b_knots_array[map_3b[itype][jtype][ktype]][2][i];
          }

          min_cut_3b[itype][jtype][ktype][2] = n3b_knots_array[map_3b[itype][jtype][ktype]][2][0];
          min_cut_3b[itype][ktype][jtype][1] = n3b_knots_array[map_3b[itype][ktype][jtype]][1][0];

          knot_spacing_3b[itype][jtype][ktype][2] =
              n3b_knots_array[map_3b[itype][jtype][ktype]][2][4] -
              n3b_knots_array[map_3b[itype][jtype][ktype]][2][3];
          knot_spacing_3b[itype][ktype][jtype][1] = knot_spacing_3b[itype][jtype][ktype][2];

          //skip next line
          txtfilereader.skip_line();

          int coeff_matrix_dim1 = n3b_coeff_array_size[map_3b[itype][jtype][ktype]][0];
          int coeff_matrix_dim2 = n3b_coeff_array_size[map_3b[itype][jtype][ktype]][1];
          int coeff_matrix_dim3 = n3b_coeff_array_size[map_3b[itype][jtype][ktype]][2];

          if (num_knots_3b_jk != coeff_matrix_dim3 + 3 + 1)
            error->all(FLERR,
                       "UF3: {}-{}-{} interaction block has incorrect knot "
                       "(NUM_OF_KNOTS_JK) and coeff (coeff_matrix_dim3) data "
                       "nknots!=ncoeffs + 3 + 1",
                       element1, element2, element3);

          if (num_knots_3b_ik != coeff_matrix_dim2 + 3 + 1)
            error->all(FLERR,
                       "UF3: {}-{}-{} interaction block has incorrect knot "
                       "(NUM_OF_KNOTS_IK) and coeff (coeff_matrix_dim2) data "
                       "nknots!=ncoeffs + 3 + 1",
                       element1, element2, element3);

          if (num_knots_3b_ij != coeff_matrix_dim1 + 3 + 1)
            error->all(FLERR,
                       "UF3: {}-{}-{} interaction block has incorrect knot "
                       "(NUM_OF_KNOTS_IJ) and coeff (coeff_matrix_dim1) data "
                       "nknots!=ncoeffs + 3 + 1",
                       element1, element2, element3);

          int coeff_matrix_elements_len = coeff_matrix_dim3;
          int key1 = map_3b[itype][jtype][ktype];
          int key2 = map_3b[itype][ktype][jtype];

          int line_count = 0;
          for (int i = 0; i < coeff_matrix_dim1; i++) {
            for (int j = 0; j < coeff_matrix_dim2; j++) {
              temp_line = txtfilereader.next_line(coeff_matrix_elements_len);
              ValueTokenizer coeff_line(temp_line);
              if ((int) coeff_line.count() != coeff_matrix_elements_len)
                error->all(FLERR,
                           "UF3: Error reading 3B potential block for {}-{}-{}\n"
                           "Expected {} numbers on {}th line of the block but found {} "
                           "numbers",
                           element1, element2, element3, coeff_matrix_elements_len, line_count + 8,
                           coeff_line.count());

              for (int k = 0; k < coeff_matrix_dim3; k++) {
                n3b_coeff_array[key1][i][j][k] = coeff_line.next_double();
              }
              line_count += 1;
            }
          }

          for (int i = 0; i < coeff_matrix_dim1; i++) {
            for (int j = 0; j < coeff_matrix_dim2; j++) {
              for (int k = 0; k < coeff_matrix_dim3; k++) {
                n3b_coeff_array[key2][j][i][k] = n3b_coeff_array[key1][i][j][k];
              }
            }
          }

          setflag_3b[itype][jtype][ktype] = setflag_3b[itype][ktype][jtype] = 1;
        }
      }
    }    // if #UF3 POT
  }    //while
  fclose(fp);

  //Set interaction of atom types of the same elements
  for (int i = 1; i < np1; i++) {
    for (int j = 1; j < np1; j++) {
      if (setflag[i][j] != 1) {
        //i-j interaction not set
        //maybe i-j is mapped to some other atom type interaction?
        int i_mapped_to = element_to_type[elements[i]];
        int j_mapped_to = element_to_type[elements[j]];
        if ((i_mapped_to == i) && (j_mapped_to == j))
          //i-j is not mapped to some other atom type ie interaction is missing on file
          error->all(FLERR, "UF3: Potential for interaction {}-{} ie {}-{} not found in {}",
                     i, j, elements[i_mapped_to], elements[j_mapped_to], potf_name);
        cut_2b[i][j] = cut_2b[i_mapped_to][j_mapped_to];
        n2b_knots_array_size[i][j] = n2b_knots_array_size[i_mapped_to][j_mapped_to];
        n2b_coeff_array_size[i][j] = n2b_coeff_array_size[i_mapped_to][j_mapped_to];
        knot_spacing_type_2b[i][j] = knot_spacing_type_2b[i_mapped_to][j_mapped_to];
        knot_spacing_2b[i][j] = knot_spacing_2b[i_mapped_to][j_mapped_to];
        for (int knot_no = 0; knot_no < max_num_knots_2b; knot_no++)
          n2b_knots_array[i][j][knot_no] = n2b_knots_array[i_mapped_to][j_mapped_to][knot_no];
        for (int coeff_no = 0; coeff_no < max_num_coeff_2b; coeff_no++)
          n2b_coeff_array[i][j][coeff_no] = n2b_coeff_array[i_mapped_to][j_mapped_to][coeff_no];
        setflag[i][j] = 1;
      }
    }
  }

  if (pot_3b) {
    for (int i = 1; i < np1; i++) {
      for (int j = 1; j < np1; j++) {
        for (int k = 1; k < np1; k++) {
          if (setflag_3b[i][j][k] != 1) {
            //i-j-k interaction not set
            //maybe i-j-k is mapped to some other atom type interaction?
            const int i_mapped_to = element_to_type[elements[i]];
            const int j_mapped_to = element_to_type[elements[j]];
            const int k_mapped_to = element_to_type[elements[k]];
            if ((i_mapped_to == i) && (j_mapped_to == j) && (k_mapped_to == k))
              error->all(FLERR, "UF3: Potential for interaction {}-{}-{} ie {}-{}-{} not found in {}",
                         i, j, k, elements[i_mapped_to], elements[j_mapped_to],
                         elements[k_mapped_to], potf_name);
            if (setflag_3b[i_mapped_to][j_mapped_to][k_mapped_to] != 1)
              error->all(FLERR, "UF3: Interaction {}-{}-{} mapped to {}-{}-{} but not found in {}",
                         i, j, k, i_mapped_to, j_mapped_to, k_mapped_to, potf_name);
            cut_3b_list[i][j] = std::max(cut_3b_list[i_mapped_to][j_mapped_to], cut_3b_list[i][j]);
            cut_3b[i][j][k] = cut_3b[i_mapped_to][j_mapped_to][k_mapped_to];
            knot_spacing_type_3b[i][j][k] =
                knot_spacing_type_3b[i_mapped_to][j_mapped_to][k_mapped_to];
            knot_spacing_3b[i][j][k][0] = knot_spacing_3b[i_mapped_to][j_mapped_to][k_mapped_to][0];
            knot_spacing_3b[i][j][k][1] = knot_spacing_3b[i_mapped_to][j_mapped_to][k_mapped_to][1];
            knot_spacing_3b[i][j][k][2] = knot_spacing_3b[i_mapped_to][j_mapped_to][k_mapped_to][2];
            const int key = map_3b[i][j][k];
            const int mapped_to_key = map_3b[i_mapped_to][j_mapped_to][k_mapped_to];
            auto n3b_knots_size = n3b_knots_array_size[key];
            auto n3b_knots = n3b_knots_array[key];
            auto n3b_coeff_size = n3b_coeff_array_size[key];
            auto n3b_coeff = n3b_coeff_array[key];
            auto n3b_knots_size_mapped = n3b_knots_array_size[mapped_to_key];
            auto n3b_knots_mapped = n3b_knots_array[mapped_to_key];
            auto n3b_coeff_size_mapped = n3b_coeff_array_size[mapped_to_key];
            auto n3b_coeff_mapped = n3b_coeff_array[mapped_to_key];
            n3b_knots_size[0] = n3b_knots_size_mapped[0];
            n3b_knots_size[1] = n3b_knots_size_mapped[1];
            n3b_knots_size[2] = n3b_knots_size_mapped[2];
            n3b_coeff_size[0] = n3b_coeff_size_mapped[0];
            n3b_coeff_size[1] = n3b_coeff_size_mapped[1];
            n3b_coeff_size[2] = n3b_coeff_size_mapped[2];
            min_cut_3b[i][j][k][0] = min_cut_3b[i_mapped_to][j_mapped_to][k_mapped_to][0];
            min_cut_3b[i][j][k][1] = min_cut_3b[i_mapped_to][j_mapped_to][k_mapped_to][1];
            min_cut_3b[i][j][k][2] = min_cut_3b[i_mapped_to][j_mapped_to][k_mapped_to][2];
            for (int knot_no = 0; knot_no < n3b_knots_size[0]; knot_no++) n3b_knots[0][knot_no] = n3b_knots_mapped[0][knot_no];
            for (int knot_no = 0; knot_no < n3b_knots_size[1]; knot_no++) n3b_knots[1][knot_no] = n3b_knots_mapped[1][knot_no];
            for (int knot_no = 0; knot_no < n3b_knots_size[2]; knot_no++) n3b_knots[2][knot_no] = n3b_knots_mapped[2][knot_no];
            for (int coeff1 = 0; coeff1 < n3b_coeff_size[0]; coeff1++)
              for (int coeff2 = 0; coeff2 < n3b_coeff_size[1]; coeff2++)
                for (int coeff3 = 0; coeff3 < n3b_coeff_size[2]; coeff3++)
                  n3b_coeff[coeff1][coeff2][coeff3] = n3b_coeff_mapped[coeff1][coeff2][coeff3];
            setflag_3b[i][j][k] = 1;
          }
        }
      }
    }
  }
}

//Broadcast data read from potential file to all processors
void UF3Potential::communicate()
{
  const int np1 = atom->ntypes + 1;
  MPI_Bcast(&cut_2b[0][0], np1 * np1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&n2b_knots_array_size[0][0], np1 * np1, MPI_INT, 0, world);
  MPI_Bcast(&n2b_coeff_array_size[0][0], np1 * np1, MPI_INT, 0, world);
  MPI_Bcast(&max_num_knots_2b, 1, MPI_INT, 0, world);
  MPI_Bcast(&max_num_coeff_2b, 1, MPI_INT, 0, world);

  if (pot_3b) {
    MPI_Bcast(&cut_3b_list[0][0], np1 * np1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&cut_3b[0][0][0], np1 * np1 * np1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&n3b_knots_array_size[0][0], tot_interaction_count_3b * 3, MPI_INT, 0, world);
    MPI_Bcast(&n3b_coeff_array_size[0][0], tot_interaction_count_3b * 3, MPI_INT, 0, world);
    MPI_Bcast(&max_num_knots_3b, 1, MPI_INT, 0, world);
    MPI_Bcast(&max_num_coeff_3b, 1, MPI_INT, 0, world);
  }

  if (comm->me != 0) {
    memory->destroy(n2b_knots_array);
    memory->destroy(n2b_coeff_array);
    memory->create(n2b_knots_array, np1, np1, max_num_knots_2b, "pair:n2b_knots_array");
    memory->create(n2b_coeff_array, np1, np1, max_num_coeff_2b, "pair:n2b_coeff_array");

    if (pot_3b) {
      memory->destroy(n3b_knots_array);
      memory->destroy(n3b_coeff_array);
      memory->create(n3b_knots_array, tot_interaction_count_3b, 3,
                     max_num_knots_3b, "pair:n3b_knots_array");
      memory->create(n3b_coeff_array, tot_interaction_count_3b, max_num_coeff_3b, max_num_coeff_3b,
                     max_num_coeff_3b, "pair:n3b_coeff_array");
    }
  }

  MPI_Bcast(&knot_spacing_type_2b[0][0], np1 * np1, MPI_INT, 0, world);
  MPI_Bcast(&knot_spacing_2b[0][0], np1 * np1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&n2b_knots_array[0][0][0], np1 * np1 * max_num_knots_2b, MPI_DOUBLE, 0, world);
  MPI_Bcast(&n2b_coeff_array[0][0][0], np1 * np1 * max_num_coeff_2b, MPI_DOUBLE, 0, world);
  MPI_Bcast(&setflag[0][0], np1 * np1, MPI_INT, 0, world);

  if (pot_3b) {
    MPI_Bcast(&knot_spacing_type_3b[0][0][0], np1 * np1 * np1, MPI_INT, 0, world);
    MPI_Bcast(&knot_spacing_3b[0][0][0][0], np1 * np1 * np1 * 3, MPI_DOUBLE, 0, world);
    MPI_Bcast(&n3b_knots_array[0][0][0], tot_interaction_count_3b * 3 * max_num_knots_3b, MPI_DOUBLE, 0, world);
    MPI_Bcast(&setflag_3b[0][0][0], np1 * np1 * np1, MPI_INT, 0, world);
    MPI_Bcast(&min_cut_3b[0][0][0][0], np1 * np1 * np1 * 3, MPI_DOUBLE, 0, world);
    MPI_Bcast(&n3b_coeff_array[0][0][0][0],
              tot_interaction_count_3b * max_num_coeff_3b * max_num_coeff_3b * max_num_coeff_3b,
              MPI_DOUBLE, 0, world);
  }
}




void UF3Potential::create_bsplines()
{
  const int ntypes = atom->ntypes;
  bsplines_created = 1;
  const int spacing_type = knot_spacing_type_2b[1][1];
  for (int i = 1; i < ntypes + 1; i++) {
    for (int j = 1; j < ntypes + 1; j++) {
      if (setflag[i][j] != 1) error->all(FLERR, "UF3: 2-body {}-{} interaction not set",i, j);
    }
  }
  if (pot_3b) {
    for (int i = 1; i < ntypes + 1; i++) {
      for (int j = 1; j < ntypes + 1; j++) {
        for (int k = 1; k < ntypes + 1; k++) {
          if (setflag_3b[i][j][k] != 1) error->all(FLERR, "UF3: 3-body UF {}-{}-{} interaction not set", i, j, k);
          if (spacing_type != knot_spacing_type_3b[i][j][k])
            error->all(FLERR,
                       "UF3: In the current version the knot spacing type, "
                       "for all interactions needs to be same. For {}-{}-{} "
                       "i.e. {}-{}-{} interaction expected {}, but found {}",
                       i, j, k, elements[i], elements[j], elements[k], spacing_type,
                       knot_spacing_type_3b[i][j][k]);
        }
      }
    }
  }

  if (spacing_type) {
    get_starting_index_2b_ptr = &UF3Potential::get_starting_index_nonuniform_2b;
    if (pot_3b) get_starting_index_3b_ptr = &UF3Potential::get_starting_index_nonuniform_3b;
  } else {
    get_starting_index_2b_ptr = &UF3Potential::get_starting_index_uniform_2b;
    if (pot_3b) get_starting_index_3b_ptr = &UF3Potential::get_starting_index_uniform_3b;
  }

  create_cached_constants_2b();
  if (pot_3b) create_cached_constants_3b();
}

int UF3Potential::get_starting_index_2b(int i, int j, double r)
{
  return (this->*get_starting_index_2b_ptr)(i, j, r);
}

int UF3Potential::get_starting_index_3b(int i, int j, int k, double r, int knot_dim)
{
  return (this->*get_starting_index_3b_ptr)(i, j, k, r, knot_dim);
}

int UF3Potential::get_starting_index_uniform_2b(int i, int j, double r)
{
  return 3 + (int) ((r - n2b_knots_array[i][j][0]) / (knot_spacing_2b[i][j]));
}

int UF3Potential::get_starting_index_uniform_3b(int i, int j, int k, double r, int knot_dim)
{
  return 3 +
      (int) (((r - n3b_knots_array[map_3b[i][j][k]][knot_dim][0]) /
              knot_spacing_3b[i][j][k][knot_dim]));
}

int UF3Potential::get_starting_index_nonuniform_2b(int i, int j, double r)
{
  for (int l = 3; l < n2b_knots_array_size[i][j] - 1; ++l) {
    if ((n2b_knots_array[i][j][l] <= r) && (r < n2b_knots_array[i][j][l + 1])) return l;
  }
  return -1;
}

int UF3Potential::get_starting_index_nonuniform_3b(int i, int j, int k, double r, int knot_dim)
{
  for (int l = 3; l < n3b_knots_array_size[map_3b[i][j][k]][knot_dim] - 1; ++l) {
    if ((n3b_knots_array[map_3b[i][j][k]][knot_dim][l] <= r) &&
        (r < n3b_knots_array[map_3b[i][j][k]][knot_dim][l + 1]))
      return l;
  }
  return -1;
}

void UF3Potential::create_cached_constants_2b()
{
  const int np1 = atom->ntypes + 1;
  memory->destroy(cached_constants_2b);
  memory->destroy(cached_constants_2b_deri);
  memory->create(cached_constants_2b, np1, np1, max_num_coeff_2b, 16, "uf3:cached_constants_2b");
  memory->create(cached_constants_2b_deri, np1, np1, max_num_coeff_2b - 1, 9, "uf3:cached_constants_2b_deri");

  for (int i = 1; i < np1; i++) {
    for (int j = 1; j < np1; j++) {
      for (int l = 0; l < n2b_coeff_array_size[i][j]; l++) {
        uf3_bspline_basis3 bspline_basis(lmp, &n2b_knots_array[i][j][l], n2b_coeff_array[i][j][l]);
        for (int cc = 0; cc < 16; cc++) cached_constants_2b[i][j][l][cc] = bspline_basis.constants[cc];
      }
    }
  }

  for (int i = 1; i < np1; i++) {
    for (int j = 1; j < np1; j++) {
      //initialize coeff and knots for derivative
      double *knots_for_deri = nullptr;
      auto n2b_knots_size_ij = n2b_knots_array_size[i][j];
      memory->create(knots_for_deri, n2b_knots_size_ij - 2, "pair:knots_for_deri");
      for (int l = 1; l < n2b_knots_size_ij - 1; l++) knots_for_deri[l - 1] = n2b_knots_array[i][j][l];
      double *coeff_for_deri = nullptr;
      memory->create(coeff_for_deri, n2b_knots_size_ij - 1, "pair:coeff_for_deri");
      for (int l = 0; l < n2b_coeff_array_size[i][j] - 1; l++) {
        double dntemp = 3 / (n2b_knots_array[i][j][l + 4] - n2b_knots_array[i][j][l + 1]);
        coeff_for_deri[l] = (n2b_coeff_array[i][j][l + 1] - n2b_coeff_array[i][j][l]) * dntemp;
      }
      for (int l = 0; l < n2b_coeff_array_size[i][j] - 1; l++) {
        uf3_bspline_basis2 bspline_basis_deri(lmp, &knots_for_deri[l], coeff_for_deri[l]);
        for (int cc = 0; cc < 9; cc++) {
          cached_constants_2b_deri[i][j][l][cc] = bspline_basis_deri.constants[cc];
        }
      }
      memory->destroy(knots_for_deri);
      memory->destroy(coeff_for_deri);
    }
  }
}

void UF3Potential::create_cached_constants_3b()
{
  const int ntypes = atom->ntypes;
  memory->destroy(coeff_for_der_jk);
  memory->destroy(coeff_for_der_ik);
  memory->destroy(coeff_for_der_ij);
  memory->destroy(cached_constants_3b);
  memory->destroy(cached_constants_3b_deri);

  memory->create(coeff_for_der_jk, tot_interaction_count_3b, max_num_coeff_3b, max_num_coeff_3b,
                 max_num_coeff_3b, "pair:coeff_for_der_jk");

  memory->create(coeff_for_der_ik, tot_interaction_count_3b, max_num_coeff_3b, max_num_coeff_3b,
                 max_num_coeff_3b, "pair:coeff_for_der_ik");

  memory->create(coeff_for_der_ij, tot_interaction_count_3b, max_num_coeff_3b, max_num_coeff_3b,
                 max_num_coeff_3b, "pair:coeff_for_der_ij");

  memory->create(cached_constants_3b, tot_interaction_count_3b, 3, max_num_coeff_3b, 16,
                 "pair:cached_constants_3b");

  memory->create(cached_constants_3b_deri, tot_interaction_count_3b, 3, max_num_coeff_3b - 1, 9,
                 "pair:cached_constants_3b_deri");

  for (int i = 1; i < ntypes + 1; i++) {
    for (int j = 1; j < ntypes + 1; j++) {
      for (int k = 1; k < ntypes + 1; k++) {
        const int map_to = map_3b[i][j][k];
        for (int l = 0; l < n3b_knots_array_size[map_to][2] - 4; l++) {
          uf3_bspline_basis3 bspline_basis_ij(lmp, &n3b_knots_array[map_to][2][l], 1);
          for (int cc = 0; cc < 16; cc++)
            cached_constants_3b[map_to][0][l][cc] = bspline_basis_ij.constants[cc];
        }
        for (int l = 0; l < n3b_knots_array_size[map_to][1] - 4; l++) {
          uf3_bspline_basis3 bspline_basis_ik(lmp, &n3b_knots_array[map_to][1][l], 1);
          for (int cc = 0; cc < 16; cc++)
            cached_constants_3b[map_to][1][l][cc] = bspline_basis_ik.constants[cc];
        }
        for (int l = 0; l < n3b_knots_array_size[map_to][0] - 4; l++) {
          uf3_bspline_basis3 bspline_basis_jk(lmp, &n3b_knots_array[map_to][0][l], 1);
          for (int cc = 0; cc < 16; cc++)
            cached_constants_3b[map_to][2][l][cc] = bspline_basis_jk.constants[cc];
        }
      }
    }
  }

  for (int i = 1; i < ntypes + 1; i++) {
    for (int j = 1; j < ntypes + 1; j++) {
      for (int k = 1; k < ntypes + 1; k++) {
        const int map_to = map_3b[i][j][k];
        double **knots_for_der = nullptr;

        //n3b_knots_array_size[map_to][0] for jk knot vector --> always largest
        memory->create(knots_for_der, 3, n3b_knots_array_size[map_to][0] - 1, "pair:knots_for_der");

        //--deri_basis_jk
        for (int l = 1; l < n3b_knots_array_size[map_to][0] - 1; l++)
          knots_for_der[0][l - 1] = n3b_knots_array[map_to][0][l];

        for (int l = 0; l < n3b_coeff_array_size[map_to][0]; l++) {
          for (int m = 0; m < n3b_coeff_array_size[map_to][1]; m++) {
            for (int n = 0; n < n3b_coeff_array_size[map_to][2] - 1; n++) {
              double dntemp =
                  3 / (n3b_knots_array[map_to][0][n + 4] - n3b_knots_array[map_to][0][n + 1]);
              coeff_for_der_jk[map_to][l][m][n] =
                  ((n3b_coeff_array[map_to][l][m][n + 1] - n3b_coeff_array[map_to][l][m][n]) *
                   dntemp);
            }
          }
        }

        //--deri_basis_ik
        for (int l = 1; l < n3b_knots_array_size[map_to][1] - 1; l++)
          knots_for_der[1][l - 1] = n3b_knots_array[map_to][1][l];

        for (int l = 0; l < n3b_coeff_array_size[map_to][0]; l++) {
          for (int m = 0; m < n3b_coeff_array_size[map_to][1] - 1; m++) {
            double dntemp =
                3 / (n3b_knots_array[map_to][1][m + 4] - n3b_knots_array[map_to][1][m + 1]);
            for (int n = 0; n < n3b_coeff_array_size[map_to][2]; n++) {
              coeff_for_der_ik[map_to][l][m][n] =
                  ((n3b_coeff_array[map_to][l][m + 1][n] - n3b_coeff_array[map_to][l][m][n]) *
                   dntemp);
            }
          }
        }

        //--deri_basis_ij
        for (int l = 1; l < n3b_knots_array_size[map_to][2] - 1; l++)
          knots_for_der[2][l - 1] = n3b_knots_array[map_to][2][l];

        for (int l = 0; l < n3b_coeff_array_size[map_to][0] - 1; l++) {
          double dntemp =
              3 / (n3b_knots_array[map_to][2][l + 4] - n3b_knots_array[map_to][2][l + 1]);
          for (int m = 0; m < n3b_coeff_array_size[map_to][1]; m++) {
            for (int n = 0; n < n3b_coeff_array_size[map_to][2]; n++) {
              coeff_for_der_ij[map_to][l][m][n] =
                  ((n3b_coeff_array[map_to][l + 1][m][n] - n3b_coeff_array[map_to][l][m][n]) *
                   dntemp);
            }
          }
        }

        for (int l = 0; l < n3b_coeff_array_size[map_to][0] - 1; l++) {
          uf3_bspline_basis2 bspline_basis_deri_ij(lmp, &knots_for_der[2][l], 1);
          for (int cc = 0; cc < 9; cc++) {
            cached_constants_3b_deri[map_to][0][l][cc] = bspline_basis_deri_ij.constants[cc];
          }
        }

        for (int l = 0; l < n3b_coeff_array_size[map_to][1] - 1; l++) {
          uf3_bspline_basis2 bspline_basis_deri_ik(lmp, &knots_for_der[1][l], 1);
          for (int cc = 0; cc < 9; cc++) {
            cached_constants_3b_deri[map_to][1][l][cc] = bspline_basis_deri_ik.constants[cc];
          }
        }

        for (int l = 0; l < n3b_coeff_array_size[map_to][2] - 1; l++) {
          uf3_bspline_basis2 bspline_basis_deri_jk(lmp, &knots_for_der[0][l], 1);
          for (int cc = 0; cc < 9; cc++) {
            cached_constants_3b_deri[map_to][2][l][cc] = bspline_basis_deri_jk.constants[cc];
          }
        }

        memory->destroy(knots_for_der);
      }
    }
  }
}


/* ---------------------------------------------------------------------- */

double UF3Potential::memory_usage()
{
  const int ntypes = atom->ntypes;

  double bytes = (double) (ntypes + 1) * (ntypes + 1) * (ntypes + 1) * sizeof(int); //***setflag_3b
  bytes += (double) (ntypes + 1) * (ntypes + 1) * sizeof(int);    //knot_spacing_type_2b
  bytes += (double) (ntypes + 1) * (ntypes + 1) * (ntypes + 1) * sizeof(int);    //knot_spacing_type_3b
  bytes += (double) (ntypes + 1) * (ntypes + 1) * sizeof(double);    //cut
  bytes += (double) (ntypes + 1) * (ntypes + 1) * (ntypes + 1) * sizeof(double);    //***cut_3b
  bytes += (double) (ntypes + 1) * (ntypes + 1) * sizeof(double);    //cut_3b_list
  bytes += (double) (ntypes + 1) * (ntypes + 1) * (ntypes + 1) * 3 * sizeof(double);    //min_cut_3b
  bytes += (double) (ntypes + 1) * (ntypes + 1) * sizeof(double);    //knot_spacing_2b
  bytes += (double) (ntypes + 1) * (ntypes + 1) * (ntypes + 1) * sizeof(double);    //knot_spacing_3b
  bytes += (double) (ntypes + 1) * (ntypes + 1) * max_num_knots_2b * sizeof(double);    //n2b_knots_array
  bytes += (double) (ntypes + 1) * (ntypes + 1) * max_num_coeff_2b * sizeof(double);    //n2b_coeff_array
  bytes += (double) (ntypes + 1) * (ntypes + 1) * sizeof(int);    //n2b_knots_array_size
  bytes += (double) (ntypes + 1) * (ntypes + 1) * sizeof(int);    //n2b_coeff_array_size
  bytes += (double) (ntypes + 1) * (ntypes + 1) * max_num_coeff_2b * 16 * sizeof(double);    //cached_constants_2b,
  bytes += (double) (ntypes + 1) * (ntypes + 1) * (max_num_coeff_2b - 1) * 9 * sizeof(double);    //cached_constants_2b_deri

  if (pot_3b) {
    bytes += (double) (ntypes + 1) * (ntypes + 1) * (ntypes + 1) * sizeof(int);    //map_3b
    bytes += (double) tot_interaction_count_3b * 3 * max_num_knots_3b * sizeof(double);    //n3b_knots_array
    bytes += (double) tot_interaction_count_3b * max_num_coeff_3b * max_num_coeff_3b * max_num_coeff_3b * sizeof(double);    //n3b_coeff_array

    bytes += (double) tot_interaction_count_3b * 3 * sizeof(int);    //n3b_knots_array_size
    bytes += (double) tot_interaction_count_3b * 3 * sizeof(int);    //n3b_coeff_array_size

    bytes += (double) tot_interaction_count_3b * max_num_coeff_3b * max_num_coeff_3b *
        max_num_coeff_3b * 3 * sizeof(double);    //coeff_for_der_jk coeff_for_der_ik coeff_for_der_ij

    bytes += (double) tot_interaction_count_3b * 3 * max_num_coeff_3b * 16 * sizeof(double);    //cached_constants_3b
    bytes += (double) tot_interaction_count_3b * 3 * (max_num_coeff_3b - 1) * 9 * sizeof(double);    //cached_constants_3b_deri
  }

  bytes += (double) 4 * sizeof(int);     //max_num_knots_2b, max_num_coeff_2b, max_num_knots_3b, max_num_coeff_3b
  bytes += (double) 1 * sizeof(bool);    //pot_3b

  return bytes;
}



