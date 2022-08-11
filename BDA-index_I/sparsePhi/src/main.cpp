/**
 * @file    main.cpp
 * @section LICENCE
 *
 * This file is part of EM-SparsePhi v0.1.0
 * See: http://www.cs.helsinki.fi/group/pads/
 *
 * Copyright (C) 2016
 *   Juha Karkkainen <juha.karkkainen (at) cs.helsinki.fi>
 *   Dominik Kempa <dominik.kempa (at) gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 **/

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <ctime>
#include <string>
#include <sstream>
#include <getopt.h>
#include <unistd.h>

#include "uint40.hpp"
#include "uint48.hpp"
#include "em_sparse_phi_src/em_sparse_phi.hpp"


char *program_name;

void usage(int status) {
  printf(

"Usage: %s [OPTION]... FILE\n"
"Construct the LCP array for text stored in FILE.\n"
"\n"
"Mandatory arguments to long options are mandatory for short options too.\n"
"  -h, --help              display this help and exit\n"
"  -i, --intsize=SIZE      use integers of SIZE bytes (default: 5). Currently\n"
"                          supported values are 4, 5, 6, and 8\n"
"  -m, --mem=MEM           use MEM MiB of RAM for computation (default: 3584)\n"
"  -o, --output=OUTFILE    specify output filename (default: FILE.lcpX, where\n"
"                          X = integer size, see the -i flag)\n"
"  -s, --sa=SUFARRAY       specify the location of the suffix array of FILE\n"
"                          (default: FILE.saX, X = integer size, see -i flag)\n",
    program_name);

  std::exit(status);
}

bool file_exists(std::string fname) {
  std::FILE *f = std::fopen(fname.c_str(), "r");
  bool ret = (f != NULL);
  if (f != NULL) std::fclose(f);

  return ret;
}

template<typename int_type>
std::string intToStr(int_type x) {
  std::stringstream ss;
  ss << x;
  return ss.str();
}

int main(int argc, char **argv) {
  srand(time(0) + getpid());
  program_name = argv[0];

  static struct option long_options[] = {
    {"help",    no_argument,       NULL, 'h'},
    {"intsize", required_argument, NULL, 'i'},
    {"mem",     required_argument, NULL, 'm'},
    {"output",  required_argument, NULL, 'o'},
    {"sa",      required_argument, NULL, 's'},
    {NULL, 0, NULL, 0}
  };

  std::uint64_t ram_use = 3584UL << 20;
  std::uint64_t int_size = 5;
  std::string sa_filename("");
  std::string out_filename("");

  // Parse command-line options.
  int c;
  while ((c = getopt_long(argc, argv, "hi:m:o:s:", long_options, NULL)) != -1) {
    switch(c) {
      case 'h':
        usage(EXIT_FAILURE);
      case 'i':
        int_size = std::atol(optarg);
        if (!(int_size == 4 || int_size == 5 || int_size == 6 || int_size == 8)) {
          fprintf(stderr, "Error: invalid int size (%lu)\n\n", int_size);
          usage(EXIT_FAILURE);
        }
        break;
      case 'm':
        ram_use = std::atol(optarg) << 20;
        if (ram_use == 0) {
          fprintf(stderr, "Error: invalid RAM limit (%lu)\n\n", ram_use);
          usage(EXIT_FAILURE);
        }
        break;
      case 'o':
        out_filename = std::string(optarg);
        break;
      case 's':
        sa_filename = std::string(optarg);
        break;
      default:
        usage(EXIT_FAILURE);
    }
  }

  if (optind >= argc) {
    fprintf(stderr, "Error: FILE not provided\n\n");
    usage(EXIT_FAILURE);
  }

  // Parse the text filename.
  std::string text_filename = std::string(argv[optind++]);
  if (optind < argc) {
    fprintf(stderr, "Warning: multiple input files provided. "
    "Only the first will be processed.\n");
  }

  // Set default SA and output filenames (if not provided).
  if (sa_filename.empty()) sa_filename = text_filename + ".sa" + intToStr(int_size);
  if (out_filename.empty()) out_filename = text_filename + ".lcp" + intToStr(int_size);

  // Check for the existence of text and suffix array.
  if (!file_exists(text_filename)) {
    fprintf(stderr, "Error: input file (%s) does not exist\n\n",
        text_filename.c_str());
    usage(EXIT_FAILURE);
  }
  if (!file_exists(sa_filename)) {
    fprintf(stderr, "Error: suffix array (%s) does not exist\n\n",
        sa_filename.c_str());
    usage(EXIT_FAILURE);
  }

  if (file_exists(out_filename)) {
    // Output file exists, should we proceed?
    char *line = NULL;
    std::uint64_t buflen = 0;
    std::int64_t len = 0L;

    do {
      printf("Output file (%s) exists. Overwrite? [y/n]: ",
          out_filename.c_str());
      if ((len = getline(&line, &buflen, stdin)) == -1) {
        printf("\nError: failed to read answer\n\n");
        std::fflush(stdout);
        usage(EXIT_FAILURE);
      }
    } while (len != 2 || (line[0] != 'y' && line[0] != 'n'));

    if (line[0] == 'n') {
      free(line);
      std::exit(EXIT_FAILURE);
    }
    free(line);
  }

  // Run the EM LCP array construction algorithm.
  if (int_size == 4) em_sparse_phi<std::uint32_t>(text_filename, sa_filename, out_filename, ram_use);
  else if (int_size == 5) em_sparse_phi<uint40>(text_filename, sa_filename, out_filename, ram_use);
  else if (int_size == 6) em_sparse_phi<uint48>(text_filename, sa_filename, out_filename, ram_use);
  else em_sparse_phi<std::uint64_t>(text_filename, sa_filename, out_filename, ram_use);
}
