/*
 * Copyright (c) 2020 Maxime Schmitt <maxime.schmitt@manchester.ac.uk>
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors
 * may be used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include "fdtd.h"
#include "fdtd_common.h"
#include "initialize.h"
#include "time_measurement.h"
#include <getopt.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>

static struct option opt_options[] = {
    {"one-dimensional", no_argument, 0, '1'},
    {"two-dimensional", no_argument, 0, '2'},
    {"three-dimensional", no_argument, 0, '3'},
    {"setup-id", required_argument, 0, 's'},
    {"size-x", required_argument, 0, 'x'},
    {"size-y", required_argument, 0, 'y'},
    {"size-z", required_argument, 0, 'z'},
    {"output", optional_argument, 0, 'o'},
    {"courant-friedrichs-levy-condition", required_argument, 0, 'c'},
    {"smallest-wavelength", required_argument, 0, 'w'},
    {"cpml-absorbing-thickness", required_argument, 0, 'a'},
    {"stop-sim-time", required_argument, 0, 't'},
    {"num-iterations", required_argument, 0, 'i'},
    {"help", no_argument, 0, 'h'},
    {"quiet", no_argument, 0, 'q'},
    {0, 0, 0, 0}};

static const char options[] = ":123s:x:y:z:o:c:w:a:t:i:hq";

static const char help_string[] =
    "Options:"
    "\n  -1 --one-dimensional     : 1D solver"
    "\n  -2 --one-dimensional     : 2D solver"
    "\n  -3 --one-dimensional     : 3D solver"
    "\n  -s --setup-id            : Predefined problem identifier"
    "\n                             1D : 0 - Half air half water, "
    "left-to-right gaussian"
    "\n                             2D : 0 - West air east water, "
    "west-to-east gaussian"
    "\n                                  1 - Air with high "
    "permittivity centered object, west pulse"
    "\n                                  2 - Centered gaussian "
    "excitation in free space"
    "\n                             3D : 0 - Half air half water, "
    "west-to-east gaussian"
    "\n                                  1 - Air with high permittivity "
    "centered object"
    "\n  -x --size-x              : Size of the domain (e.g. 0.00001)"
    "\n  -y --size-y              : Size of the domain (e.g. 0.00001)"
    "\n  -z --size-z              : Size of the domain (e.g. 0.00001)"
    "\n  -o --output              : Select the output file name"
    "\n  -c --courant-friedrichs-levy-condition : Stability value"
    "\n  -w --smallest-wavelength : Smallest wavelength in the simulation"
    "\n  -a --cpml-absorbing-thickness  : Size of the absorbing boundary wall"
    "\n  -t --stop-sim-time       : Stop the simulation when the time is "
    "reached"
    "\n  -i --num-iterations      : Stop the simulation after the specified "
    "amount of solver iterations"
    "\n  -h --help                : Print this help"
    "\n  -q --quiet               : Do not print information to the user from "
    "inside the main kernel";

#define default_domain_size float_cst(0.00001)
#define default_cpml_width 20
#define default_smallest_wavelength float_cst(450e-9)
#define default_Sc float_cst(-1.)
#define default_end_time float_cst(-1.)
#define default_iteration_count 400

int main(int argc, char **argv) {
  float_type domain_size[3] = {default_domain_size, default_domain_size,
                               default_domain_size};
  char *output_filename = NULL;
  float_type Sc = default_Sc;
  float_type smallest_wavelength = default_smallest_wavelength;
  unsigned border_cpml_width = default_cpml_width;
  unsigned dimension = 1;
  unsigned setup_id = 0; // Default to the first one for each the dimension
  float_type end_time = default_end_time;
  size_t num_iterations = default_iteration_count;
  bool verbose = true;

  while (true) {
    int sscanf_return;
    int optchar = getopt_long(argc, argv, options, opt_options, NULL);
    if (optchar == -1)
      break;
    switch (optchar) {
    case '1':
      dimension = 1;
      break;
    case '2':
      dimension = 2;
      break;
    case '3':
      dimension = 3;
      break;
    case 'o':
      if (optarg != NULL && optarg[0] != '-' && optarg[0] != '\0') {
        output_filename = optarg;
      } else {
        output_filename = "gridData.dat";
        if (optarg != NULL)
          optind--;
      }
      break;
    case 'q':
      verbose = false;
      break;
    case 'x':
#if float_type == double
      sscanf_return = sscanf(optarg, "%lf", &domain_size[0]);
#else
      sscanf_return = sscanf(optarg, "%f", &domain_size[0]);
#endif
      if (sscanf_return == EOF || sscanf_return == 0 ||
          domain_size[0] < float_cst(0.)) {
        fprintf(stderr,
                "Please enter a positive floating point number for the domain "
                "size instead of \"-%c %s\"\n",
                optchar, optarg);
        domain_size[0] = default_domain_size;
      }
      break;
    case 'y':
#if float_type == double
      sscanf_return = sscanf(optarg, "%lf", &domain_size[1]);
#else
      sscanf_return = sscanf(optarg, "%f", &domain_size[1]);
#endif
      if (sscanf_return == EOF || sscanf_return == 0 ||
          domain_size[1] < float_cst(0.)) {
        fprintf(stderr,
                "Please enter a positive floating point number for the domain "
                "size instead of \"-%c %s\"\n",
                optchar, optarg);
        domain_size[1] = default_domain_size;
      }
      break;
    case 'z':
#if float_type == double
      sscanf_return = sscanf(optarg, "%lf", &domain_size[2]);
#else
      sscanf_return = sscanf(optarg, "%f", &domain_size[2]);
#endif
      if (sscanf_return == EOF || sscanf_return == 0 ||
          domain_size[2] < float_cst(0.)) {
        fprintf(stderr,
                "Please enter a positive floating point number for the domain "
                "size instead of \"-%c %s\"\n",
                optchar, optarg);
        domain_size[2] = default_domain_size;
      }
      break;
    case 'a':
      sscanf_return = sscanf(optarg, "%u", &border_cpml_width);
      if (sscanf_return == EOF || sscanf_return == 0) {
        fprintf(stderr,
                "Please enter a positive integer for the cpml thickness "
                "instead of \"-%c %s\"\n",
                optchar, optarg);
        border_cpml_width = default_cpml_width;
      }
      break;
    case 's':
      sscanf_return = sscanf(optarg, "%u", &setup_id);
      if (sscanf_return == EOF || sscanf_return == 0) {
        fprintf(stderr,
                "Please enter a positive integer for the setup id instead of "
                "\"-%c %s\"\n",
                optchar, optarg);
        setup_id = 0;
      }
      break;
    case 'i':
      sscanf_return = sscanf(optarg, "%zu", &num_iterations);
      if (sscanf_return == EOF || sscanf_return == 0) {
        fprintf(stderr,
                "Please enter a positive integer for the number of iterations "
                "instead of "
                "\"-%c %s\"\n",
                optchar, optarg);
        num_iterations = 0;
      }
      break;
    case 'c':
#if float_type == double
      sscanf_return = sscanf(optarg, "%lf", &Sc);
#else
      sscanf_return = sscanf(optarg, "%f", &Sc);
#endif
      if (sscanf_return == EOF || sscanf_return == 0 || Sc < float_cst(0.)) {
        fprintf(
            stderr,
            "Please enter a positive floating point number for the "
            "Courant-Friedrichs-Levy stability value instead of \"-%c %s\"\n",
            optchar, optarg);
        Sc = default_Sc;
      }
      break;
    case 't':
#if float_type == double
      sscanf_return = sscanf(optarg, "%lf", &end_time);
#else
      sscanf_return = sscanf(optarg, "%f", &end_time);
#endif
      if (sscanf_return == EOF || sscanf_return == 0 ||
          end_time < float_cst(0.)) {
        fprintf(stderr,
                "Please enter a positive floating point number for the "
                "simulation end time instead of \"-%c %s\"\n",
                optchar, optarg);
        end_time = default_end_time;
      }
      break;
    case 'w':
#if float_type == double
      sscanf_return = sscanf(optarg, "%lf", &smallest_wavelength);
#else
      sscanf_return = sscanf(optarg, "%f", &smallest_wavelength);
#endif
      if (sscanf_return == EOF || sscanf_return == 0 ||
          smallest_wavelength < float_cst(0.)) {
        fprintf(
            stderr,
            "Please enter a positive floating point number for the "
            "Courant-Friedrichs-Levy stability value instead of \"-%c %s\"\n",
            optchar, optarg);
        smallest_wavelength = default_smallest_wavelength;
      }
      break;
    case 'h':
      printf("Usage: %s <options>\n%s\n", argv[0], help_string);
      return EXIT_SUCCESS;
    }
  }

  unsigned initialize_setup_id;
  switch (dimension) {
  case 1: {
    if (setup_id >= last_1D_setup) {
      fprintf(stderr,
              "The input setup id %u does not map to any available 1D setup\n",
              setup_id);
      exit(EXIT_FAILURE);
    }
    initialize_setup_id = setup_id;
    if (Sc == default_Sc) {
      Sc = float_cst(1.);
    }
  } break;
  case 2: {
    if (setup_id >= last_2D_setup - last_1D_setup - 1) {
      fprintf(stderr,
              "The input setup id %u does not map to any available 2D setup\n",
              setup_id);
      exit(EXIT_FAILURE);
    }
    initialize_setup_id = setup_id + last_1D_setup + 1;
    if (Sc == default_Sc) {
      Sc = float_cst(1.) / sqrt(float_cst(3.));
    }
  } break;
  case 3: {
    if (setup_id >= last_3D_setup - last_2D_setup - 1) {
      fprintf(stderr,
              "The input setup id %u does not map to any available 3D setup\n",
              setup_id);
      exit(EXIT_FAILURE);
    }
    initialize_setup_id = setup_id + last_2D_setup + 1;
    if (Sc == default_Sc) {
      Sc = float_cst(1.) / sqrt(float_cst(4.));
    }
  } break;
  }

  struct fdtd fdtd =
      initializeFdtd_cmpl(initialize_setup_id, domain_size, Sc,
                          smallest_wavelength, border_cpml_width);

  float_type stop_time;
  if (end_time > float_cst(0.)) {
    stop_time = end_time;
  } else {
    stop_time = (float_type)num_iterations * get_time_step_fdtd(&fdtd);
  }
  time_measure startTime, endTime;
  get_current_time(&startTime);
  run_fdtd(&fdtd, stop_time, verbose);
  get_current_time(&endTime);
  fprintf(stdout, "Kernel time %.4fs\n",
          measuring_difftime(startTime, endTime));
  if (output_filename) {
    dump_fdtd(&fdtd, output_filename, dump_ez);
  }
  free_fdtd(&fdtd);

  return EXIT_SUCCESS;
}
