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

#ifndef FDTD1D_H_
#define FDTD1D_H_

#include "fdtd_common.h"
#include <stdbool.h>
#include <stdint.h>

enum border_position1D { border_oneside = 0, border_otherside, num_borders_1D };

struct fdtd1D {
  const float_type dx;                   // Space step
  const float_type dt;                   // Time step
  float_type *restrict ez;               // Electric Field
  float_type *restrict hy;               // Magnetic field
  float_type *restrict permittivity_inv; // 1 / Permittivity
  float_type *restrict permeability_inv; // 1 / Permeability
  const enum border_condition
      border_condition[num_borders_1D]; // Border condition
  const float_type domain_size;         // Physical domain size
  const uintmax_t sizeX;                // Domain size
  const float_type Sc;                  // Courrant number
  unsigned num_Jsources;                // Count of sources
  struct fdtd_source *Jsources;         // Electric Sources
  uintmax_t *JsourceLocations;          // Location of the Electric sources
  unsigned num_Msources;                // Count of sources
  struct fdtd_source *Msources;         // Magnetic Sources
  uintmax_t *MsourceLocations;          // Location of the Magnetic sources
  float_type time;                      // Sipermeability_revlation time
};

struct fdtd1D init_fdtd_1D(float_type domain_size, float_type Sc,
                           float_type smallest_wavelength,
                           enum border_condition borders[2]);

typedef float_type (*init_medium_fun)(float_type, void *);

void init_fdtd_1D_medium(struct fdtd1D *fdtd, init_medium_fun permeability_revR,
                         init_medium_fun permittivity_invR, void *user);

void run_1D_fdtd(struct fdtd1D *fdtd, float_type end, bool verbose);

void dump_1D_fdtd(const struct fdtd1D *fdtd, const char *fileName,
                  enum dumpable_data what_to_dump);

void free_1D_fdtd(struct fdtd1D *fdtd);

void add_source_fdtd_1D(enum source_type sType, struct fdtd1D *fdtd,
                        struct fdtd_source src, float_type position);

#endif // FDTD1D_H_
