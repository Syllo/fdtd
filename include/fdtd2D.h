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

#ifndef FDTD2D_H_
#define FDTD2D_H_

#include <stdbool.h>
#include <stdint.h>

#include "fdtd_common.h"

#define arrayOffset2D(sizex, sizey, x, y) ((sizey * x) + y)
#define VLA_2D_definition(type, size1, size2, name, ptr)                       \
  type(*restrict name)[size2] = ptr;

#define VLA_2D_size(type, size1, size2) sizeof(type[size1][size2])

enum border_position2D {
  border_south = 0,
  border_north,
  border_east,
  border_west,
  num_borders_2D
};

struct fdtd2D {
  const float_type dx;    // Space step
  const float_type dy;    // Space step
  const float_type dt;    // Time step
  void *ez;               // Electric Field
  void *hx;               // Magnetic field
  void *hy;               // Magnetic field
  void *permittivity_inv; // 1 / Permittivity
  void *permeability_inv; // 1 / Permeability
  // The psi are discrete unknowns used to update the fields e and h with CPML
  // absorbing boundaries
  void *psi_hx_y[2];            // hx psi boundary normal to y (east & west)
  void *psi_hy_x[2];            // hy psi boundary normal to x (north & south)
  void *psi_ez[num_borders_2D]; // ez psi boundary normal to x and y
  // CPML discrete unknown for b and c are CPML constants that are used to
  // update psi
  float_type *bx, *by;
  float_type *cx, *cy;
  const uintmax_t cpml_thickness; // Absorbing CPML border thickness
  const enum border_condition
      border_condition[num_borders_2D]; // Border condition
  const float_type domain_size[2];      // Physical domain size
  const uintmax_t sizeX;                // Domain size
  const uintmax_t sizeY;                // Domain size
  const float_type Sc;                  // Courrant number
  unsigned num_Jsources;                // Count of sources
  struct fdtd_source *Jsources;         // Electric Sources
  void *JsourceLocations;               // Location of the Electric sources
  unsigned num_Msources;                // Count of sources
  struct fdtd_source *Msources;         // Magnetic Sources
  void *MsourceLocations;               // Location of the Magnetic sources
  float_type time;                      // Simulation current time
};

struct fdtd2D init_fdtd_2D(float_type domain_size[2], float_type Sc,
                           float_type smallest_wavelength,
                           enum border_condition borders[num_borders_2D]);

struct fdtd2D init_fdtd_2D_cpml(float_type domain_size[2], float_type Sc,
                                float_type smallest_wavelength,
                                enum border_condition borders[num_borders_2D],
                                uintmax_t cpml_thickness);

typedef float_type (*init_medium_fun_2D)(float_type, float_type, void *);

void init_fdtd_2D_medium(struct fdtd2D *fdtd,
                         init_medium_fun_2D permeability_revR,
                         init_medium_fun_2D permittivity_invR, void *user);

void run_2D_fdtd(struct fdtd2D *fdtd, float_type end);

void dump_2D_fdtd(const struct fdtd2D *fdtd, const char *fileName,
                  enum dumpable_data what_to_dump);

struct fdtd_source gaussian_source(float_type delay, float_type peak_time,
                                   float_type peak_val);

void free_2D_fdtd(struct fdtd2D *fdtd);

void add_source_fdtd_2D(enum source_type sType, struct fdtd2D *fdtd,
                        struct fdtd_source src, float_type positionX,
                        float_type positionY);

#endif // FDTD2D_H_
