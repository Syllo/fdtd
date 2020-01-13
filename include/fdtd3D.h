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

#ifndef FDTD3D_H_
#define FDTD3D_H_

#include "fdtd_common.h"
#include <stdbool.h>
#include <stdint.h>

float_type gaussian_pulse_val(float_type time, struct fdtd_source *src);

#define VLA_3D_definition(type, size1, size2, size3, name, ptr)                \
  type(*restrict name)[size2][size3] = ptr;
#define VLA_3D_size(type, size1, size2, size3) sizeof(type[size1][size2][size3])

enum border_position3D {
  border_front = 0,
  border_back,
  border_top,
  border_bottom,
  border_right,
  border_left,
  num_borders_3D,
};

struct fdtd3D {
  const float_type dx;    // Space step
  const float_type dy;    // Space step
  const float_type dz;    // Space step
  const float_type dt;    // Time step
  void *hx;               // Magnetic field
  void *hy;               // Magnetic field
  void *hz;               // Magnetic field
  void *ex;               // Electric Field
  void *ey;               // Electric Field
  void *ez;               // Electric Field
  void *permittivity_inv; // 1 / Permittivity
  void *permeability_inv; // 1 / Permeability
  // The psi are discrete unknowns used to update the fields e and h with CPML
  // absorbing boundaries
  void *psi_hx_y[2]; // hx psi boundary normal to y (left & right)
  void *psi_hx_z[2]; // hx psi boundary normal to z (front & back)
  void *psi_hy_x[2]; // hy psi boundary normal to x (bottom & top)
  void *psi_hy_z[2]; // hy psi boundary normal to z (front & back)
  void *psi_hz_x[2]; // hz psi boundary normal to x (bottom & top)
  void *psi_hz_y[2]; // hz psi boundary normal to y (left & right)
  void *psi_ex_y[2]; // ex psi boundary normal to y (left & right)
  void *psi_ex_z[2]; // ex psi boundary normal to z (front & back)
  void *psi_ey_x[2]; // ey psi boundary normal to x (bottom & top)
  void *psi_ey_z[2]; // ey psi boundary normal to z (front & back)
  void *psi_ez_x[2]; // ez psi boundary normal to x (bottom & top)
  void *psi_ez_y[2]; // ez psi boundary normal to y (left & right)
  // CPML discrete unknown for b and c are CPML constants that are used to
  // update psi
  float_type *restrict bx, *restrict by, *restrict bz;
  float_type *restrict cx, *restrict cy, *restrict cz;
  const uintmax_t cpml_thickness; // Absorbing CPML border thickness
  const enum border_condition
      border_condition[num_borders_3D]; // Border condition
  const float_type domain_size[3];      // Physical domain size
  const uintmax_t sizeX;                // Domain size
  const uintmax_t sizeY;                // Domain size
  const uintmax_t sizeZ;                // Domain size
  const float_type Sc;                  // Courrant number
  unsigned num_Jsources;                // Count of sources
  struct fdtd_source *Jsources;         // Electric Sources
  void *JsourceLocations;               // Location of the Electric sources
  unsigned num_Msources;                // Count of sources
  struct fdtd_source *Msources;         // Magnetic Sources
  void *MsourceLocations;               // Location of the Magnetic sources
  float_type time;                      // Simulation current time
};

struct fdtd3D init_fdtd_3D(float_type domain_size[3], float_type Sc,
                           float_type smallest_wavelength,
                           enum border_condition borders[num_borders_3D]);

struct fdtd3D init_fdtd_3D_cpml(float_type domain_size[3], float_type Sc,
                                float_type smallest_wavelength,
                                enum border_condition borders[num_borders_3D],
                                uintmax_t cpml_thickness);

typedef float_type (*init_medium_fun_3D)(float_type, float_type, float_type,
                                         void *);

void init_fdtd_3D_medium(struct fdtd3D *fdtd,
                         init_medium_fun_3D permeability_invR,
                         init_medium_fun_3D permittivity_invR, void *user);

void run_3D_fdtd(struct fdtd3D *fdtd, float_type end, bool verbose);

void dump_3D_fdtd(const struct fdtd3D *fdtd, const char *fileName,
                  enum dumpable_data what_to_dump);

void free_3D_fdtd(struct fdtd3D *fdtd);

void add_source_fdtd_3D(enum source_type sType, struct fdtd3D *fdtd,
                        struct fdtd_source src, float_type positionX,
                        float_type positionY, float_type positionZ);

#endif // FDTD3D_H_
