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

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>

#include "fdtd.h"
#include "fdtd1D.h"
#include "fdtd_common.h"

static void update_electric_field(struct fdtd1D *fdtd) {
  float_type dtdx = fdtd->dt / fdtd->dx;

  for (uintmax_t i = 1; i < fdtd->sizeX; ++i) {
    fdtd->ez[i] = fdtd->ez[i] + (fdtd->hy[i] - fdtd->hy[i - 1]) * dtdx *
                                    fdtd->permittivity_inv[i];
  }
}

static void border_condition_electric(struct fdtd1D *fdtd) {
  for (enum border_position1D i = border_oneside; i < num_borders_1D; ++i) {
    if (fdtd->border_condition[i] == border_perfect_electric_conductor) {
      uintmax_t pos = i * (fdtd->sizeX - 1);
      fdtd->ez[pos] = float_cst(0.);
    }
  }
}

static void update_magnetic_field(struct fdtd1D *fdtd) {
  float_type dtdx = fdtd->dt / fdtd->dx;

  for (uintmax_t i = 0; i < fdtd->sizeX - 1; ++i) {
    fdtd->hy[i] = fdtd->hy[i] + (fdtd->ez[i + 1] - fdtd->ez[i]) * dtdx *
                                    fdtd->permeability_inv[i];
  }
}

static void border_condition_magnetic(struct fdtd1D *fdtd) {
  for (enum border_position1D i = border_oneside; i < num_borders_1D; ++i) {
    if (fdtd->border_condition[i] == border_perfect_magnetic_conductor) {
      uintmax_t pos = i * (fdtd->sizeX - 1);
      fdtd->hy[pos] = float_cst(0.);
    }
  }
}

static void apply_M_sources(struct fdtd1D *fdtd) {
  for (uintmax_t i = 0; i < fdtd->num_Msources; ++i) {
    switch (fdtd->Msources[i].type) {
    case source_gaussian_pulse:
      fdtd->hy[fdtd->MsourceLocations[i]] +=
          gaussian_pulse_val(fdtd->time, &fdtd->Msources[i]);
      break;
    default:
      break;
    }
  }
}

static void apply_J_sources(struct fdtd1D *fdtd) {
  for (uintmax_t i = 0; i < fdtd->num_Jsources; ++i) {
    switch (fdtd->Jsources[i].type) {
    case source_gaussian_pulse:
      fdtd->ez[fdtd->JsourceLocations[i]] +=
          gaussian_pulse_val(fdtd->time, &fdtd->Jsources[i]);
      break;
    default:
      break;
    }
  }
}

void init_fdtd_1D_medium(struct fdtd1D *fdtd, init_medium_fun permeability_invR,
                         init_medium_fun permittivity_invR, void *user) {
  float_type pos = float_cst(0.);
  for (uintmax_t i = 0; i < fdtd->sizeX; ++i, pos += fdtd->dx) {
    fdtd->permeability_inv[i] =
        float_cst(1.) / (permeability_invR(pos, user) * mu0);
    fdtd->permittivity_inv[i] =
        float_cst(1.) / (permittivity_invR(pos, user) * eps0);
  }
}

// Permittivity of metal as free space
struct fdtd1D init_fdtd_1D(float_type domain_size, float_type Sc,
                           float_type smallest_wavelength,
                           enum border_condition borders[num_borders_1D]) {

  float_type dx = smallest_wavelength / float_cst(20.);
  float_type dt = dx * Sc / c_light;
  float_type sizeXf = ceil(domain_size / dx);
  uintmax_t sizeX = (uintmax_t)sizeXf;
  struct fdtd1D fdtd = {
      .dx = dx,
      .dt = dt,
      .ez = calloc(sizeX, sizeof(float_type)),
      .hy = calloc(sizeX, sizeof(float_type)),
      .permittivity_inv = calloc(sizeX, sizeof(float_type)),
      .permeability_inv = calloc(sizeX, sizeof(float_type)),
      .border_condition = {[border_oneside] = borders[border_oneside],
                           [border_otherside] = borders[border_otherside]},
      .domain_size = domain_size,
      .sizeX = sizeX,
      .Sc = Sc,
      .num_Jsources = 0,
      .Jsources = NULL,
      .JsourceLocations = NULL,
      .num_Msources = 0,
      .Msources = NULL,
      .MsourceLocations = NULL,
      .time = 0,
  };
  fprintf(stderr, "Dt %e Dx %e (%.0f)\n", dt, dx, sizeXf);

  return fdtd;
}

void run_1D_fdtd(struct fdtd1D *fdtd, float_type end_time) {
  fprintf(stderr, "It will take %.0f iterations\n",
          (end_time - fdtd->time) / fdtd->dt);
  /*exit(1);*/
  for (; fdtd->time < end_time; fdtd->time += fdtd->dt) {

    update_magnetic_field(fdtd);
    apply_M_sources(fdtd);
    border_condition_magnetic(fdtd);

    update_electric_field(fdtd);
    apply_J_sources(fdtd);
    border_condition_electric(fdtd);

    /*fprintf(stderr, "Time %e/%e\n", fdtd->time, end_time);*/
  }
}

void dump_1D_fdtd(const struct fdtd1D *fdtd, const char *fileName,
                  enum dumpable_data what_to_dump) {
  FILE *out = fopen(fileName, "w");
  float_type *array;
  switch (what_to_dump) {
  case dump_ez:
    array = fdtd->ez;
    break;
  case dump_hy:
    array = fdtd->hy;
    break;
  case dump_permeability:
    array = fdtd->permeability_inv;
    break;
  case dump_permittivity:
    array = fdtd->permittivity_inv;
    break;
  default:
    fprintf(stderr, "Dump of \"%s\" not available for 1D fdtd\n",
            dumpable_data_name[what_to_dump]);
    exit(EXIT_FAILURE);
  }
  for (uintmax_t i = 0; i < fdtd->sizeX; ++i) {
    fprintf(out, "%e %e\n", (float_type)i * fdtd->dx, array[i]);
  }
  fclose(out);
}

void free_1D_fdtd(struct fdtd1D *fdtd) {
  free(fdtd->ez);
  free(fdtd->hy);
  free(fdtd->permittivity_inv);
  free(fdtd->permeability_inv);
  free(fdtd->JsourceLocations);
  free(fdtd->Jsources);
  free(fdtd->MsourceLocations);
  free(fdtd->Msources);
}

void add_source_fdtd_1D(enum source_type sType, struct fdtd1D *fdtd,
                        struct fdtd_source src, float_type position) {
  float_type posX = ceil(position / fdtd->dx);
  switch (sType) {
  case source_electric:
    fdtd->num_Jsources++;
    fdtd->JsourceLocations =
        realloc(fdtd->JsourceLocations,
                fdtd->num_Jsources * sizeof(*fdtd->JsourceLocations));
    fdtd->JsourceLocations[fdtd->num_Jsources - 1] = (uintmax_t)posX;
    fdtd->Jsources =
        realloc(fdtd->Jsources, fdtd->num_Jsources * sizeof(*fdtd->Jsources));
    fdtd->Jsources[fdtd->num_Jsources - 1] = src;
    break;
  case source_magnetic:
    fdtd->num_Msources++;
    fdtd->MsourceLocations =
        realloc(fdtd->MsourceLocations,
                fdtd->num_Msources * sizeof(*fdtd->MsourceLocations));
    fdtd->MsourceLocations[fdtd->num_Msources - 1] = (uintmax_t)posX;
    fdtd->Msources =
        realloc(fdtd->Msources, fdtd->num_Msources * sizeof(*fdtd->Msources));
    fdtd->Msources[fdtd->num_Msources - 1] = src;
    break;
  }
}
