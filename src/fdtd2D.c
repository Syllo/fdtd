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

#include "fdtd2D.h"
#include "fdtd.h"
#include "fdtd_common.h"
#include "time_measurement.h"
#include <inttypes.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>

extern double rand_skip_percent;
extern bool sort_skip;
extern bool interpolate;

struct dataPosDeviation {
  double previous_val;
  double error;
  size_t posX, posY;
};

static int compareDeviation(const void *a, const void *b) {
  const struct dataPosDeviation *data1 = (const struct dataPosDeviation *)a;
  const struct dataPosDeviation *data2 = (const struct dataPosDeviation *)b;
  if (data1->error < data2->error) {
    return -1;
  } else {
    if (data1->error > data2->error) {
      return 1;
    } else {
      return 0;
    }
  }
}

static void update_electric_field(struct fdtd2D *fdtd) {
  VLA_2D_definition(float_type, fdtd->sizeX, fdtd->sizeY, hx, fdtd->hx);
  VLA_2D_definition(float_type, fdtd->sizeX, fdtd->sizeY, hy, fdtd->hy);
  VLA_2D_definition(float_type, fdtd->sizeX, fdtd->sizeY, ez, fdtd->ez);
  VLA_2D_definition(float_type, fdtd->sizeX, fdtd->sizeY, permittivity_inv,
                    fdtd->permittivity_inv);

  float_type _dx = float_cst(1.) / fdtd->dx;
  float_type _dy = float_cst(1.) / fdtd->dy;

  struct dataPosDeviation(*dpd)[fdtd->sizeY - 1] =
      malloc(sizeof(struct dataPosDeviation[fdtd->sizeX - 1][fdtd->sizeY - 1]));

  for (uintmax_t i = 1; i < fdtd->sizeX; ++i) {
    for (uintmax_t j = 1; j < fdtd->sizeY; ++j) {
      if (i > fdtd->cpml_thickness + 1 &&
          i < fdtd->sizeX - fdtd->cpml_thickness - 1 &&
          j > fdtd->cpml_thickness + 1 &&
          j < fdtd->sizeY - fdtd->cpml_thickness - 1 && !sort_skip &&
          !interpolate && drand48() < rand_skip_percent)
        continue;
      if (sort_skip) {
        dpd[i - 1][j - 1].previous_val = ez[i][j];
      }
      ez[i][j] = ez[i][j] + ((hy[i][j] - hy[i - 1][j]) * _dx -
                             (hx[i][j] - hx[i][j - 1]) * _dy) *
                                fdtd->dt * permittivity_inv[i][j];
      if (sort_skip) {
        dpd[i - 1][j - 1].error = dpd[i - 1][j - 1].previous_val - ez[i][j];
        dpd[i - 1][j - 1].error *= dpd[i - 1][j - 1].error;
        dpd[i - 1][j - 1].posX = i;
        dpd[i - 1][j - 1].posY = j;
      }
    }
  }
  if (interpolate) {
    for (uintmax_t i = fdtd->cpml_thickness + 1;
         i < fdtd->sizeX - fdtd->cpml_thickness - 1; ++i) {
      for (uintmax_t j = fdtd->cpml_thickness + 1;
           j < fdtd->sizeY - fdtd->cpml_thickness - 1; ++j) {
        if (drand48() < rand_skip_percent) {
          ez[i][j] =
              (ez[i - 1][j] + ez[i + 1][j] + ez[i][j - 1] + ez[i][j + 1]) / 4.;
        }
      }
    }
  }
  if (sort_skip) {
    qsort(dpd, (fdtd->sizeX - 1) * (fdtd->sizeY - 1),
          sizeof(struct dataPosDeviation), compareDeviation);
    double threshold_d =
        ceil(((fdtd->sizeX - 1) * (fdtd->sizeY - 1)) * rand_skip_percent);
    size_t threshold = (size_t)threshold_d;
    size_t whereAmI = 0;
    for (size_t i = 0; whereAmI < threshold && i < fdtd->sizeX - 1; ++i) {
      for (size_t j = 0; whereAmI < threshold && j < fdtd->sizeY - 1;
           ++j, whereAmI++) {
        // simulate skipping the computation for the values that have lowest
        // update derivative
        ez[dpd[i][j].posX][dpd[i][j].posY] = dpd[i][j].previous_val;
      }
    }
  }
  free(dpd);
}

static void update_electric_cpml(struct fdtd2D *fdtd) {
  VLA_2D_definition(float_type, fdtd->cpml_thickness, fdtd->sizeY, psi_ez_south,
                    fdtd->psi_ez[border_south]);
  VLA_2D_definition(float_type, fdtd->cpml_thickness, fdtd->sizeY, psi_ez_north,
                    fdtd->psi_ez[border_north]);
  VLA_2D_definition(float_type, fdtd->sizeX, fdtd->cpml_thickness, psi_ez_east,
                    fdtd->psi_ez[border_east]);
  VLA_2D_definition(float_type, fdtd->sizeX, fdtd->cpml_thickness, psi_ez_west,
                    fdtd->psi_ez[border_west]);
  VLA_2D_definition(float_type, fdtd->sizeX, fdtd->sizeY, hx, fdtd->hx);
  VLA_2D_definition(float_type, fdtd->sizeX, fdtd->sizeY, hy, fdtd->hy);
  VLA_2D_definition(float_type, fdtd->sizeX, fdtd->sizeY, ez, fdtd->ez);
  VLA_2D_definition(float_type, fdtd->sizeX, fdtd->sizeY, permittivity_inv,
                    fdtd->permittivity_inv);

  float_type _dy = float_cst(1.) / fdtd->dy;
  float_type _dx = float_cst(1.) / fdtd->dx;

  if (fdtd->border_condition[border_south] & border_cpml) { // ez_x
    for (uintmax_t i = 0; i < fdtd->cpml_thickness; ++i) {
      for (uintmax_t j = 0; j < fdtd->sizeY; ++j) {
        psi_ez_south[i][j] = fdtd->bx[i] * psi_ez_south[i][j] +
                             fdtd->cx[i] * (hy[1 + i][j] - hy[i][j]) * _dx;
        ez[1 + i][j] = ez[1 + i][j] + fdtd->dt * permittivity_inv[1 + i][j] *
                                          psi_ez_south[i][j];
      }
    }
  }
  if (fdtd->border_condition[border_north] & border_cpml) { // ez_x
    for (uintmax_t i = 0; i < fdtd->cpml_thickness; ++i) {
      for (uintmax_t j = 0; j < fdtd->sizeY; ++j) {
        psi_ez_north[i][j] =
            fdtd->bx[i] * psi_ez_north[i][j] +
            fdtd->cx[i] *
                (hy[fdtd->sizeX - 1 - i][j] - hy[fdtd->sizeX - 2 - i][j]) * _dx;
        ez[fdtd->sizeX - 1 - i][j] =
            ez[fdtd->sizeX - 1 - i][j] +
            fdtd->dt * permittivity_inv[fdtd->sizeX - 1 - i][j] *
                psi_ez_north[i][j];
      }
    }
  }
  if (fdtd->border_condition[border_west] & border_cpml) { // ez_y
    for (uintmax_t i = 0; i < fdtd->sizeX; ++i) {
      for (uintmax_t j = 0; j < fdtd->cpml_thickness; ++j) {
        psi_ez_west[i][j] = fdtd->by[j] * psi_ez_west[i][j] +
                            fdtd->cy[j] * (hx[i][1 + j] - hx[i][j]) * _dy;
        ez[i][1 + j] = ez[i][1 + j] - fdtd->dt * permittivity_inv[i][1 + j] *
                                          psi_ez_west[i][j];
      }
    }
  }
  if (fdtd->border_condition[border_east] & border_cpml) { // ez_y
    for (uintmax_t i = 0; i < fdtd->sizeX; ++i) {
      for (uintmax_t j = 0; j < fdtd->cpml_thickness; ++j) {
        psi_ez_east[i][j] =
            fdtd->by[j] * psi_ez_east[i][j] +
            fdtd->cy[j] *
                (hx[i][fdtd->sizeY - 1 - j] - hx[i][fdtd->sizeY - 2 - j]) * _dy;
        ez[i][fdtd->sizeY - 1 - j] =
            ez[i][fdtd->sizeY - 1 - j] -
            fdtd->dt * permittivity_inv[i][fdtd->sizeY - 1 - j] *
                psi_ez_east[i][j];
      }
    }
  }
}

static void border_condition_electric(struct fdtd2D *fdtd) {
  VLA_2D_definition(float_type, fdtd->sizeX, fdtd->sizeY, ez, fdtd->ez);
  for (enum border_position2D i = border_south; i < num_borders_2D; ++i) {
    if (fdtd->border_condition[i] & border_perfect_electric_conductor) {
      switch (i) {
      case border_south:
        for (uintmax_t j = 0; j < fdtd->sizeY; ++j)
          ez[0][j] = float_cst(0.);
        break;
      case border_north:
        for (uintmax_t j = 0; j < fdtd->sizeY; ++j)
          ez[fdtd->sizeX - 1][j] = float_cst(0.);
        break;
      case border_east:
        for (uintmax_t j = 0; j < fdtd->sizeX; ++j)
          ez[j][fdtd->sizeY - 1] = float_cst(0.);
        break;
      case border_west:
        for (uintmax_t j = 0; j < fdtd->sizeX; ++j)
          ez[j][0] = float_cst(0.);
        break;
      default:
        fprintf(stderr, "Error while processing the magnetic border "
                        "conditions: unknown border\n");
        exit(EXIT_FAILURE);
      }
    }
  }
}

static void update_magnetic_field(struct fdtd2D *fdtd) {
  VLA_2D_definition(float_type, fdtd->sizeX, fdtd->sizeY, hx, fdtd->hx);
  VLA_2D_definition(float_type, fdtd->sizeX, fdtd->sizeY, hy, fdtd->hy);
  VLA_2D_definition(float_type, fdtd->sizeX, fdtd->sizeY, ez, fdtd->ez);
  VLA_2D_definition(float_type, fdtd->sizeX, fdtd->sizeY, permeability_inv,
                    fdtd->permeability_inv);

  float_type _dy = float_cst(1.) / fdtd->dy;
  float_type _dx = float_cst(1.) / fdtd->dx;
  struct dataPosDeviation(*dpd)[fdtd->sizeY - 1] =
      malloc(sizeof(struct dataPosDeviation[fdtd->sizeX - 1][fdtd->sizeY - 1]));

  for (uintmax_t i = 0; i < fdtd->sizeX - 1; ++i) {
    for (uintmax_t j = 0; j < fdtd->sizeY - 1; ++j) {
      if (i > fdtd->cpml_thickness + 1 &&
          i < fdtd->sizeX - fdtd->cpml_thickness - 1 &&
          j > fdtd->cpml_thickness + 1 &&
          j < fdtd->sizeY - fdtd->cpml_thickness - 1 && !sort_skip &&
          !interpolate && drand48() < rand_skip_percent)
        continue;
      if (sort_skip) {
        dpd[i][j].previous_val = hx[i][j];
      }
      hx[i][j] = hx[i][j] + (ez[i][j] - ez[i][j + 1]) * _dy * fdtd->dt *
                                permeability_inv[i][j];
      if (sort_skip) {
        dpd[i][j].error = dpd[i][j].previous_val - hx[i][j];
        dpd[i][j].error *= dpd[i][j].error;
        dpd[i][j].posX = i;
        dpd[i][j].posY = j;
      }
    }
  }
  if (interpolate) {
    for (uintmax_t i = fdtd->cpml_thickness + 1;
         i < fdtd->sizeX - fdtd->cpml_thickness - 1; ++i) {
      for (uintmax_t j = fdtd->cpml_thickness + 1;
           j < fdtd->sizeY - fdtd->cpml_thickness - 1; ++j) {
        if (drand48() < rand_skip_percent) {
          hx[i][j] =
              (hx[i - 1][j] + hx[i + 1][j] + hx[i][j - 1] + hx[i][j + 1]) / 4.;
        }
      }
    }
  }
  if (sort_skip) {
    qsort(dpd, (fdtd->sizeX - 1) * (fdtd->sizeY - 1),
          sizeof(struct dataPosDeviation), compareDeviation);
    double threshold_d =
        ceil(((fdtd->sizeX - 1) * (fdtd->sizeY - 1)) * rand_skip_percent);
    size_t threshold = (size_t)threshold_d;
    size_t whereAmI = 0;
    for (size_t i = 0; whereAmI < threshold && i < fdtd->sizeX - 1; ++i) {
      for (size_t j = 0; whereAmI < threshold && j < fdtd->sizeY - 1;
           ++j, whereAmI++) {
        // simulate skipping the computation for the values that have lowest
        // update derivative
        hx[dpd[i][j].posX][dpd[i][j].posY] = dpd[i][j].previous_val;
      }
    }
  }
  for (uintmax_t i = 0; i < fdtd->sizeX - 1; ++i) {
    for (uintmax_t j = 0; j < fdtd->sizeY - 1; ++j) {
      if (i > fdtd->cpml_thickness + 1 &&
          i < fdtd->sizeX - fdtd->cpml_thickness - 1 &&
          j > fdtd->cpml_thickness + 1 &&
          j < fdtd->sizeY - fdtd->cpml_thickness - 1 && !sort_skip &&
          !interpolate && drand48() < rand_skip_percent)
        continue;
      if (sort_skip) {
        dpd[i][j].previous_val = hy[i][j];
      }
      hy[i][j] = hy[i][j] + (ez[i + 1][j] - ez[i][j]) * _dx * fdtd->dt *
                                permeability_inv[i][j];
      if (sort_skip) {
        dpd[i][j].error = dpd[i][j].previous_val - hy[i][j];
        dpd[i][j].error *= dpd[i][j].error;
        dpd[i][j].posX = i;
        dpd[i][j].posY = j;
      }
    }
  }
  if (interpolate) {
    for (uintmax_t i = fdtd->cpml_thickness + 1;
         i < fdtd->sizeX - fdtd->cpml_thickness - 1; ++i) {
      for (uintmax_t j = fdtd->cpml_thickness + 1;
           j < fdtd->sizeY - fdtd->cpml_thickness - 1; ++j) {
        if (drand48() < rand_skip_percent) {
          hy[i][j] =
              (hy[i - 1][j] + hy[i + 1][j] + hy[i][j - 1] + hy[i][j + 1]) / 4.;
        }
      }
    }
  }
  if (sort_skip) {
    qsort(dpd, (fdtd->sizeX - 1) * (fdtd->sizeY - 1),
          sizeof(struct dataPosDeviation), compareDeviation);
    double threshold_d =
        ceil(((fdtd->sizeX - 1) * (fdtd->sizeY - 1)) * rand_skip_percent);
    size_t threshold = (size_t)threshold_d;
    size_t whereAmI = 0;
    for (size_t i = 0; whereAmI < threshold && i < fdtd->sizeX - 1; ++i) {
      for (size_t j = 0; whereAmI < threshold && j < fdtd->sizeY - 1;
           ++j, whereAmI++) {
        // simulate skipping the computation for the values that have lowest
        // update derivative
        hy[dpd[i][j].posX][dpd[i][j].posY] = dpd[i][j].previous_val;
      }
    }
  }
  free(dpd);
}

static void update_magnetic_cpml(struct fdtd2D *fdtd) {
  VLA_2D_definition(float_type, fdtd->sizeX, fdtd->cpml_thickness, psi_hx_west,
                    fdtd->psi_hx_y[0]);
  VLA_2D_definition(float_type, fdtd->sizeX, fdtd->cpml_thickness, psi_hx_east,
                    fdtd->psi_hx_y[1]);
  VLA_2D_definition(float_type, fdtd->cpml_thickness, fdtd->sizeY, psi_hy_south,
                    fdtd->psi_hy_x[0]);
  VLA_2D_definition(float_type, fdtd->cpml_thickness, fdtd->sizeY, psi_hy_north,
                    fdtd->psi_hy_x[1]);
  VLA_2D_definition(float_type, fdtd->sizeX, fdtd->sizeY, hx, fdtd->hx);
  VLA_2D_definition(float_type, fdtd->sizeX, fdtd->sizeY, hy, fdtd->hy);
  VLA_2D_definition(float_type, fdtd->sizeX, fdtd->sizeY, ez, fdtd->ez);
  VLA_2D_definition(float_type, fdtd->sizeX, fdtd->sizeY, permeability_inv,
                    fdtd->permeability_inv);

  float_type _dy = float_cst(1.) / fdtd->dy;
  float_type _dx = float_cst(1.) / fdtd->dx;

  if (fdtd->border_condition[border_west] & border_cpml) { // hx_y
    for (uintmax_t i = 0; i < fdtd->sizeX; ++i) {
      for (uintmax_t j = 0; j < fdtd->cpml_thickness; ++j) {
        psi_hx_west[i][j] = fdtd->by[j] * psi_hx_west[i][j] +
                            fdtd->cy[j] * (ez[i][j] - ez[i][j + 1]) * _dy;
        hx[i][j] =
            hx[i][j] + fdtd->dt * permeability_inv[i][j] * psi_hx_west[i][j];
      }
    }
  }
  if (fdtd->border_condition[border_east] & border_cpml) { // hx_y
    for (uintmax_t i = 0; i < fdtd->sizeX; ++i) {
      for (uintmax_t j = 0; j < fdtd->cpml_thickness; ++j) {
        psi_hx_east[i][j] =
            fdtd->by[j] * psi_hx_east[i][j] +
            fdtd->cy[j] *
                (ez[i][fdtd->sizeY - 2 - j] - ez[i][fdtd->sizeY - 1 - j]) * _dy;
        hx[i][fdtd->sizeY - 2 - j] =
            hx[i][fdtd->sizeY - 2 - j] +
            fdtd->dt * permeability_inv[i][fdtd->sizeY - 2 - j] *
                psi_hx_east[i][j];
      }
    }
  }
  if (fdtd->border_condition[border_south] & border_cpml) { // hy_x
    for (uintmax_t i = 0; i < fdtd->cpml_thickness; ++i) {
      for (uintmax_t j = 0; j < fdtd->sizeY; ++j) {
        psi_hy_south[i][j] = fdtd->bx[i] * psi_hy_south[i][j] +
                             fdtd->cx[i] * (ez[i + 1][j] - ez[i][j]) * _dx;
        hy[i][j] =
            hy[i][j] + fdtd->dt * permeability_inv[i][j] * psi_hy_south[i][j];
      }
    }
  }
  if (fdtd->border_condition[border_north] & border_cpml) { // hy_x
    for (uintmax_t i = 0; i < fdtd->cpml_thickness; ++i) {
      for (uintmax_t j = 0; j < fdtd->sizeY; ++j) {
        psi_hy_north[i][j] =
            fdtd->bx[i] * psi_hy_north[i][j] +
            fdtd->cx[i] *
                (ez[fdtd->sizeX - 1 - i][j] - ez[fdtd->sizeX - 2 - i][j]) * _dx;
        hy[fdtd->sizeX - 2 - i][j] =
            hy[fdtd->sizeX - 2 - i][j] +
            fdtd->dt * permeability_inv[fdtd->sizeX - 2 - i][j] *
                psi_hy_north[i][j];
      }
    }
  }
}

static void border_condition_magnetic(struct fdtd2D *fdtd) {
  VLA_2D_definition(float_type, fdtd->sizeX, fdtd->sizeY, hx, fdtd->hx);
  VLA_2D_definition(float_type, fdtd->sizeX, fdtd->sizeY, hy, fdtd->hy);
  for (enum border_position2D i = border_south; i < num_borders_2D; ++i) {
    if (fdtd->border_condition[i] & border_perfect_electric_conductor) {
      switch (i) {
      case border_south:
        for (uintmax_t j = 0; j < fdtd->sizeY; ++j) {
          hx[0][j] = float_cst(0.);
          hy[0][j] = float_cst(0.);
        }
        break;
      case border_north:
        for (uintmax_t j = 0; j < fdtd->sizeY; ++j) {
          hx[fdtd->sizeX - 2][j] = float_cst(0.);
          hy[fdtd->sizeX - 2][j] = float_cst(0.);
          hx[fdtd->sizeX - 1][j] = float_cst(0.);
          hy[fdtd->sizeX - 1][j] = float_cst(0.);
        }
        break;
      case border_east:
        for (uintmax_t j = 0; j < fdtd->sizeX; ++j) {
          hx[j][fdtd->sizeY - 2] = float_cst(0.);
          hy[j][fdtd->sizeY - 2] = float_cst(0.);
          hx[j][fdtd->sizeY - 1] = float_cst(0.);
          hy[j][fdtd->sizeY - 1] = float_cst(0.);
        }
        break;
      case border_west:
        for (uintmax_t j = 0; j < fdtd->sizeX; ++j) {
          hx[j][0] = float_cst(0.);
          hy[j][0] = float_cst(0.);
        }
        break;
      default:
        fprintf(stderr, "Error while processing the electric border "
                        "conditions: unknown border\n");
        exit(EXIT_FAILURE);
      }
    }
  }
}

static void apply_M_sources(struct fdtd2D *fdtd) {
  VLA_2D_definition(float_type, fdtd->sizeX, fdtd->sizeY, hx, fdtd->hx);
  VLA_2D_definition(float_type, fdtd->sizeX, fdtd->sizeY, hy, fdtd->hy);
  VLA_2D_definition(uintmax_t, fdtd->num_Msources, 2, MsourceLocations,
                    fdtd->MsourceLocations);
  for (uintmax_t i = 0; i < fdtd->num_Msources; ++i) {
    switch (fdtd->Msources[i].type) {
    case source_gaussian_pulse:
      hx[MsourceLocations[i][0]][MsourceLocations[i][1]] +=
          gaussian_pulse_val(fdtd->time, &fdtd->Msources[i]);
      hy[MsourceLocations[i][0]][MsourceLocations[i][1]] +=
          gaussian_pulse_val(fdtd->time, &fdtd->Msources[i]);
      break;
    default:
      break;
    }
  }
}

static void apply_J_sources(struct fdtd2D *fdtd) {
  VLA_2D_definition(float_type, fdtd->sizeX, fdtd->sizeY, ez, fdtd->ez);
  VLA_2D_definition(uintmax_t, fdtd->num_Jsources, 2, JsourceLocations,
                    fdtd->JsourceLocations);
  for (uintmax_t i = 0; i < fdtd->num_Jsources; ++i) {
    switch (fdtd->Jsources[i].type) {
    case source_gaussian_pulse:
      ez[JsourceLocations[i][0]][JsourceLocations[i][1]] -=
          gaussian_pulse_val(fdtd->time, &fdtd->Jsources[i]);
      break;
    default:
      break;
    }
  }
}

void init_fdtd_2D_medium(struct fdtd2D *fdtd,
                         init_medium_fun_2D permeability_invR,
                         init_medium_fun_2D permittivity_invR, void *user) {
  VLA_2D_definition(float_type, fdtd->sizeX, fdtd->sizeY, permittivity_inv,
                    fdtd->permittivity_inv);
  VLA_2D_definition(float_type, fdtd->sizeX, fdtd->sizeY, permeability_inv,
                    fdtd->permeability_inv);
  float_type posX = float_cst(0.);
  for (uintmax_t i = 0; i < fdtd->sizeX; ++i, posX += fdtd->dx) {
    float_type posY = float_cst(0.);
    for (uintmax_t j = 0; j < fdtd->sizeY; ++j, posY += fdtd->dy) {
      permeability_inv[i][j] =
          float_cst(1.) / (permeability_invR(posX, posY, user) * mu0);
      permittivity_inv[i][j] =
          float_cst(1.) / (permittivity_invR(posX, posY, user) * eps0);
    }
  }
}

struct fdtd2D init_fdtd_2D(float_type domain_size[2], float_type Sc,
                           float_type smallest_wavelength,
                           enum border_condition borders[num_borders_2D]) {
  return init_fdtd_2D_cpml(domain_size, Sc, smallest_wavelength, borders, 0);
}

struct fdtd2D init_fdtd_2D_cpml(float_type domain_size[2], float_type Sc,
                                float_type smallest_wavelength,
                                enum border_condition borders[num_borders_2D],
                                uintmax_t cpml_thickness) {

  float_type dx = smallest_wavelength / float_cst(20.);
  float_type dy = dx;
  float_type dt = dx * Sc / c_light;
  float_type sizeXf = floor(domain_size[0] / dx);
  float_type sizeYf = floor(domain_size[1] / dy);
  uintmax_t sizeX = (uintmax_t)sizeXf;
  uintmax_t sizeY = (uintmax_t)sizeYf;
  float_type Sc_max = float_cst(1.) / sqrt(float_cst(2.));
  if (Sc > Sc_max)
    fprintf(stderr,
            "The value of Sc is too high, the simulation may be unstable. "
            "Please use a value lesser or equal to %.5f\n",
            Sc_max);

  struct fdtd2D fdtd = {
      .dx = dx,
      .dy = dy,
      .dt = dt,
      .ez = calloc(1, VLA_2D_size(float_type, sizeX, sizeY)),
      .hx = calloc(1, VLA_2D_size(float_type, sizeX, sizeY)),
      .hy = calloc(1, VLA_2D_size(float_type, sizeX, sizeY)),
      .permeability_inv = calloc(1, VLA_2D_size(float_type, sizeX, sizeY)),
      .permittivity_inv = calloc(1, VLA_2D_size(float_type, sizeX, sizeY)),
      .psi_hx_y = {NULL, NULL},
      .psi_hy_x = {NULL, NULL},
      .psi_ez = {NULL, NULL, NULL, NULL},
      .bx = NULL,
      .by = NULL,
      .cx = NULL,
      .cy = NULL,
      .cpml_thickness = cpml_thickness,
      .border_condition = {[border_south] = borders[border_south],
                           [border_north] = borders[border_north],
                           [border_east] = borders[border_east],
                           [border_west] = borders[border_west]},
      .domain_size = {domain_size[0], domain_size[1]},
      .sizeX = sizeX,
      .sizeY = sizeY,
      .Sc = Sc,
      .num_Jsources = 0,
      .Jsources = NULL,
      .JsourceLocations = NULL,
      .num_Msources = 0,
      .Msources = NULL,
      .MsourceLocations = NULL,
      .time = float_cst(0.),
  };
  if (cpml_thickness > 0 && fdtd.border_condition[border_south] & border_cpml) {
    fdtd.psi_hy_x[0] =
        calloc(1, VLA_2D_size(float_type, cpml_thickness, sizeY));
    fdtd.psi_ez[border_south] =
        calloc(1, VLA_2D_size(float_type, cpml_thickness, sizeY));
  }
  if (cpml_thickness > 0 && fdtd.border_condition[border_north] & border_cpml) {
    fdtd.psi_hy_x[1] =
        calloc(1, VLA_2D_size(float_type, cpml_thickness, sizeY));
    fdtd.psi_ez[border_north] =
        calloc(1, VLA_2D_size(float_type, cpml_thickness, sizeY));
  }
  if (cpml_thickness > 0 && fdtd.border_condition[border_west] & border_cpml) {
    fdtd.psi_hx_y[0] =
        calloc(1, VLA_2D_size(float_type, sizeX, cpml_thickness));
    fdtd.psi_ez[border_west] =
        calloc(1, VLA_2D_size(float_type, sizeX, cpml_thickness));
  }
  if (cpml_thickness > 0 && fdtd.border_condition[border_east] & border_cpml) {
    fdtd.psi_hx_y[1] =
        calloc(1, VLA_2D_size(float_type, sizeX, cpml_thickness));
    fdtd.psi_ez[border_east] =
        calloc(1, VLA_2D_size(float_type, sizeX, cpml_thickness));
  }
  fdtd.bx = calloc(cpml_thickness, sizeof(*fdtd.bx));
  fdtd.by = fdtd.bx;
  fdtd.cx = calloc(cpml_thickness, sizeof(*fdtd.bx));
  fdtd.cy = fdtd.cx;
  float_type alpha_max = float_cst(2.) * M_PI * eps0 * fdtd.dx * float_cst(0.1);
  float_type sigma_max = float_cst(0.8) * (polynomial_taper_order + 1) /
                         (fdtd.dx * sqrt(mu0 / eps0));
  for (uintmax_t d = 0; d < cpml_thickness; ++d) {
    // store from border of domain to interface CPML / simulation medium
    fdtd.bx[cpml_thickness - d - 1] =
        b(d, cpml_thickness - 1, fdtd.dt, alpha_max, sigma_max);
    fdtd.cx[cpml_thickness - d - 1] =
        c(d, cpml_thickness - 1, fdtd.dt, alpha_max, sigma_max);
  }
  fprintf(stderr, "Dt %e Dx %e Dy %e (%.0fx%.0f)\n", dt, dx, dy, sizeXf,
          sizeYf);

  return fdtd;
}

void run_2D_fdtd(struct fdtd2D *fdtd, float_type end_time, bool verbose) {
  const double num_iter_d = ceil((end_time - fdtd->time) / fdtd->dt);
  double print_interval_d;
  if (num_iter_d >= 10.) {
    double divide = 1.;
    do {
      print_interval_d = ceil(num_iter_d / divide);
      divide = divide + 1.;
    } while (num_iter_d / print_interval_d < 10.);
  } else {
    print_interval_d = 1.;
  }
  const size_t print_interval = (size_t)print_interval_d;
  const double percent_increment = 100. / (num_iter_d / print_interval_d);
  const size_t inter_print = print_interval - 1;
  size_t iter_count = 0;
  double percentage = percent_increment;
  time_measure tstart_chunk, tend_chunk;
  get_current_time(&tstart_chunk);
  for (; fdtd->time < end_time; fdtd->time += fdtd->dt) {
    update_magnetic_field(fdtd);
    apply_M_sources(fdtd);
    update_magnetic_cpml(fdtd);
    border_condition_magnetic(fdtd);

    update_electric_field(fdtd);
    apply_J_sources(fdtd);
    update_electric_cpml(fdtd);
    border_condition_electric(fdtd);

    iter_count = iter_count == inter_print ? 0 : iter_count + 1;
    if (verbose && iter_count == 0) {
      get_current_time(&tend_chunk);
      double difference = measuring_difftime(tstart_chunk, tend_chunk);
      printf("%.0f%% -- t=%e dt=%e tend=%e (%zu iter in %.3fs)\n", percentage,
             fdtd->time, fdtd->dt, end_time, print_interval, difference);
      percentage += percent_increment;
      tstart_chunk = tend_chunk;
    }
  }
}

void dump_2D_fdtd(const struct fdtd2D *fdtd, const char *fileName,
                  enum dumpable_data what_to_dump) {
  FILE *out = fopen(fileName, "w");
  VLA_2D_definition(float_type, fdtd->sizeX, fdtd->sizeY, data, fdtd->ez);
  switch (what_to_dump) {
  case dump_ez:
    break;
  case dump_hx:
    data = fdtd->hx;
    break;
  case dump_hy:
    data = fdtd->hy;
    break;
  case dump_permeability:
    data = fdtd->permeability_inv;
    break;
  case dump_permittivity:
    data = fdtd->permittivity_inv;
    break;
  default:
    fprintf(stderr, "Dump of \"%s\" not available for 2D fdtd\n",
            dumpable_data_name[what_to_dump]);
    exit(EXIT_FAILURE);
  }
  for (uintmax_t i = 0; i < fdtd->sizeX; ++i) {
    for (uintmax_t j = 1; j < fdtd->sizeY; ++j) {
      fprintf(out, "%e %e %e\n", (float_type)i * fdtd->dx,
              (float_type)j * fdtd->dy, data[i][j]);
    }
  }
  fclose(out);
}

void free_2D_fdtd(struct fdtd2D *fdtd) {
  free(fdtd->ez);
  free(fdtd->hx);
  free(fdtd->hy);
  free(fdtd->permittivity_inv);
  free(fdtd->permeability_inv);
  free(fdtd->JsourceLocations);
  free(fdtd->Jsources);
  free(fdtd->MsourceLocations);
  free(fdtd->Msources);
  free(fdtd->psi_hx_y[0]);
  free(fdtd->psi_hx_y[1]);
  free(fdtd->psi_hy_x[0]);
  free(fdtd->psi_hy_x[1]);
  for (enum border_position2D bd = 0; bd < num_borders_2D; ++bd) {
    free(fdtd->psi_ez[bd]);
  }
  free(fdtd->bx);
  free(fdtd->cx);
}

void add_source_fdtd_2D(enum source_type sType, struct fdtd2D *fdtd,
                        struct fdtd_source src, float_type positionX,
                        float_type positionY) {
  float_type posXf = floor(positionX / fdtd->dx);
  float_type posYf = floor(positionY / fdtd->dy);
  uintmax_t posX = (uintmax_t)posXf;
  uintmax_t posY = (uintmax_t)posYf;
  if (posX >= fdtd->sizeX) {
    fprintf(stderr, "add_source_fdtd_2D: adding source outside of the X "
                    "dimension boundaries\n");
    exit(EXIT_FAILURE);
  }
  if (posY >= fdtd->sizeY) {
    fprintf(stderr, "add_source_fdtd_2D: adding source outside of the Y "
                    "dimension boundaries\n");
    exit(EXIT_FAILURE);
  }
  switch (sType) {
  case source_electric: {
    fdtd->num_Jsources++;
    fdtd->JsourceLocations = realloc(
        fdtd->JsourceLocations, VLA_2D_size(uintmax_t, fdtd->num_Jsources, 2));
    VLA_2D_definition(uintmax_t, fdtd->num_Jsources, 2, JsourceLocations,
                      fdtd->JsourceLocations);
    JsourceLocations[fdtd->num_Jsources - 1][0] = posX;
    JsourceLocations[fdtd->num_Jsources - 1][1] = posY;
    fdtd->Jsources =
        realloc(fdtd->Jsources, fdtd->num_Jsources * sizeof(*fdtd->Jsources));
    fdtd->Jsources[fdtd->num_Jsources - 1] = src;
    break;
  }
  case source_magnetic:
    fdtd->num_Msources++;
    fdtd->MsourceLocations = realloc(
        fdtd->MsourceLocations, VLA_2D_size(uintmax_t, fdtd->num_Msources, 2));
    VLA_2D_definition(uintmax_t, fdtd->num_Msources, 2, MsourceLocations,
                      fdtd->MsourceLocations);
    MsourceLocations[fdtd->num_Msources - 1][0] = posX;
    MsourceLocations[fdtd->num_Msources - 1][1] = posY;
    fdtd->Msources =
        realloc(fdtd->Msources, fdtd->num_Msources * sizeof(*fdtd->Msources));
    fdtd->Msources[fdtd->num_Msources - 1] = src;
    break;
  }
}
