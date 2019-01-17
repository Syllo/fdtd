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
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>

#include "fdtd.h"
#include "fdtd3D.h"
#include "fdtd_common.h"

static void update_electric_field(struct fdtd3D *fdtd) {
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, hx,
                    fdtd->hx);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, hy,
                    fdtd->hy);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, hz,
                    fdtd->hz);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, ex,
                    fdtd->ex);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, ey,
                    fdtd->ey);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, ez,
                    fdtd->ez);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ,
                    permittivity_inv, fdtd->permittivity_inv);

  float_type _dx = float_cst(1.) / fdtd->dx;
  float_type _dy = float_cst(1.) / fdtd->dy;
  float_type _dz = float_cst(1.) / fdtd->dz;

  for (uintmax_t i = 1; i < fdtd->sizeX; ++i) {
    for (uintmax_t j = 1; j < fdtd->sizeY; ++j) {
      for (uintmax_t k = 1; k < fdtd->sizeZ; ++k) {
        ex[i][j][k] = ex[i][j][k] + ((hz[i][j][k] - hz[i][j - 1][k]) * _dy -
                                     (hy[i][j][k] - hy[i][j][k - 1]) * _dz) *
                                        fdtd->dt * permittivity_inv[i][j][k];
      }
    }
  }
  for (uintmax_t i = 1; i < fdtd->sizeX; ++i) {
    for (uintmax_t j = 1; j < fdtd->sizeY; ++j) {
      for (uintmax_t k = 1; k < fdtd->sizeZ; ++k) {
        ey[i][j][k] = ey[i][j][k] + ((hx[i][j][k] - hx[i][j][k - 1]) * _dz -
                                     (hz[i][j][k] - hz[i - 1][j][k]) * _dx) *
                                        fdtd->dt * permittivity_inv[i][j][k];
      }
    }
  }
  for (uintmax_t i = 1; i < fdtd->sizeX; ++i) {
    for (uintmax_t j = 1; j < fdtd->sizeY; ++j) {
      for (uintmax_t k = 1; k < fdtd->sizeZ; ++k) {
        ez[i][j][k] = ez[i][j][k] + ((hy[i][j][k] - hy[i - 1][j][k]) * _dx -
                                     (hx[i][j][k] - hx[i][j - 1][k]) * _dy) *
                                        fdtd->dt * permittivity_inv[i][j][k];
      }
    }
  }
}

static void update_electric_cpml(struct fdtd3D *fdtd) {
  // Ex
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->cpml_thickness, fdtd->sizeZ,
                    psi_ex_left, fdtd->psi_ex_y[0]);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->cpml_thickness, fdtd->sizeZ,
                    psi_ex_right, fdtd->psi_ex_y[1]);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->cpml_thickness,
                    psi_ex_front, fdtd->psi_ex_z[0]);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->cpml_thickness,
                    psi_ex_back, fdtd->psi_ex_z[1]);
  // Ey
  VLA_3D_definition(float_type, fdtd->cpml_thickness, fdtd->sizeY, fdtd->sizeZ,
                    psi_ey_bottom, fdtd->psi_ey_x[0]);
  VLA_3D_definition(float_type, fdtd->cpml_thickness, fdtd->sizeY, fdtd->sizeZ,
                    psi_ey_top, fdtd->psi_ey_x[1]);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->cpml_thickness,
                    psi_ey_front, fdtd->psi_ey_z[0]);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->cpml_thickness,
                    psi_ey_back, fdtd->psi_ey_z[1]);

  // Ez
  VLA_3D_definition(float_type, fdtd->cpml_thickness, fdtd->sizeY, fdtd->sizeZ,
                    psi_ez_bottom, fdtd->psi_ez_x[0]);
  VLA_3D_definition(float_type, fdtd->cpml_thickness, fdtd->sizeY, fdtd->sizeZ,
                    psi_ez_top, fdtd->psi_ez_x[1]);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->cpml_thickness, fdtd->sizeZ,
                    psi_ez_left, fdtd->psi_ez_y[0]);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->cpml_thickness, fdtd->sizeZ,
                    psi_ez_right, fdtd->psi_ez_y[1]);
  // Sim data
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, hx,
                    fdtd->hx);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, hy,
                    fdtd->hy);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, hz,
                    fdtd->hz);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, ex,
                    fdtd->ex);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, ey,
                    fdtd->ey);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, ez,
                    fdtd->ez);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ,
                    permittivity_inv, fdtd->permittivity_inv);
  // Useful constants
  float_type _dx = float_cst(1.) / fdtd->dx;
  float_type _dy = float_cst(1.) / fdtd->dy;
  float_type _dz = float_cst(1.) / fdtd->dz;

  if (fdtd->border_condition[border_left] & border_cpml) { // ex_y & ez_y
    for (uintmax_t i = 0; i < fdtd->sizeX; ++i) {
      for (uintmax_t j = 0; j < fdtd->cpml_thickness; ++j) {
        for (uintmax_t k = 0; k < fdtd->sizeZ; ++k) {
          psi_ex_left[i][j][k] =
              fdtd->by[j] * psi_ex_left[i][j][k] +
              fdtd->cy[j] * (hz[i][1 + j][k] - hz[i][j][k]) * _dy;
          ex[i][1 + j][k] =
              ex[i][1 + j][k] +
              fdtd->dt * permittivity_inv[i][1 + j][k] * psi_ex_left[i][j][k];
          psi_ez_left[i][j][k] =
              fdtd->by[j] * psi_ez_left[i][j][k] +
              fdtd->cy[j] * (hx[i][1 + j][k] - hx[i][j][k]) * _dy;
          ez[i][1 + j][k] =
              ez[i][1 + j][k] -
              fdtd->dt * permittivity_inv[i][1 + j][k] * psi_ez_left[i][j][k];
        }
      }
    }
  }
  if (fdtd->border_condition[border_right] & border_cpml) { // ex_y & ez_y
    for (uintmax_t i = 0; i < fdtd->sizeX; ++i) {
      for (uintmax_t j = 0; j < fdtd->cpml_thickness; ++j) {
        for (uintmax_t k = 0; k < fdtd->sizeZ; ++k) {
          psi_ex_right[i][j][k] = fdtd->by[j] * psi_ex_right[i][j][k] +
                                  fdtd->cy[j] *
                                      (hz[i][fdtd->sizeY - 1 - j][k] -
                                       hz[i][fdtd->sizeY - 2 - j][k]) *
                                      _dy;
          ex[i][fdtd->sizeY - 1 - j][k] =
              ex[i][fdtd->sizeY - 1 - j][k] +
              fdtd->dt * permittivity_inv[i][fdtd->sizeY - 1 - j][k] *
                  psi_ex_right[i][j][k];
          psi_ez_right[i][j][k] = fdtd->by[j] * psi_ez_right[i][j][k] +
                                  fdtd->cy[j] *
                                      (hx[i][fdtd->sizeY - 1 - j][k] -
                                       hx[i][fdtd->sizeY - 2 - j][k]) *
                                      _dy;
          ez[i][fdtd->sizeY - 1 - j][k] =
              ez[i][fdtd->sizeY - 1 - j][k] -
              fdtd->dt * permittivity_inv[i][fdtd->sizeY - 1 - j][k] *
                  psi_ez_right[i][j][k];
        }
      }
    }
  }
  if (fdtd->border_condition[border_front] & border_cpml) { // ex_z & ey_z
    for (uintmax_t i = 0; i < fdtd->sizeX; ++i) {
      for (uintmax_t j = 0; j < fdtd->sizeY; ++j) {
        for (uintmax_t k = 0; k < fdtd->cpml_thickness; ++k) {
          psi_ex_front[i][j][k] =
              fdtd->bz[k] * psi_ex_front[i][j][k] +
              fdtd->cz[k] * (hy[i][j][1 + k] - hy[i][j][k]) * _dz;
          ex[i][j][1 + k] =
              ex[i][j][1 + k] -
              fdtd->dt * permittivity_inv[i][j][1 + k] * psi_ex_front[i][j][k];
          psi_ey_front[i][j][k] =
              fdtd->bz[k] * psi_ey_front[i][j][k] +
              fdtd->cz[k] * (hx[i][j][1 + k] - hx[i][j][k]) * _dz;
          ey[i][j][1 + k] =
              ey[i][j][1 + k] +
              fdtd->dt * permittivity_inv[i][j][1 + k] * psi_ey_front[i][j][k];
        }
      }
    }
  }
  if (fdtd->border_condition[border_back] & border_cpml) { // ex_z & ey_z
    for (uintmax_t i = 0; i < fdtd->sizeX; ++i) {
      for (uintmax_t j = 0; j < fdtd->sizeY; ++j) {
        for (uintmax_t k = 0; k < fdtd->cpml_thickness; ++k) {
          psi_ex_back[i][j][k] = fdtd->bz[k] * psi_ex_back[i][j][k] +
                                 fdtd->cz[k] *
                                     (hy[i][j][fdtd->sizeZ - 1 - k] -
                                      hy[i][j][fdtd->sizeZ - 2 - k]) *
                                     _dz;
          ex[i][j][fdtd->sizeZ - 1 - k] =
              ex[i][j][fdtd->sizeZ - 1 - k] -
              fdtd->dt * permittivity_inv[i][j][fdtd->sizeZ - 1 - k] *
                  psi_ex_back[i][j][k];
          psi_ey_back[i][j][k] = fdtd->bz[k] * psi_ey_back[i][j][k] +
                                 fdtd->cz[k] *
                                     (hx[i][j][fdtd->sizeZ - 1 - k] -
                                      hx[i][j][fdtd->sizeZ - 2 - k]) *
                                     _dz;
          ey[i][j][fdtd->sizeZ - 1 - k] =
              ey[i][j][fdtd->sizeZ - 1 - k] +
              fdtd->dt * permittivity_inv[i][j][fdtd->sizeZ - 1 - k] *
                  psi_ey_back[i][j][k];
        }
      }
    }
  }
  if (fdtd->border_condition[border_bottom] & border_cpml) { // ey_x & ez_x
    for (uintmax_t i = 0; i < fdtd->cpml_thickness; ++i) {
      for (uintmax_t j = 0; j < fdtd->sizeY; ++j) {
        for (uintmax_t k = 0; k < fdtd->sizeZ; ++k) {
          psi_ey_bottom[i][j][k] =
              fdtd->bx[i] * psi_ey_bottom[i][j][k] +
              fdtd->cx[i] * (hz[1 + i][j][k] - hz[i][j][k]) * _dx;
          ey[1 + i][j][k] =
              ey[1 + i][j][k] -
              fdtd->dt * permittivity_inv[1 + i][j][k] * psi_ey_bottom[i][j][k];
          psi_ez_bottom[i][j][k] =
              fdtd->bx[i] * psi_ez_bottom[i][j][k] +
              fdtd->cx[i] * (hy[1 + i][j][k] - hy[i][j][k]) * _dx;
          ez[1 + i][j][k] =
              ez[1 + i][j][k] +
              fdtd->dt * permittivity_inv[1 + i][j][k] * psi_ez_bottom[i][j][k];
        }
      }
    }
  }
  if (fdtd->border_condition[border_top] & border_cpml) { // ey_x & ez_x
    for (uintmax_t i = 0; i < fdtd->cpml_thickness; ++i) {
      for (uintmax_t j = 0; j < fdtd->sizeY; ++j) {
        for (uintmax_t k = 0; k < fdtd->sizeZ; ++k) {
          psi_ey_top[i][j][k] = fdtd->bx[i] * psi_ey_top[i][j][k] +
                                fdtd->cx[i] *
                                    (hz[fdtd->sizeX - 1 - i][j][k] -
                                     hz[fdtd->sizeX - 2 - i][j][k]) *
                                    _dx;
          ey[fdtd->sizeX - 1 - i][j][k] =
              ey[fdtd->sizeX - 1 - i][j][k] -
              fdtd->dt * permittivity_inv[fdtd->sizeX - 1 - i][j][k] *
                  psi_ey_top[i][j][k];
          psi_ez_top[i][j][k] = fdtd->bx[i] * psi_ez_top[i][j][k] +
                                fdtd->cx[i] *
                                    (hy[fdtd->sizeX - 1 - i][j][k] -
                                     hy[fdtd->sizeX - 2 - i][j][k]) *
                                    _dx;
          ez[fdtd->sizeX - 1 - i][j][k] =
              ez[fdtd->sizeX - 1 - i][j][k] +
              fdtd->dt * permittivity_inv[fdtd->sizeX - 1 - i][j][k] *
                  psi_ez_top[i][j][k];
        }
      }
    }
  }
}

static void border_condition_electric(struct fdtd3D *fdtd) {
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, ex,
                    fdtd->ex);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, ey,
                    fdtd->ey);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, ez,
                    fdtd->ez);

  for (enum border_position3D i = border_front; i < num_borders_3D; ++i) {
    if (fdtd->border_condition[i] == border_perfect_electric_conductor) {
      switch (i) {
      case border_front:
        for (uintmax_t j = 0; j < fdtd->sizeX; ++j)
          for (uintmax_t k = 0; k < fdtd->sizeY; ++k) {
            ex[j][k][0] = float_cst(0.);
            ey[j][k][0] = float_cst(0.);
            ez[j][k][0] = float_cst(0.);
          }
        break;
      case border_back:
        for (uintmax_t j = 0; j < fdtd->sizeX; ++j)
          for (uintmax_t k = 0; k < fdtd->sizeY; ++k) {
            ex[j][k][fdtd->sizeZ - 1] = float_cst(0.);
            ey[j][k][fdtd->sizeZ - 1] = float_cst(0.);
            ez[j][k][fdtd->sizeZ - 1] = float_cst(0.);
          }
        break;
      case border_top:
        for (uintmax_t j = 0; j < fdtd->sizeX; ++j)
          for (uintmax_t k = 0; k < fdtd->sizeZ; ++k) {
            ex[fdtd->sizeX - 1][j][k] = float_cst(0.);
            ey[fdtd->sizeX - 1][j][k] = float_cst(0.);
            ez[fdtd->sizeX - 1][j][k] = float_cst(0.);
          }
        break;
      case border_bottom:
        for (uintmax_t j = 0; j < fdtd->sizeX; ++j)
          for (uintmax_t k = 0; k < fdtd->sizeZ; ++k) {
            ex[0][j][k] = float_cst(0.);
            ey[0][j][k] = float_cst(0.);
            ez[0][j][k] = float_cst(0.);
          }
        break;
      case border_right:
        for (uintmax_t j = 0; j < fdtd->sizeX; ++j)
          for (uintmax_t k = 0; k < fdtd->sizeZ; ++k) {
            ex[j][fdtd->sizeY - 1][k] = float_cst(0.);
            ey[j][fdtd->sizeY - 1][k] = float_cst(0.);
            ez[j][fdtd->sizeY - 1][k] = float_cst(0.);
          }
        break;
      case border_left:
        for (uintmax_t j = 0; j < fdtd->sizeX; ++j)
          for (uintmax_t k = 0; k < fdtd->sizeZ; ++k) {
            ex[j][0][k] = float_cst(0.);
            ey[j][0][k] = float_cst(0.);
            ez[j][0][k] = float_cst(0.);
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

static void update_magnetic_field(struct fdtd3D *fdtd) {
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, hx,
                    fdtd->hx);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, hy,
                    fdtd->hy);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, hz,
                    fdtd->hz);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, ex,
                    fdtd->ex);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, ey,
                    fdtd->ey);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, ez,
                    fdtd->ez);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ,
                    permeability_inv, fdtd->permeability_inv);

  float_type _dy = float_cst(1.) / fdtd->dy;
  float_type _dx = float_cst(1.) / fdtd->dx;
  float_type _dz = float_cst(1.) / fdtd->dz;

  for (uintmax_t i = 0; i < fdtd->sizeX - 1; ++i) {
    for (uintmax_t j = 0; j < fdtd->sizeY - 1; ++j) {
      for (uintmax_t k = 0; k < fdtd->sizeZ - 1; ++k) {
        hx[i][j][k] = hx[i][j][k] + ((ey[i][j][k + 1] - ey[i][j][k]) * _dz -
                                     (ez[i][j + 1][k] - ez[i][j][k]) * _dy) *
                                        fdtd->dt * permeability_inv[i][j][k];
      }
    }
  }
  for (uintmax_t i = 0; i < fdtd->sizeX - 1; ++i) {
    for (uintmax_t j = 0; j < fdtd->sizeY - 1; ++j) {
      for (uintmax_t k = 0; k < fdtd->sizeZ - 1; ++k) {
        hy[i][j][k] = hy[i][j][k] + ((ez[i + 1][j][k] - ez[i][j][k]) * _dx -
                                     (ex[i][j][k + 1] - ex[i][j][k]) * _dz) *
                                        fdtd->dt * permeability_inv[i][j][k];
      }
    }
  }
  for (uintmax_t i = 0; i < fdtd->sizeX - 1; ++i) {
    for (uintmax_t j = 0; j < fdtd->sizeY - 1; ++j) {
      for (uintmax_t k = 0; k < fdtd->sizeZ - 1; ++k) {
        hz[i][j][k] = hz[i][j][k] + ((ex[i][j + 1][k] - ex[i][j][k]) * _dy -
                                     (ey[i + 1][j][k] - ey[i][j][k]) * _dx) *
                                        fdtd->dt * permeability_inv[i][j][k];
      }
    }
  }
}

static void update_magnetic_cpml(struct fdtd3D *fdtd) {
  // Hx
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->cpml_thickness, fdtd->sizeZ,
                    psi_hx_left, fdtd->psi_hx_y[0]);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->cpml_thickness, fdtd->sizeZ,
                    psi_hx_right, fdtd->psi_hx_y[1]);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->cpml_thickness,
                    psi_hx_front, fdtd->psi_hx_z[0]);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->cpml_thickness,
                    psi_hx_back, fdtd->psi_hx_z[1]);
  // Hy
  VLA_3D_definition(float_type, fdtd->cpml_thickness, fdtd->sizeY, fdtd->sizeZ,
                    psi_hy_bottom, fdtd->psi_hy_x[0]);
  VLA_3D_definition(float_type, fdtd->cpml_thickness, fdtd->sizeY, fdtd->sizeZ,
                    psi_hy_top, fdtd->psi_hy_x[1]);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->cpml_thickness,
                    psi_hy_front, fdtd->psi_hy_z[0]);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->cpml_thickness,
                    psi_hy_back, fdtd->psi_hy_z[1]);

  // Hz
  VLA_3D_definition(float_type, fdtd->cpml_thickness, fdtd->sizeY, fdtd->sizeZ,
                    psi_hz_bottom, fdtd->psi_hz_x[0]);
  VLA_3D_definition(float_type, fdtd->cpml_thickness, fdtd->sizeY, fdtd->sizeZ,
                    psi_hz_top, fdtd->psi_hz_x[1]);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->cpml_thickness, fdtd->sizeZ,
                    psi_hz_left, fdtd->psi_hz_y[0]);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->cpml_thickness, fdtd->sizeZ,
                    psi_hz_right, fdtd->psi_hz_y[1]);
  // Sim data
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, hx,
                    fdtd->hx);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, hy,
                    fdtd->hy);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, hz,
                    fdtd->hz);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, ex,
                    fdtd->ex);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, ey,
                    fdtd->ey);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, ez,
                    fdtd->ez);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ,
                    permeability_inv, fdtd->permeability_inv);
  // Useful constants
  float_type _dx = float_cst(1.) / fdtd->dx;
  float_type _dy = float_cst(1.) / fdtd->dy;
  float_type _dz = float_cst(1.) / fdtd->dz;

  if (fdtd->border_condition[border_left] & border_cpml) { // hx_y & hz_y
    for (uintmax_t i = 0; i < fdtd->sizeX; ++i) {
      for (uintmax_t j = 0; j < fdtd->cpml_thickness; ++j) {
        for (uintmax_t k = 0; k < fdtd->sizeZ; ++k) {
          psi_hx_left[i][j][k] =
              fdtd->by[j] * psi_hx_left[i][j][k] +
              fdtd->cy[j] * (ez[i][j + 1][k] - ez[i][j][k]) * _dy;
          hx[i][j][k] = hx[i][j][k] - fdtd->dt * permeability_inv[i][j][k] *
                                          psi_hx_left[i][j][k];
          psi_hz_left[i][j][k] =
              fdtd->by[j] * psi_hz_left[i][j][k] +
              fdtd->cy[j] * (ex[i][j + 1][k] - ex[i][j][k]) * _dy;
          hz[i][j][k] = hz[i][j][k] + fdtd->dt * permeability_inv[i][j][k] *
                                          psi_hz_left[i][j][k];
        }
      }
    }
  }
  if (fdtd->border_condition[border_right] & border_cpml) { // hx_y & hz_y
    for (uintmax_t i = 0; i < fdtd->sizeX; ++i) {
      for (uintmax_t j = 0; j < fdtd->cpml_thickness; ++j) {
        for (uintmax_t k = 0; k < fdtd->sizeZ; ++k) {
          psi_hx_right[i][j][k] = fdtd->by[j] * psi_hx_right[i][j][k] +
                                  fdtd->cy[j] *
                                      (ez[i][fdtd->sizeY - 1 - j][k] -
                                       ez[i][fdtd->sizeY - 2 - j][k]) *
                                      _dy;
          hx[i][fdtd->sizeY - 2 - j][k] =
              hx[i][fdtd->sizeY - 2 - j][k] -
              fdtd->dt * permeability_inv[i][fdtd->sizeY - 2 - j][k] *
                  psi_hx_right[i][j][k];
          psi_hz_right[i][j][k] = fdtd->by[j] * psi_hz_right[i][j][k] +
                                  fdtd->cy[j] *
                                      (ex[i][fdtd->sizeY - 1 - j][k] -
                                       ex[i][fdtd->sizeY - 2 - j][k]) *
                                      _dy;
          hz[i][fdtd->sizeY - 2 - j][k] =
              hz[i][fdtd->sizeY - 2 - j][k] +
              fdtd->dt * permeability_inv[i][fdtd->sizeY - 2 - j][k] *
                  psi_hz_right[i][j][k];
        }
      }
    }
  }
  if (fdtd->border_condition[border_front] & border_cpml) { // hx_z & hy_z
    for (uintmax_t i = 0; i < fdtd->sizeX; ++i) {
      for (uintmax_t j = 0; j < fdtd->sizeY; ++j) {
        for (uintmax_t k = 0; k < fdtd->cpml_thickness; ++k) {
          psi_hx_front[i][j][k] =
              fdtd->bz[k] * psi_hx_front[i][j][k] +
              fdtd->cz[k] * (ey[i][j][k + 1] - ey[i][j][k]) * _dz;
          hx[i][j][k] = hx[i][j][k] + fdtd->dt * permeability_inv[i][j][k] *
                                          psi_hx_front[i][j][k];
          psi_hy_front[i][j][k] =
              fdtd->bz[k] * psi_hy_front[i][j][k] +
              fdtd->cz[k] * (ex[i][j][k + 1] - ex[i][j][k]) * _dz;
          hy[i][j][k] = hy[i][j][k] - fdtd->dt * permeability_inv[i][j][k] *
                                          psi_hy_front[i][j][k];
        }
      }
    }
  }
  if (fdtd->border_condition[border_back] & border_cpml) { // hx_z & hy_z
    for (uintmax_t i = 0; i < fdtd->sizeX; ++i) {
      for (uintmax_t j = 0; j < fdtd->sizeY; ++j) {
        for (uintmax_t k = 0; k < fdtd->cpml_thickness; ++k) {
          psi_hx_back[i][j][k] = fdtd->bz[k] * psi_hx_back[i][j][k] +
                                 fdtd->cz[k] *
                                     (ey[i][j][fdtd->sizeZ - 1 - k] -
                                      ey[i][j][fdtd->sizeZ - 2 - k]) *
                                     _dz;
          hx[i][j][fdtd->sizeZ - 2 - k] =
              hx[i][j][fdtd->sizeZ - 2 - k] +
              fdtd->dt * permeability_inv[i][j][fdtd->sizeZ - 2 - k] *
                  psi_hx_back[i][j][k];
          psi_hy_back[i][j][k] = fdtd->bz[k] * psi_hy_back[i][j][k] +
                                 fdtd->cz[k] *
                                     (ex[i][j][fdtd->sizeZ - 1 - k] -
                                      ex[i][j][fdtd->sizeZ - 2 - k]) *
                                     _dz;
          hy[i][j][fdtd->sizeZ - 2 - k] =
              hy[i][j][fdtd->sizeZ - 2 - k] -
              fdtd->dt * permeability_inv[i][j][fdtd->sizeZ - 2 - k] *
                  psi_hy_back[i][j][k];
        }
      }
    }
  }
  if (fdtd->border_condition[border_bottom] & border_cpml) { // hy_x & hz_x
    for (uintmax_t i = 0; i < fdtd->cpml_thickness; ++i) {
      for (uintmax_t j = 0; j < fdtd->sizeY; ++j) {
        for (uintmax_t k = 0; k < fdtd->sizeZ; ++k) {
          psi_hy_bottom[i][j][k] =
              fdtd->bx[i] * psi_hy_bottom[i][j][k] +
              fdtd->cx[i] * (ez[i + 1][j][k] - ez[i][j][k]) * _dx;
          hy[i][j][k] = hy[i][j][k] + fdtd->dt * permeability_inv[i][j][k] *
                                          psi_hy_bottom[i][j][k];
          psi_hz_bottom[i][j][k] =
              fdtd->bx[i] * psi_hz_bottom[i][j][k] +
              fdtd->cx[i] * (ey[i + 1][j][k] - ey[i][j][k]) * _dx;
          hz[i][j][k] = hz[i][j][k] - fdtd->dt * permeability_inv[i][j][k] *
                                          psi_hz_bottom[i][j][k];
        }
      }
    }
  }
  if (fdtd->border_condition[border_top] & border_cpml) { // hy_x & hz_x
    for (uintmax_t i = 0; i < fdtd->cpml_thickness; ++i) {
      for (uintmax_t j = 0; j < fdtd->sizeY; ++j) {
        for (uintmax_t k = 0; k < fdtd->sizeZ; ++k) {
          psi_hy_top[i][j][k] = fdtd->bx[i] * psi_hy_top[i][j][k] +
                                fdtd->cx[i] *
                                    (ez[fdtd->sizeX - 1 - i][j][k] -
                                     ez[fdtd->sizeX - 2 - i][j][k]) *
                                    _dx;
          hy[fdtd->sizeX - 2 - i][j][k] =
              hy[fdtd->sizeX - 2 - i][j][k] +
              fdtd->dt * permeability_inv[fdtd->sizeX - 2 - i][j][k] *
                  psi_hy_top[i][j][k];
          psi_hz_top[i][j][k] = fdtd->bx[i] * psi_hz_top[i][j][k] +
                                fdtd->cx[i] *
                                    (ey[fdtd->sizeX - 1 - i][j][k] -
                                     ey[fdtd->sizeX - 2 - i][j][k]) *
                                    _dx;
          hz[fdtd->sizeX - 2 - i][j][k] =
              hz[fdtd->sizeX - 2 - i][j][k] -
              fdtd->dt * permeability_inv[fdtd->sizeX - 2 - i][j][k] *
                  psi_hz_top[i][j][k];
        }
      }
    }
  }
}

static void border_condition_magnetic(struct fdtd3D *fdtd) {
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, hx,
                    fdtd->hx);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, hy,
                    fdtd->hy);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, hz,
                    fdtd->hz);
  for (enum border_position3D i = border_front; i < num_borders_3D; ++i) {
    if (fdtd->border_condition[i] == border_perfect_electric_conductor) {
      switch (i) {
      case border_front:
        for (uintmax_t j = 0; j < fdtd->sizeX; ++j)
          for (uintmax_t k = 0; k < fdtd->sizeY; ++k) {
            hx[j][k][0] = float_cst(0.);
            hy[j][k][0] = float_cst(0.);
            hz[j][k][0] = float_cst(0.);
          }
        break;
      case border_back:
        for (uintmax_t j = 0; j < fdtd->sizeX; ++j)
          for (uintmax_t k = 0; k < fdtd->sizeY; ++k) {
            hx[j][k][fdtd->sizeZ - 1] = float_cst(0.);
            hy[j][k][fdtd->sizeZ - 1] = float_cst(0.);
            hz[j][k][fdtd->sizeZ - 1] = float_cst(0.);
          }
        break;
      case border_top:
        for (uintmax_t j = 0; j < fdtd->sizeX; ++j)
          for (uintmax_t k = 0; k < fdtd->sizeZ; ++k) {
            hx[fdtd->sizeX - 1][j][k] = float_cst(0.);
            hy[fdtd->sizeX - 1][j][k] = float_cst(0.);
            hz[fdtd->sizeX - 1][j][k] = float_cst(0.);
          }
        break;
      case border_bottom:
        for (uintmax_t j = 0; j < fdtd->sizeX; ++j)
          for (uintmax_t k = 0; k < fdtd->sizeZ; ++k) {
            hx[0][j][k] = float_cst(0.);
            hy[0][j][k] = float_cst(0.);
            hz[0][j][k] = float_cst(0.);
          }
        break;
      case border_right:
        for (uintmax_t j = 0; j < fdtd->sizeX; ++j)
          for (uintmax_t k = 0; k < fdtd->sizeZ; ++k) {
            hx[j][fdtd->sizeY - 1][k] = float_cst(0.);
            hy[j][fdtd->sizeY - 1][k] = float_cst(0.);
            hz[j][fdtd->sizeY - 1][k] = float_cst(0.);
          }
        break;
      case border_left:
        for (uintmax_t j = 0; j < fdtd->sizeX; ++j)
          for (uintmax_t k = 0; k < fdtd->sizeZ; ++k) {
            hx[j][0][k] = float_cst(0.);
            hy[j][0][k] = float_cst(0.);
            hz[j][0][k] = float_cst(0.);
          }
        break;
      default:
        fprintf(stderr, "Error while processing the magnetic border "
                        "conditions: unknown border\n");
        exit(EXIT_FAILURE);
      }
    }
  }
}

static void apply_M_sources(struct fdtd3D *fdtd) {
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, hx,
                    fdtd->hx);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, hy,
                    fdtd->hy);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, hz,
                    fdtd->hz);
  VLA_2D_definition(uintmax_t, fdtd->num_Msources, 3, MsourceLocations,
                    fdtd->MsourceLocations);
  for (uintmax_t i = 0; i < fdtd->num_Msources; ++i) {
    switch (fdtd->Msources[i].type) {
    case source_gaussian_pulse:
      hy[MsourceLocations[i][0]][MsourceLocations[i][1]]
        [MsourceLocations[i][2]] +=
          gaussian_pulse_val(fdtd->time, &fdtd->Msources[i]);
      hx[MsourceLocations[i][0]][MsourceLocations[i][1]]
        [MsourceLocations[i][2]] +=
          gaussian_pulse_val(fdtd->time, &fdtd->Msources[i]);
      hz[MsourceLocations[i][0]][MsourceLocations[i][1]]
        [MsourceLocations[i][2]] +=
          gaussian_pulse_val(fdtd->time, &fdtd->Msources[i]);
      break;
    default:
      break;
    }
  }
}

static void apply_J_sources(struct fdtd3D *fdtd) {
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, ex,
                    fdtd->ex);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, ey,
                    fdtd->ey);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, ez,
                    fdtd->ez);
  VLA_2D_definition(uintmax_t, fdtd->num_Jsources, 3, JsourceLocations,
                    fdtd->JsourceLocations);
  for (uintmax_t i = 0; i < fdtd->num_Jsources; ++i) {
    switch (fdtd->Jsources[i].type) {
    case source_gaussian_pulse:
      ex[JsourceLocations[i][0]][JsourceLocations[i][1]]
        [JsourceLocations[i][1]] +=
          gaussian_pulse_val(fdtd->time, &fdtd->Jsources[i]);
      ey[JsourceLocations[i][0]][JsourceLocations[i][1]]
        [JsourceLocations[i][1]] +=
          gaussian_pulse_val(fdtd->time, &fdtd->Jsources[i]);
      ez[JsourceLocations[i][0]][JsourceLocations[i][1]]
        [JsourceLocations[i][1]] +=
          gaussian_pulse_val(fdtd->time, &fdtd->Jsources[i]);
      break;
    default:
      break;
    }
  }
}

void init_fdtd_3D_medium(struct fdtd3D *fdtd,
                         init_medium_fun_3D permeability_invR,
                         init_medium_fun_3D permittivity_invR, void *user) {
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ,
                    permittivity_inv, fdtd->permittivity_inv);
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ,
                    permeability_inv, fdtd->permeability_inv);
  float_type posX = float_cst(0.);
  for (uintmax_t i = 0; i < fdtd->sizeX; ++i, posX += fdtd->dx) {
    float_type posY = float_cst(0.);
    for (uintmax_t j = 0; j < fdtd->sizeY; ++j, posY += fdtd->dy) {
      float_type posZ = float_cst(0.);
      for (uintmax_t k = 0; k < fdtd->sizeZ; ++k, posZ += fdtd->dz) {
        permeability_inv[i][j][k] =
            float_cst(1.) / (permeability_invR(posX, posY, posZ, user) * mu0);
        permittivity_inv[i][j][k] =
            float_cst(1.) / (permittivity_invR(posX, posY, posZ, user) * eps0);
      }
    }
  }
}

struct fdtd3D init_fdtd_3D(float_type domain_size[3], float_type Sc,
                           float_type smallest_wavelength,
                           enum border_condition borders[num_borders_3D]) {
  return init_fdtd_3D_cpml(domain_size, Sc, smallest_wavelength, borders, 0);
}

struct fdtd3D init_fdtd_3D_cpml(float_type domain_size[3], float_type Sc,
                                float_type smallest_wavelength,
                                enum border_condition borders[num_borders_3D],
                                uintmax_t cpml_thickness) {

  float_type dx = smallest_wavelength / float_cst(20.);
  float_type dy = dx;
  float_type dz = dx;
  float_type dt = dx * Sc / c_light;
  float_type sizeXf = ceil(domain_size[0] / dx);
  float_type sizeYf = ceil(domain_size[1] / dy);
  float_type sizeZf = ceil(domain_size[2] / dy);
  uintmax_t sizeX = (uintmax_t)sizeXf;
  uintmax_t sizeY = (uintmax_t)sizeYf;
  uintmax_t sizeZ = (uintmax_t)sizeZf;
  float_type Sc_max = float_cst(1.) / sqrt(float_cst(3.));
  if (Sc > Sc_max)
    fprintf(stderr,
            "The value of Sc is too high, the simulation may be unstable. "
            "Please use a value lesser or equal to %.5f\n",
            Sc_max);

  struct fdtd3D fdtd = {
      .dx = dx,
      .dy = dy,
      .dz = dz,
      .dt = dt,
      .hx = calloc(1, VLA_3D_size(float_type, sizeX, sizeY, sizeZ)),
      .hy = calloc(1, VLA_3D_size(float_type, sizeX, sizeY, sizeZ)),
      .hz = calloc(1, VLA_3D_size(float_type, sizeX, sizeY, sizeZ)),
      .ex = calloc(1, VLA_3D_size(float_type, sizeX, sizeY, sizeZ)),
      .ey = calloc(1, VLA_3D_size(float_type, sizeX, sizeY, sizeZ)),
      .ez = calloc(1, VLA_3D_size(float_type, sizeX, sizeY, sizeZ)),
      .permittivity_inv =
          calloc(1, VLA_3D_size(float_type, sizeX, sizeY, sizeZ)),
      .permeability_inv =
          calloc(1, VLA_3D_size(float_type, sizeX, sizeY, sizeZ)),
      .psi_hx_y = {NULL, NULL},
      .psi_hx_z = {NULL, NULL},
      .psi_hy_x = {NULL, NULL},
      .psi_hy_z = {NULL, NULL},
      .psi_hz_x = {NULL, NULL},
      .psi_hz_y = {NULL, NULL},
      .psi_ex_y = {NULL, NULL},
      .psi_ex_z = {NULL, NULL},
      .psi_ey_x = {NULL, NULL},
      .psi_ey_z = {NULL, NULL},
      .psi_ez_x = {NULL, NULL},
      .psi_ez_y = {NULL, NULL},
      .bx = NULL,
      .by = NULL,
      .bz = NULL,
      .cx = NULL,
      .cy = NULL,
      .cz = NULL,
      .cpml_thickness = cpml_thickness,
      .border_condition =
          {
              [border_front] = borders[border_front],
              [border_back] = borders[border_back],
              [border_top] = borders[border_top],
              [border_bottom] = borders[border_bottom],
              [border_right] = borders[border_right],
              [border_left] = borders[border_left],
          },
      .domain_size = {domain_size[0], domain_size[1], domain_size[2]},
      .sizeX = sizeX,
      .sizeY = sizeY,
      .sizeZ = sizeZ,
      .Sc = Sc,
      .num_Jsources = 0,
      .Jsources = NULL,
      .JsourceLocations = NULL,
      .num_Msources = 0,
      .Msources = NULL,
      .MsourceLocations = NULL,
      .time = float_cst(0.),
  };

  if (cpml_thickness > 0 && fdtd.border_condition[border_front] & border_cpml) {
    fdtd.psi_hx_z[0] =
        calloc(1, VLA_3D_size(float_type, sizeX, sizeY, cpml_thickness));
    fdtd.psi_hy_z[0] =
        calloc(1, VLA_3D_size(float_type, sizeX, sizeY, cpml_thickness));
    fdtd.psi_ex_z[0] =
        calloc(1, VLA_3D_size(float_type, sizeX, sizeY, cpml_thickness));
    fdtd.psi_ey_z[0] =
        calloc(1, VLA_3D_size(float_type, sizeX, sizeY, cpml_thickness));
  }
  if (cpml_thickness > 0 && fdtd.border_condition[border_back] & border_cpml) {
    fdtd.psi_hx_z[1] =
        calloc(1, VLA_3D_size(float_type, sizeX, sizeY, cpml_thickness));
    fdtd.psi_hy_z[1] =
        calloc(1, VLA_3D_size(float_type, sizeX, sizeY, cpml_thickness));
    fdtd.psi_ex_z[1] =
        calloc(1, VLA_3D_size(float_type, sizeX, sizeY, cpml_thickness));
    fdtd.psi_ey_z[1] =
        calloc(1, VLA_3D_size(float_type, sizeX, sizeY, cpml_thickness));
  }
  if (cpml_thickness > 0 &&
      fdtd.border_condition[border_bottom] & border_cpml) {
    fdtd.psi_hy_x[0] =
        calloc(1, VLA_3D_size(float_type, cpml_thickness, sizeY, sizeZ));
    fdtd.psi_hz_x[0] =
        calloc(1, VLA_3D_size(float_type, cpml_thickness, sizeY, sizeZ));
    fdtd.psi_ey_x[0] =
        calloc(1, VLA_3D_size(float_type, cpml_thickness, sizeY, sizeZ));
    fdtd.psi_ez_x[0] =
        calloc(1, VLA_3D_size(float_type, cpml_thickness, sizeY, sizeZ));
  }
  if (cpml_thickness > 0 && fdtd.border_condition[border_top] & border_cpml) {
    fdtd.psi_hy_x[1] =
        calloc(1, VLA_3D_size(float_type, cpml_thickness, sizeY, sizeZ));
    fdtd.psi_hz_x[1] =
        calloc(1, VLA_3D_size(float_type, cpml_thickness, sizeY, sizeZ));
    fdtd.psi_ey_x[1] =
        calloc(1, VLA_3D_size(float_type, cpml_thickness, sizeY, sizeZ));
    fdtd.psi_ez_x[1] =
        calloc(1, VLA_3D_size(float_type, cpml_thickness, sizeY, sizeZ));
  }
  if (cpml_thickness > 0 && fdtd.border_condition[border_left] & border_cpml) {
    fdtd.psi_hx_y[0] =
        calloc(1, VLA_3D_size(float_type, sizeX, cpml_thickness, sizeZ));
    fdtd.psi_hz_y[0] =
        calloc(1, VLA_3D_size(float_type, sizeX, cpml_thickness, sizeZ));
    fdtd.psi_ex_y[0] =
        calloc(1, VLA_3D_size(float_type, sizeX, cpml_thickness, sizeZ));
    fdtd.psi_ez_y[0] =
        calloc(1, VLA_3D_size(float_type, sizeX, cpml_thickness, sizeZ));
  }
  if (cpml_thickness > 0 && fdtd.border_condition[border_right] & border_cpml) {
    fdtd.psi_hx_y[1] =
        calloc(1, VLA_3D_size(float_type, sizeX, cpml_thickness, sizeZ));
    fdtd.psi_hz_y[1] =
        calloc(1, VLA_3D_size(float_type, sizeX, cpml_thickness, sizeZ));
    fdtd.psi_ex_y[1] =
        calloc(1, VLA_3D_size(float_type, sizeX, cpml_thickness, sizeZ));
    fdtd.psi_ez_y[1] =
        calloc(1, VLA_3D_size(float_type, sizeX, cpml_thickness, sizeZ));
  }
  fdtd.bx = calloc(cpml_thickness, sizeof(*fdtd.bx));
  fdtd.by = fdtd.bx;
  fdtd.bz = fdtd.bx;
  fdtd.cx = calloc(cpml_thickness, sizeof(*fdtd.bx));
  fdtd.cy = fdtd.cx;
  fdtd.cz = fdtd.cx;
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

  fprintf(stderr, "Dt %e Dx %e Dy %e Dz %e (%.0fx%.0fx%.0f)\n", dt, dx, dy, dz,
          sizeXf, sizeYf, sizeZf);

  return fdtd;
}

void run_3D_fdtd(struct fdtd3D *fdtd, float_type end_time) {
  fprintf(stderr, "It will take %.0f iterations\n",
          (end_time - fdtd->time) / fdtd->dt);
  /*exit(1);*/
  for (; fdtd->time < end_time; fdtd->time += fdtd->dt) {

    update_magnetic_field(fdtd);
    apply_M_sources(fdtd);
    update_magnetic_cpml(fdtd);
    border_condition_magnetic(fdtd);

    update_electric_field(fdtd);
    apply_J_sources(fdtd);
    update_electric_cpml(fdtd);
    border_condition_electric(fdtd);

    /*fprintf(stderr, "Time %e/%e\n", fdtd->time, end_time);*/
  }
}

void dump_3D_fdtd(const struct fdtd3D *fdtd, const char *fileName,
                  enum dumpable_data what_to_dump) {
  FILE *out = fopen(fileName, "w");
  VLA_3D_definition(float_type, fdtd->sizeX, fdtd->sizeY, fdtd->sizeZ, data,
                    fdtd->ez);
  switch (what_to_dump) {
  case dump_ex:
    break;
  case dump_ey:
    data = fdtd->ey;
    break;
  case dump_ez:
    data = fdtd->ez;
    break;
  case dump_hx:
    data = fdtd->hx;
    break;
  case dump_hy:
    data = fdtd->hy;
    break;
  case dump_hz:
    data = fdtd->hz;
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
      for (uintmax_t k = 1; k < fdtd->sizeZ; ++k) {
        fprintf(out, "%e %e %e %e\n", (float_type)i * fdtd->dx,
                (float_type)j * fdtd->dy, (float_type)k * fdtd->dt,
                data[i][j][k]);
      }
    }
  }
  fclose(out);
}

void free_3D_fdtd(struct fdtd3D *fdtd) {
  free(fdtd->ex);
  free(fdtd->ey);
  free(fdtd->ez);
  free(fdtd->hx);
  free(fdtd->hy);
  free(fdtd->hz);
  free(fdtd->permittivity_inv);
  free(fdtd->permeability_inv);
  free(fdtd->JsourceLocations);
  free(fdtd->Jsources);
  free(fdtd->MsourceLocations);
  free(fdtd->Msources);
  free(fdtd->psi_hx_z[0]);
  free(fdtd->psi_hx_z[1]);
  free(fdtd->psi_hx_y[0]);
  free(fdtd->psi_hx_y[1]);
  free(fdtd->psi_hy_z[0]);
  free(fdtd->psi_hy_z[1]);
  free(fdtd->psi_hy_x[0]);
  free(fdtd->psi_hy_x[1]);
  free(fdtd->psi_hz_x[0]);
  free(fdtd->psi_hz_x[1]);
  free(fdtd->psi_hz_y[0]);
  free(fdtd->psi_hz_y[1]);
  free(fdtd->psi_ex_z[0]);
  free(fdtd->psi_ex_z[1]);
  free(fdtd->psi_ex_y[0]);
  free(fdtd->psi_ex_y[1]);
  free(fdtd->psi_ey_z[0]);
  free(fdtd->psi_ey_z[1]);
  free(fdtd->psi_ey_x[0]);
  free(fdtd->psi_ey_x[1]);
  free(fdtd->psi_ez_x[0]);
  free(fdtd->psi_ez_x[1]);
  free(fdtd->psi_ez_y[0]);
  free(fdtd->psi_ez_y[1]);
  free(fdtd->bx);
  free(fdtd->cx);
}

void add_source_fdtd_3D(enum source_type sType, struct fdtd3D *fdtd,
                        struct fdtd_source src, float_type positionX,
                        float_type positionY, float_type positionZ) {
  float_type posX = ceil(positionX / fdtd->dx);
  float_type posY = ceil(positionY / fdtd->dy);
  float_type posZ = ceil(positionZ / fdtd->dz);
  switch (sType) {
  case source_electric: {
    fdtd->num_Jsources++;
    fdtd->JsourceLocations = realloc(
        fdtd->JsourceLocations, VLA_2D_size(uintmax_t, fdtd->num_Jsources, 3));
    VLA_2D_definition(uintmax_t, fdtd->num_Jsources, 3, JsourceLocations,
                      fdtd->JsourceLocations);
    JsourceLocations[fdtd->num_Jsources - 1][0] = (uintmax_t)posX;
    JsourceLocations[fdtd->num_Jsources - 1][1] = (uintmax_t)posY;
    JsourceLocations[fdtd->num_Jsources - 1][2] = (uintmax_t)posZ;
    fdtd->Jsources =
        realloc(fdtd->Jsources, fdtd->num_Jsources * sizeof(*fdtd->Jsources));
    fdtd->Jsources[fdtd->num_Jsources - 1] = src;
    break;
  }
  case source_magnetic:
    fdtd->num_Msources++;
    fdtd->MsourceLocations = realloc(
        fdtd->MsourceLocations, VLA_2D_size(uintmax_t, fdtd->num_Msources, 3));
    VLA_2D_definition(uintmax_t, fdtd->num_Msources, 3, MsourceLocations,
                      fdtd->MsourceLocations);
    MsourceLocations[fdtd->num_Msources - 1][0] = (uintmax_t)posX;
    MsourceLocations[fdtd->num_Msources - 1][1] = (uintmax_t)posY;
    MsourceLocations[fdtd->num_Msources - 1][2] = (uintmax_t)posZ;
    fdtd->Msources =
        realloc(fdtd->Msources, fdtd->num_Msources * sizeof(*fdtd->Msources));
    fdtd->Msources[fdtd->num_Msources - 1] = src;
    break;
  }
}
