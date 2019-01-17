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

#include <initialize.h>
#include <stdio.h>
#include <stdlib.h>

#include "fdtd.h"
#include "fdtd_common.h"

struct two_medium_data {
  float_type permeability1, permeability2;
  float_type permittivity1, permittivity2;
  float_type switch_location;
};

struct middle_object_2D {
  float_type permeability_medium, permeability_object;
  float_type permittivity_medium, permittivity_object;
  float_type object_center[2];
  float_type object_dimensions[2];
};

struct middle_object_3D {
  float_type permeability_medium, permeability_object;
  float_type permittivity_medium, permittivity_object;
  float_type object_center[3];
  float_type object_dimensions[3];
};

static float_type init_permeability_two_parts_1D(float_type pos, void *user) {
  struct two_medium_data *tpd = (struct two_medium_data *)user;
  if (pos < tpd->switch_location) {
    return tpd->permeability1;
  } else {
    return tpd->permeability2;
  }
}

static float_type init_permittivity_two_parts_1D(float_type pos, void *user) {
  struct two_medium_data *tpd = (struct two_medium_data *)user;
  if (pos < tpd->switch_location) {
    return tpd->permittivity1;
  } else {
    return tpd->permittivity2;
  }
}

static float_type init_permeability_object_2D(float_type posX, float_type posY,
                                              void *user) {
  struct middle_object_2D *mo = (struct middle_object_2D *)user;
  if (posX < mo->object_center[0] - mo->object_dimensions[0] / float_cst(2.) ||
      posX > mo->object_center[0] + mo->object_dimensions[0] / float_cst(2.)) {
    return mo->permeability_medium;
  } else {
    if (posY <
            mo->object_center[1] - mo->object_dimensions[1] / float_cst(2.) ||
        posY >
            mo->object_center[1] + mo->object_dimensions[1] / float_cst(2.)) {
      return mo->permeability_medium;
    } else {
      return mo->permeability_object;
    }
  }
}

static float_type init_permittivity_object_2D(float_type posX, float_type posY,
                                              void *user) {
  struct middle_object_2D *mo = (struct middle_object_2D *)user;
  if (posX < mo->object_center[0] - mo->object_dimensions[0] / float_cst(2.) ||
      posX > mo->object_center[0] + mo->object_dimensions[0] / float_cst(2.)) {
    return mo->permittivity_medium;
  } else {
    if (posY <
            mo->object_center[1] - mo->object_dimensions[1] / float_cst(2.) ||
        posY >
            mo->object_center[1] + mo->object_dimensions[1] / float_cst(2.)) {
      return mo->permittivity_medium;
    } else {
      return mo->permittivity_object;
    }
  }
}

static float_type init_permeability_object_3D(float_type posX, float_type posY,
                                              float_type posZ, void *user) {
  struct middle_object_3D *mo = (struct middle_object_3D *)user;
  if (posX < mo->object_center[0] - mo->object_dimensions[0] / float_cst(2.) ||
      posX > mo->object_center[0] + mo->object_dimensions[0] / float_cst(2.)) {
    return mo->permeability_medium;
  } else {
    if (posY <
            mo->object_center[1] - mo->object_dimensions[1] / float_cst(2.) ||
        posY >
            mo->object_center[1] + mo->object_dimensions[1] / float_cst(2.)) {
      return mo->permeability_medium;
    } else {
      if (posZ <
              mo->object_center[2] - mo->object_dimensions[2] / float_cst(2.) ||
          posZ >
              mo->object_center[2] + mo->object_dimensions[2] / float_cst(2.)) {
        return mo->permeability_medium;
      } else {
        return mo->permeability_object;
      }
    }
  }
}

static float_type init_permittivity_object_3D(float_type posX, float_type posY,
                                              float_type posZ, void *user) {
  struct middle_object_3D *mo = (struct middle_object_3D *)user;
  if (posX < mo->object_center[0] - mo->object_dimensions[0] / float_cst(2.) ||
      posX > mo->object_center[0] + mo->object_dimensions[0] / float_cst(2.)) {
    return mo->permittivity_medium;
  } else {
    if (posY <
            mo->object_center[1] - mo->object_dimensions[1] / float_cst(2.) ||
        posY >
            mo->object_center[1] + mo->object_dimensions[1] / float_cst(2.)) {
      return mo->permittivity_medium;
    } else {
      if (posZ <
              mo->object_center[2] - mo->object_dimensions[2] / float_cst(2.) ||
          posZ >
              mo->object_center[2] + mo->object_dimensions[2] / float_cst(2.)) {
        return mo->permittivity_medium;
      } else {
        return mo->permittivity_object;
      }
    }
  }
}

static struct fdtd initializeFdtd1D(unsigned setupID, float_type domain_size,
                                    float_type Sc,
                                    float_type smallest_wavelength) {
  enum setup1D sID = (enum setup1D)setupID;
  switch (sID) {
  case half_air_half_water_1D: {
    enum border_condition bc[2] = {border_perfect_electric_conductor,
                                   border_perfect_magnetic_conductor};
    struct fdtd1D fdtd = init_fdtd_1D(domain_size, Sc, smallest_wavelength, bc);
    struct two_medium_data tmd = {.permittivity1 = float_cst(1.00058986),
                                  .permittivity2 = float_cst(78.4),
                                  .permeability1 = float_cst(1.00000037),
                                  .permeability2 = float_cst(0.999992),
                                  .switch_location =
                                      domain_size / float_cst(2.)};
    init_fdtd_1D_medium(&fdtd, init_permeability_two_parts_1D,
                        init_permittivity_two_parts_1D, &tmd);
    struct fdtd_source src = gaussian_source(
        float_cst(25.) * fdtd.dt, float_cst(3.) * fdtd.dt, float_cst(1.e-2));
    add_source_fdtd_1D(source_magnetic, &fdtd, src, float_cst(0.));
    struct fdtd retval = {.oneDim = fdtd, .type = fdtd_one_dim};
    return retval;
  }
  default:
    fprintf(stderr, "The specified 1D setup ID is does not exist\n");
    exit(EXIT_FAILURE);
  }
}

static struct fdtd initializeFdtd2D(unsigned setupID, float_type *domain_size,
                                    float_type Sc,
                                    float_type smallest_wavelength,
                                    uintmax_t cpml_thickness) {
  enum setup2D sID = (enum setup2D)setupID;
  switch (sID) {
  case object_high_permitivity_in_air_west_gaussian_pulse_centered_2D: {
    enum border_condition bc[num_borders_2D] = {
        [border_south] = border_perfect_electric_conductor | border_cpml,
        [border_north] = border_perfect_electric_conductor | border_cpml,
        [border_east] = border_perfect_electric_conductor,
        [border_west] = border_perfect_electric_conductor};
    struct fdtd2D fdtd = init_fdtd_2D_cpml(domain_size, Sc, smallest_wavelength,
                                           bc, cpml_thickness);
    float_type smallest_dim_size = fdtd.domain_size[0] > fdtd.domain_size[1]
                                       ? fdtd.domain_size[1]
                                       : fdtd.domain_size[0];
    struct middle_object_2D mo = {
        .permittivity_medium = float_cst(1.00058986),
        .permittivity_object = float_cst(1e9),
        .permeability_medium = float_cst(1.00000037),
        .permeability_object = float_cst(1.),
        .object_center = {domain_size[0] / float_cst(2.),
                          domain_size[1] / float_cst(2.)},
        .object_dimensions = {smallest_dim_size / float_cst(2.),
                              smallest_dim_size / float_cst(2.)}};
    init_fdtd_2D_medium(&fdtd, init_permeability_object_2D,
                        init_permittivity_object_2D, &mo);

    struct fdtd_source src = gaussian_source(
        float_cst(25.) * fdtd.dt, float_cst(3.) * fdtd.dt, float_cst(100));
    for (uintmax_t j = fdtd.cpml_thickness;
         j < fdtd.sizeX - fdtd.cpml_thickness; ++j) {
      add_source_fdtd_2D(source_electric, &fdtd, src, (float_type)j * fdtd.dx,
                         fdtd.dy);
    }

    struct fdtd retval = {.twoDims = fdtd, .type = fdtd_two_dims};
    return retval;
  }
  case free_space_gaussian_exitation_centered_absorbing_border_2D: {
    enum border_condition bc[num_borders_2D] = {
        [border_south] = border_perfect_electric_conductor | border_cpml,
        [border_north] = border_perfect_electric_conductor | border_cpml,
        [border_east] = border_perfect_electric_conductor | border_cpml,
        [border_west] = border_perfect_electric_conductor | border_cpml};
    struct fdtd2D fdtd = init_fdtd_2D_cpml(domain_size, Sc, smallest_wavelength,
                                           bc, cpml_thickness);

    struct middle_object_2D mo = {
        .permittivity_medium = float_cst(1.),
        .permittivity_object = float_cst(1.),
        .permeability_medium = float_cst(1.),
        .permeability_object = float_cst(1.),
        .object_center = {float_cst(-1.), float_cst(-1.)},
        .object_dimensions = {float_cst(0.), float_cst(0.)}};
    init_fdtd_2D_medium(&fdtd, init_permeability_object_2D,
                        init_permittivity_object_2D, &mo);

    struct fdtd_source src = gaussian_source(
        float_cst(30.) * fdtd.dt, float_cst(15.) * fdtd.dt, float_cst(1.));
    add_source_fdtd_2D(source_electric, &fdtd, src,
                       fdtd.domain_size[0] / float_cst(2.),
                       fdtd.domain_size[1] / float_cst(2.));
    struct fdtd retval = {.twoDims = fdtd, .type = fdtd_two_dims};
    return retval;
  } break;
  case west_air_east_water_west_gaussian_pulse_centered_2D: {
    enum border_condition bc[num_borders_2D] = {
        [border_south] = border_perfect_electric_conductor | border_cpml,
        [border_north] = border_perfect_electric_conductor | border_cpml,
        [border_east] = border_perfect_electric_conductor,
        [border_west] = border_perfect_electric_conductor};
    struct fdtd2D fdtd = init_fdtd_2D_cpml(domain_size, Sc, smallest_wavelength,
                                           bc, cpml_thickness);
    struct middle_object_2D mo = {
        .permittivity_medium = float_cst(1.00058986),
        .permittivity_object = float_cst(1.77),
        .permeability_medium = float_cst(1.00000037),
        .permeability_object = float_cst(0.999992),
        .object_center = {domain_size[0] / float_cst(2.),
                          float_cst(3.) * domain_size[1] / float_cst(4.)},
        .object_dimensions = {domain_size[0] * float_cst(2.),
                              domain_size[1] / float_cst(2.)}};
    init_fdtd_2D_medium(&fdtd, init_permeability_object_2D,
                        init_permittivity_object_2D, &mo);

    struct fdtd_source src = gaussian_source(
        float_cst(30.) * fdtd.dt, float_cst(15.) * fdtd.dt, float_cst(1000.));
    add_source_fdtd_2D(source_electric, &fdtd, src,
                       fdtd.domain_size[0] / float_cst(2.), fdtd.dy);
    struct fdtd retval = {.twoDims = fdtd, .type = fdtd_two_dims};
    return retval;
  } break;
  default:
    fprintf(stderr, "The specified 1D setup ID is does not exist\n");
    exit(EXIT_FAILURE);
  }
}

static struct fdtd initializeFdtd3D(unsigned setupID, float_type *domain_size,
                                    float_type Sc,
                                    float_type smallest_wavelength,
                                    uintmax_t cpml_thickness) {
  enum setup3D sID = (enum setup3D)setupID;
  switch (sID) {
  case air_with_object_of_high_permitivity_half_height_centered_3D: {
    enum border_condition bc[num_borders_3D] = {
        [border_front] = border_perfect_electric_conductor,
        [border_back] = border_perfect_electric_conductor,
        [border_bottom] = border_perfect_electric_conductor,
        [border_top] = border_perfect_electric_conductor,
        [border_left] = border_perfect_electric_conductor,
        [border_right] = border_perfect_electric_conductor};
    struct fdtd3D fdtd = init_fdtd_3D_cpml(domain_size, Sc, smallest_wavelength,
                                           bc, cpml_thickness);
    struct middle_object_3D mo = {
        .permittivity_medium = float_cst(1.00058986),
        .permittivity_object = float_cst(1e9),
        .permeability_medium = float_cst(1.00000037),
        .permeability_object = float_cst(1.),
        .object_center = {domain_size[0] / float_cst(2.),
                          domain_size[1] / float_cst(2.),
                          domain_size[2] / float_cst(2.)},
        .object_dimensions = {domain_size[1] / float_cst(2.),
                              domain_size[1] / float_cst(2.),
                              domain_size[1] / float_cst(2.)}};
    init_fdtd_3D_medium(&fdtd, init_permeability_object_3D,
                        init_permittivity_object_3D, &mo);

    struct fdtd_source src = gaussian_source(
        float_cst(10.) * fdtd.dt, float_cst(5.) * fdtd.dt, float_cst(1.e-2));
    for (uintmax_t j = 0; j < fdtd.sizeY; ++j) {
      add_source_fdtd_3D(
          source_magnetic, &fdtd, src, (cpml_thickness + 2) * fdtd.dx,
          (float_type)j * fdtd.dy, (cpml_thickness + 2) * fdtd.dz);
    }

    struct fdtd retval = {.threeDims = fdtd, .type = fdtd_three_dims};
    return retval;
  }
  case half_air_half_water_3D: {
    enum border_condition bc[num_borders_3D] = {
        [border_front] = border_perfect_electric_conductor,
        [border_back] = border_perfect_electric_conductor,
        [border_bottom] = border_perfect_electric_conductor,
        [border_top] = border_perfect_electric_conductor,
        [border_left] = border_perfect_electric_conductor,
        [border_right] = border_perfect_electric_conductor};
    struct fdtd3D fdtd = init_fdtd_3D_cpml(domain_size, Sc, smallest_wavelength,
                                           bc, cpml_thickness);
    struct middle_object_3D mo = {
        .permittivity_medium = float_cst(1.00058986),
        .permittivity_object = float_cst(1.77),
        .permeability_medium = float_cst(1.00000037),
        .permeability_object = float_cst(0.999992),
        .object_center = {domain_size[0] / float_cst(2.),
                          domain_size[1] / float_cst(2.),
                          float_cst(3.) * domain_size[2] / float_cst(4.)},
        .object_dimensions = {domain_size[0] * float_cst(2.),
                              domain_size[1] * float_cst(2.),
                              domain_size[1] / float_cst(2.)}};
    init_fdtd_3D_medium(&fdtd, init_permeability_object_3D,
                        init_permittivity_object_3D, &mo);

    struct fdtd_source src = gaussian_source(
        float_cst(10.) * fdtd.dt, float_cst(5.) * fdtd.dt, float_cst(1.e-2));
    for (uintmax_t j = 0; j < fdtd.sizeY; ++j) {
      add_source_fdtd_3D(
          source_magnetic, &fdtd, src, (cpml_thickness + 2) * fdtd.dx,
          (float_type)j * fdtd.dy, (cpml_thickness + 2) * fdtd.dz);
    }
    struct fdtd retval = {.threeDims = fdtd, .type = fdtd_three_dims};
    return retval;
  } break;
  default:
    fprintf(stderr, "The specified 1D setup ID is does not exist\n");
    exit(EXIT_FAILURE);
  }
}

struct fdtd initializeFdtd(unsigned setupID, float_type *domain_size,
                           float_type Sc, float_type smallest_wavelength) {
  return initializeFdtd_cmpl(setupID, domain_size, Sc, smallest_wavelength, 0);
}

struct fdtd initializeFdtd_cmpl(unsigned setupID, float_type *domain_size,
                                float_type Sc, float_type smallest_wavelength,
                                uintmax_t cpml_thickness) {
  if (setupID < last_1D_setup) { // 1D
    return initializeFdtd1D(setupID, domain_size[0], Sc, smallest_wavelength);
  } else if (setupID < last_2D_setup) { // 2D
    return initializeFdtd2D(setupID, domain_size, Sc, smallest_wavelength,
                            cpml_thickness);
  } else if (setupID < last_3D_setup) { // 3D
    return initializeFdtd3D(setupID, domain_size, Sc, smallest_wavelength,
                            cpml_thickness);
  } else {
    fprintf(stderr, "This setup does not exist\n");
    exit(EXIT_SUCCESS);
  }
}
