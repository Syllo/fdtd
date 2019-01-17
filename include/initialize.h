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

#ifndef INITIALIZE_H_
#define INITIALIZE_H_

#include "fdtd.h"
#include "fdtd1D.h"
#include "fdtd2D.h"
#include "fdtd3D.h"
#include "fdtd_common.h"
#include <stdint.h>

enum setup1D {
  half_air_half_water_1D = 0,
  last_1D_setup,
};

enum setup2D {
  west_air_east_water_west_gaussian_pulse_centered_2D = last_1D_setup + 1,
  object_high_permitivity_in_air_west_gaussian_pulse_centered_2D,
  free_space_gaussian_exitation_centered_absorbing_border_2D,
  last_2D_setup,
};

enum setup3D {
  half_air_half_water_3D = last_2D_setup + 1,
  air_with_object_of_high_permitivity_half_height_centered_3D,
  last_3D_setup,
};

struct fdtd initializeFdtd_cmpl(unsigned setupID, float_type *domain_size,
                                float_type Sc, float_type smallest_wavelength,
                                uintmax_t cmpl_thickness);

struct fdtd initializeFdtd(unsigned setupID, float_type *domain_size,
                           float_type Sc, float_type smallest_wavelength);

#endif // INITIALIZE_H_
