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

#include "fdtd_common.h"

struct fdtd_source gaussian_source(float_type delay, float_type peak_time,
                                   float_type peak_val) {
  struct fdtd_source gs = {.type = source_gaussian_pulse,
                           .gaussian_pulse_delay = delay,
                           .gaussian_pulse_peak_time = peak_time,
                           .gaussian_peak_val = peak_val};
  return gs;
}

float_type gaussian_pulse_val(float_type time, struct fdtd_source *src) {
  /*if ((time - src->gaussian_pulse_delay) >=*/
  /*(-src->gaussian_pulse_peak_time / float_cst(2.)) &&*/
  /*(time - src->gaussian_pulse_delay) <=*/
  /*(src->gaussian_pulse_peak_time / float_cst(2.))) {*/
  float_type exponent =
      (time - src->gaussian_pulse_delay) / src->gaussian_pulse_peak_time;
  exponent *= exponent;
  return src->gaussian_peak_val * exp(-exponent);
  /*} else {*/
  /*return 0.;*/
  /*}*/
}

extern inline float_type Kappa(uintmax_t dist_from_border,
                               uintmax_t CPML_region_width);

extern inline float_type sigma(uintmax_t dist_from_border,
                               uintmax_t CPML_region_width,
                               float_type sigma_max);

extern inline float_type alpha(uintmax_t dist_from_border,
                               uintmax_t CPML_region_width,
                               float_type alpha_max);

extern inline float_type b(uintmax_t dist_from_border,
                           uintmax_t CPML_region_width, float_type dt,
                           float_type alpha_max, float_type sigma_max);

extern inline float_type c(uintmax_t dist_from_border,
                           uintmax_t CPML_region_width, float_type dt,
                           float_type alpha_max, float_type sigma_max);

const char *dumpable_data_name[num_dumpable_data] = {
    "Electric field X directed components",
    "Electric field Y directed components",
    "Electric field Z directed components",
    "Magnetic field X directed components",
    "Magnetic field Y directed components",
    "Magnetic field Z directed components",
    "Multiplicative inverse of the permittivity",
    "Multiplicative inverse of the permeability",
};
