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

#ifndef FDTD_COMMON_H_
#define FDTD_COMMON_H_

#ifdef FDTD_USE_DOUBLE
typedef double float_type;
#define float_cst(a) a
#else
typedef float float_type;
#define float_cst(a) a##f
#endif

#include <inttypes.h>
#include <tgmath.h>

#undef M_PI
#define M_PI float_cst(3.14159265358979323846)

#define mu0 (float_cst(4.) * M_PI * float_cst(1e-7))
#define eps0 (float_cst(625000.) / (float_cst(22468879468420441.) * M_PI))

#define c_light 299792458
#define c_lightf float_cst(c_light.)

// CPML constant 1 <= kappa <= 20
// Bigger values increases reflection at nornal incidence and reduces reflection
// at other angles
// WARNING: setting kappa_max at values other than 1 is not currently supported
#define kappa_max float_cst(1.)
// Typically found that 3 <= polynomial_taper_order <= 4 is nearly optimal for
// most FDTD simulations
#define polynomial_taper_order float_cst(4.)

enum border_condition {
  border_perfect_electric_conductor = 1,
  border_perfect_magnetic_conductor = 1 << 1,
  border_cpml = 1 << 2,
  /*border_warp = 1<<2,*/
};

enum fdtd_source_type {
  source_gaussian_pulse,
};

struct fdtd_source {
  enum fdtd_source_type type;
  union {
    struct {
      float_type gaussian_pulse_delay;
      float_type gaussian_pulse_peak_time;
      float_type gaussian_peak_val;
    };
  };
};

struct fdtd_source gaussian_source(float_type delay, float_type peak_time,
                                   float_type peak_val);

float_type gaussian_pulse_val(float_type time,
                              struct fdtd_source *gaussian_src);

enum source_type {
  source_magnetic,
  source_electric,
};

enum dumpable_data {
  dump_ex,           // Ez
  dump_ey,           // Ez
  dump_ez,           // Ez
  dump_hx,           // Hx
  dump_hy,           // Hy
  dump_hz,           // Hy
  dump_permittivity, // 1 / permittivity
  dump_permeability, // 1 / Permeability
  num_dumpable_data
};

extern const char *dumpable_data_name[num_dumpable_data];

// Distance 0            = interface CPML / simulation medium
// Distance region_width = simulation border
inline float_type Kappa(uintmax_t dist_from_border,
                        uintmax_t CPML_region_width) {
  float_type db = (float_type)dist_from_border;
  float_type cpml_width = (float_type)CPML_region_width;
  return float_cst(1.) + (kappa_max - float_cst(1.)) *
                             pow(db / cpml_width, polynomial_taper_order);
}

inline float_type sigma(uintmax_t dist_from_border, uintmax_t CPML_region_width,
                        float_type sigma_max) {
  float_type db = (float_type)dist_from_border;
  float_type cpml_width = (float_type)CPML_region_width;
  return sigma_max * pow(db / cpml_width, polynomial_taper_order);
}

inline float_type alpha(uintmax_t dist_from_border, uintmax_t CPML_region_width,
                        float_type alpha_max) {
  float_type db = (float_type)dist_from_border;
  float_type cpml_width = (float_type)CPML_region_width;
  return alpha_max *
         pow(float_cst(1.) - db / cpml_width, polynomial_taper_order);
}

inline float_type b(uintmax_t dist_from_border, uintmax_t CPML_region_width,
                    float_type dt, float_type alpha_max, float_type sigma_max) {
  return exp(-dt *
             (sigma(dist_from_border, CPML_region_width, sigma_max) /
                  (eps0 * Kappa(dist_from_border, CPML_region_width)) +
              alpha(dist_from_border, CPML_region_width, alpha_max) / eps0));
}

inline float_type c(uintmax_t dist_from_border, uintmax_t CPML_region_width,
                    float_type dt, float_type alpha_max, float_type sigma_max) {
  return (sigma(dist_from_border, CPML_region_width, sigma_max) /
          (sigma(dist_from_border, CPML_region_width, sigma_max) *
               Kappa(dist_from_border, CPML_region_width) +
           Kappa(dist_from_border, CPML_region_width) *
               Kappa(dist_from_border, CPML_region_width) *
               alpha(dist_from_border, CPML_region_width, alpha_max))) *
         (b(dist_from_border, CPML_region_width, dt, alpha_max, sigma_max) -
          float_cst(1.));
}

#endif // FDTD_COMMON_H_
