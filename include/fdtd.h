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

#ifndef FDTD_H_
#define FDTD_H_

#include "fdtd1D.h"
#include "fdtd2D.h"
#include "fdtd3D.h"
#include <stdio.h>
#include <stdlib.h>

enum fdtd_type {
  fdtd_one_dim,
  fdtd_two_dims,
  fdtd_three_dims,
};

struct fdtd {
  enum fdtd_type type;
  union {
    struct fdtd1D oneDim;
    struct fdtd2D twoDims;
    struct fdtd3D threeDims;
  };
};

inline void run_fdtd(struct fdtd *fdtd, float_type end_time) {
  switch (fdtd->type) {
  case fdtd_one_dim:
    run_1D_fdtd(&fdtd->oneDim, end_time);
    break;
  case fdtd_two_dims:
    run_2D_fdtd(&fdtd->twoDims, end_time);
    break;
  case fdtd_three_dims:
    run_3D_fdtd(&fdtd->threeDims, end_time);
    break;
  }
}

inline void dump_fdtd(const struct fdtd *fdtd, const char *fileName,
                      enum dumpable_data dd) {
  switch (fdtd->type) {
  case fdtd_one_dim:
    dump_1D_fdtd(&fdtd->oneDim, fileName, dd);
    break;
  case fdtd_two_dims:
    dump_2D_fdtd(&fdtd->twoDims, fileName, dd);
    break;
  case fdtd_three_dims:
    dump_3D_fdtd(&fdtd->threeDims, fileName, dd);
    break;
  }
}

inline void free_fdtd(struct fdtd *fdtd) {
  switch (fdtd->type) {
  case fdtd_one_dim:
    free_1D_fdtd(&fdtd->oneDim);
    break;
  case fdtd_two_dims:
    free_2D_fdtd(&fdtd->twoDims);
    break;
  case fdtd_three_dims:
    free_3D_fdtd(&fdtd->threeDims);
    break;
  }
}

inline float_type get_time_step_fdtd(const struct fdtd *fdtd) {
  switch (fdtd->type) {
  case fdtd_one_dim:
    return fdtd->oneDim.dt;
  case fdtd_two_dims:
    return fdtd->twoDims.dt;
  case fdtd_three_dims:
    return fdtd->threeDims.dt;
  default:
    fprintf(stderr, "get_time_step: unknown fdtd type\n");
    exit(EXIT_FAILURE);
  }
}

#endif // FDTD_H_
