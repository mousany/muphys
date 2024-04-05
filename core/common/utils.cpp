// ICON
//
// ---------------------------------------------------------------
// Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
// Contact information: icon-model.org
//
// See AUTHORS.TXT for a list of authors
// See LICENSES/ for license information
// SPDX-License-Identifier: BSD-3-Clause
// ---------------------------------------------------------------
//
#include "utils.hpp"
#include <iostream>

#if defined(MU_ENABLE_SEQ)
void utils_muphys::calc_dz(array_1d_t<real_t> &z, array_1d_t<real_t> &dz,
                           size_t &ncells, size_t &nlev) {
  dz.resize(ncells * nlev);
#else
void utils_muphys::calc_dz(real_t *z, std::unique_ptr<real_t[]> &dz,
                           size_t ncells, size_t nlev) {
  dz.reset(new real_t[ncells * nlev]());
#endif

#if defined(MU_ENABLE_SEQ)
  std::vector<std::vector<real_t>> zh(nlev + 1, std::vector<real_t>(ncells));
#else
  std::unique_ptr<real_t[]> zh(new real_t[(nlev + 1) * ncells]());
#endif

  for (size_t i = 0; i < ncells; i++) {
#if defined(MU_ENABLE_SEQ)
    zh[nlev][i]
#else
    zh[(nlev)*ncells + i]
#endif
        = (static_cast<real_t>(3.0) * z[i + (nlev - 1) * (ncells)] -
           z[i + (nlev - 2) * (ncells)]) *
          static_cast<real_t>(0.5);
  }

  // The loop is intentionally i<nlev; since we are using an unsigned integer
  // data type, when i reaches 0, and you try to decrement further, (to -1), it
  // wraps to the maximum value representable by size_t.
  for (size_t i = nlev - 1; i < nlev; --i) {
    for (size_t j = 0; j < ncells; j++) {
#if defined(MU_ENABLE_SEQ)
      zh[i][j] = static_cast<real_t>(2.0) * z[j + (i * ncells)] - zh[i + 1][j];
      dz[i * ncells + j] = -zh[i + 1][j] + zh[i][j];
#else
      zh[i * ncells + j] = static_cast<real_t>(2.0) * z[j + (i * ncells)] -
                           zh[(i + 1) * ncells + j];
      dz[i * ncells + j] = -zh[(i + 1) * ncells + j] + zh[i * ncells + j];
#endif
    }
  }
}
