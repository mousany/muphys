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
#pragma once
#include "types.hpp"

namespace utils_muphys {

#if defined(MU_ENABLE_SEQ)
void calc_dz(array_1d_t<real_t> &z, array_1d_t<real_t> &dz, size_t &ncells,
             size_t &nlev);
#elif defined(MU_ENABLE_OMP)
void calc_dz(real_t *z, std::unique_ptr<real_t[]> &dz, size_t ncells,
             size_t nlev);
#else
void calc_dz(real_t *z, real_t *&dz, size_t ncells, size_t nlev);
#endif

} // namespace utils_muphys
