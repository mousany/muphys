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

#include "../common/constants.hpp"
#include "../common/types.hpp"
#include <cmath>

namespace property {

/**
 * @brief TODO
 * @param [in] density TODO
 * @param [in] params TODO
 * @return Fall speed
 */
template <typename array_t>
TARGET real_t fall_speed(real_t density, array_t params) {
  return params[0] * pow((density + params[2]), params[1]);
}

#ifdef MU_ENABLE_OMP
template <const size_t idx> TARGET real_t fall_speed(real_t density);

template <> TARGET real_t fall_speed<0>(real_t density) {
  return static_cast<real_t>(14.58) *
         pow((density + static_cast<real_t>(1.0e-12)),
             static_cast<real_t>(0.111));
}

template <> TARGET real_t fall_speed<1>(real_t density) {
  return static_cast<real_t>(1.25) *
         pow((density + static_cast<real_t>(1.0e-12)),
             static_cast<real_t>(0.160));
}

template <> TARGET real_t fall_speed<2>(real_t density) {
  return static_cast<real_t>(57.80) *
         pow((density + static_cast<real_t>(1.0e-12)),
             static_cast<real_t>(static_cast<real_t>(0.5) /
                                 static_cast<real_t>(3.0)));
}

template <> TARGET real_t fall_speed<3>(real_t density) {
  return static_cast<real_t>(12.24) *
         pow((density + static_cast<real_t>(1.0e-08)),
             static_cast<real_t>(0.217));
}
#endif

} // namespace property
