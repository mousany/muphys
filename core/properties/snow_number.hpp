#pragma once

#include "../common/constants.hpp"
#include "../common/types.hpp"
#include <cmath>

namespace property {

constexpr real_t tmin = thermodyn::tmelt - static_cast<real_t>(40.);
constexpr real_t tmax = thermodyn::tmelt;
constexpr real_t qsmin = 2.0e-6;
constexpr real_t xa1 = -1.65e+0;
constexpr real_t xa2 = 5.45e-2;
constexpr real_t xa3 = 3.27e-4;
constexpr real_t xb1 = 1.42e+0;
constexpr real_t xb2 = 1.19e-2;
constexpr real_t xb3 = 9.60e-5;
constexpr real_t n0s0 = 8.00e+5;
constexpr real_t n0s1 =
    static_cast<real_t>(13.5) * static_cast<real_t>(5.65e+05);
constexpr real_t n0s2 = -0.107;
constexpr real_t n0s3 = 13.5;
constexpr real_t n0s4 = static_cast<real_t>(0.5) * n0s1;
constexpr real_t n0s5 = 1.e6;
constexpr real_t n0s6 = static_cast<real_t>(1.e2) * n0s1;
constexpr real_t n0s7 = 1.e9;

/**
 * @brief TODO
 * @param [in] t Temperature
 * @param [in] rho Ambient air density
 * @param [in] qs Snow specific mass
 * @return Snow number
 */
TARGET real_t snow_number(real_t t, real_t rho, real_t qs) {

  if (qs > graupel_ct::qmin) {
    real_t tc = fmax(fmin(t, tmax), tmin) - thermodyn::tmelt;
    real_t alf = pow(static_cast<real_t>(10.), (xa1 + tc * (xa2 + tc * xa3)));
    real_t bet = xb1 + tc * (xb2 + tc * xb3);
    real_t n0s =
        n0s3 *
        pow(((qs + qsmin) * rho / graupel_ct::ams),
            (static_cast<real_t>(4.0) - static_cast<real_t>(3) * bet)) /
        (alf * alf * alf);
    real_t y = exp(n0s2 * tc);
    real_t n0smn = fmax(n0s4 * y, n0s5);
    real_t n0smx = fmin(n0s6 * y, n0s7);
    return fmin(n0smx, fmax(n0smn, n0s));
  } else {
    return n0s0;
  }
}
} // namespace property
