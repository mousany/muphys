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

#if defined(MU_ENABLE_OMP)

#include "core/common/graupel.hpp"
#include <algorithm>
#include <cstring>
#include <iostream>
#include <numeric>

using namespace property;
using namespace thermo;
using namespace transition;
using namespace idx;
using namespace graupel_ct;

// struct t_qx_ptr {
//   t_qx_ptr(array_1d_t<real_t> &p_, array_1d_t<real_t> &x_) : p(p_), x(x_) {}

//   array_1d_t<real_t> &p;
//   array_1d_t<real_t> &x;
// }; // pointer vector

struct Properties {
  real_t qv;
  real_t qc;
  real_t qi;
  real_t qr;
  real_t qs;
  real_t qg;
};

/**
 * @brief TODO
 *
 * @param precip time step for integration of microphysics  (s)
 * @param zeta dt/(2dz)
 * @param vc state dependent fall speed correction
 * @param flx flux into cell from above
 * @param vt terminal velocity
 * @param q specific mass of hydrometeor
 * @param q_kp1 specific mass in next lower cell
 * @param rho density
 */
template <const size_t idx>
TARGET void precip(real_t &precip_0, real_t &precip_1, real_t &precip_2,
                   real_t zeta, real_t vc, real_t flx, real_t vt, real_t q,
                   real_t q_kp1, real_t rho) {
  real_t rho_x, flx_eff, flx_partial;
  rho_x = q * rho;
  flx_eff = (rho_x / zeta) + static_cast<real_t>(2.0) * flx;
  flx_partial = rho_x * vc * fall_speed<idx>(rho_x);
  flx_partial = std::fmin(flx_partial, flx_eff);
  precip_0 = (zeta * (flx_eff - flx_partial)) /
             ((static_cast<real_t>(1.0) + zeta * vt) * rho); // q update
  precip_1 =
      (precip_0 * rho * vt + flx_partial) * static_cast<real_t>(0.5); // flx
  rho_x = (precip_0 + q_kp1) * static_cast<real_t>(0.5) * rho;
  precip_2 = vc * fall_speed<idx>(rho_x); // vt
}

void pack_props(size_t nvec, size_t ke, std::unique_ptr<Properties[]> &props,
                size_t ivstart, size_t ivend, size_t kstart, real_t *qv,
                real_t *qc, real_t *qi, real_t *qr, real_t *qs, real_t *qg) {
  props.reset(new (std::align_val_t(64)) Properties[nvec * ke]);

  const size_t i_start = kstart * nvec;
  const size_t i_end = ke * nvec;

#pragma omp parallel for
  for (size_t oned_vec_index = i_start; oned_vec_index < i_end;
       oned_vec_index++) {
    props[oned_vec_index].qv = qv[oned_vec_index];
    props[oned_vec_index].qc = qc[oned_vec_index];
    props[oned_vec_index].qi = qi[oned_vec_index];
    props[oned_vec_index].qr = qr[oned_vec_index];
    props[oned_vec_index].qs = qs[oned_vec_index];
    props[oned_vec_index].qg = qg[oned_vec_index];
  }
}

void unpack_props(size_t nvec, size_t ke, std::unique_ptr<Properties[]> &props,
                  size_t ivstart, size_t ivend, size_t kstart, real_t *qv,
                  real_t *qc, real_t *qi, real_t *qr, real_t *qs, real_t *qg) {
  const size_t i_start = kstart * nvec;
  const size_t i_end = ke * nvec;
#pragma omp parallel for
  for (size_t oned_vec_index = i_start; oned_vec_index < i_end;
       oned_vec_index++) {
    qv[oned_vec_index] = props[oned_vec_index].qv;
    qc[oned_vec_index] = props[oned_vec_index].qc;
    qi[oned_vec_index] = props[oned_vec_index].qi;
    qr[oned_vec_index] = props[oned_vec_index].qr;
    qs[oned_vec_index] = props[oned_vec_index].qs;
    qg[oned_vec_index] = props[oned_vec_index].qg;
  }
}
// data[r, c, k]
// r: R, c: C from input, k: K constant
// r linear increase, c any

// Loop 1: iterate r, c
//  1. for each c, find the smallest r (as r*) that satisfies condition P1
//      P1: data[r, c, k] > Constant
//  2. find all r, c, that satisfies P1 and P2
//      P2: dynamically calculated from data[r, c, k]

// Loop 2: for all r, c found in Loop 1
//  1. transform data[r, c, k] by F
//      F only depends on data[r, c, k]

// Loop 3: iterate r, c
//  1. for each c, where r >= r*, update data[r, c, k] by G
//      G depends on an internal state from r-1, c and data[r+1, c, k]

TARGET void solidify(size_t oned_vec_index, bool is_sig_present, real_t dt,
                     real_t *dz, real_t *t, real_t *rho, real_t *p,
                     Properties *props, real_t qnc) {

  real_t sx2x[nx * nx] = {0.}; // conversion rates

  real_t dvsw = props[oned_vec_index].qv -
                qsat_rho(t[oned_vec_index], rho[oned_vec_index]);
  real_t qvsi = qsat_ice_rho(t[oned_vec_index], rho[oned_vec_index]);
  real_t dvsi = props[oned_vec_index].qv - qvsi;
  real_t n_snow = snow_number(t[oned_vec_index], rho[oned_vec_index],
                              props[oned_vec_index].qs);
  real_t l_snow =
      snow_lambda(rho[oned_vec_index], props[oned_vec_index].qs, n_snow);

  sx2x[lqc * nx + lqr] =
      cloud_to_rain(t[oned_vec_index], props[oned_vec_index].qc,
                    props[oned_vec_index].qr, qnc);
  sx2x[lqr * nx + lqv] = rain_to_vapor(t[oned_vec_index], rho[oned_vec_index],
                                       props[oned_vec_index].qc,
                                       props[oned_vec_index].qr, dvsw, dt);
  sx2x[lqc * nx + lqi] =
      cloud_x_ice(t[oned_vec_index], props[oned_vec_index].qc,
                  props[oned_vec_index].qi, dt);
  sx2x[lqi * nx + lqc] = -std::fmin(sx2x[lqc * nx + lqi], 0.0);
  sx2x[lqc * nx + lqi] = std::fmax(sx2x[lqc * nx + lqi], 0.0);
  sx2x[lqc * nx + lqs] =
      cloud_to_snow(t[oned_vec_index], props[oned_vec_index].qc,
                    props[oned_vec_index].qs, n_snow, l_snow);
  sx2x[lqc * nx + lqg] =
      cloud_to_graupel(t[oned_vec_index], rho[oned_vec_index],
                       props[oned_vec_index].qc, props[oned_vec_index].qg);

  real_t ice_dep = 0.0;
  real_t eta = 0.0;

  if (t[oned_vec_index] < tmelt) {
    real_t n_ice = ice_number(t[oned_vec_index], rho[oned_vec_index]);
    real_t m_ice = ice_mass(props[oned_vec_index].qi, n_ice);
    real_t x_ice = ice_sticking(t[oned_vec_index]);

    if (is_sig_present) {
      eta =
          deposition_factor(t[oned_vec_index],
                            qvsi); // neglect cloud depth cor. from gcsp_graupel
      sx2x[lqv * nx + lqi] =
          vapor_x_ice(props[oned_vec_index].qi, m_ice, eta, dvsi, dt);
      sx2x[lqi * nx + lqv] = -std::fmin(sx2x[lqv * nx + lqi], 0.0);
      sx2x[lqv * nx + lqi] = std::fmax(sx2x[lqv * nx + lqi], 0.0);
      ice_dep = std::fmin(sx2x[lqv * nx + lqi], dvsi / dt);

      sx2x[lqi * nx + lqs] =
          deposition_auto_conversion(props[oned_vec_index].qi, m_ice, ice_dep);
      sx2x[lqi * nx + lqs] =
          sx2x[lqi * nx + lqs] +
          ice_to_snow(props[oned_vec_index].qi, n_snow, l_snow, x_ice);
      sx2x[lqi * nx + lqg] = ice_to_graupel(
          rho[oned_vec_index], props[oned_vec_index].qr,
          props[oned_vec_index].qg, props[oned_vec_index].qi, x_ice);
      sx2x[lqs * nx + lqg] =
          snow_to_graupel(t[oned_vec_index], rho[oned_vec_index],
                          props[oned_vec_index].qc, props[oned_vec_index].qs);
      sx2x[lqr * nx + lqg] = rain_to_graupel(
          t[oned_vec_index], rho[oned_vec_index], props[oned_vec_index].qc,
          props[oned_vec_index].qr, props[oned_vec_index].qi,
          props[oned_vec_index].qs, m_ice, dvsw, dt);
    }
    sx2x[lqv * nx + lqi] =
        sx2x[lqv * nx + lqi] +
        ice_deposition_nucleation(t[oned_vec_index], props[oned_vec_index].qc,
                                  props[oned_vec_index].qi, n_ice, dvsi, dt);
  } else {
    sx2x[lqc * nx + lqr] =
        sx2x[lqc * nx + lqr] + sx2x[lqc * nx + lqs] + sx2x[lqc * nx + lqg];
    sx2x[lqc * nx + lqs] = 0.0;
    sx2x[lqc * nx + lqg] = 0.0;
    ice_dep = 0.0;
    eta = 0.0;
  }

  if (is_sig_present) {
    real_t dvsw0 =
        props[oned_vec_index].qv - qsat_rho(tmelt, rho[oned_vec_index]);
    sx2x[lqv * nx + lqs] =
        vapor_x_snow(t[oned_vec_index], p[oned_vec_index], rho[oned_vec_index],
                     props[oned_vec_index].qs, n_snow, l_snow, eta, ice_dep,
                     dvsw, dvsi, dvsw0, dt);
    sx2x[lqs * nx + lqv] = -std::fmin(sx2x[lqv * nx + lqs], 0.0);
    sx2x[lqv * nx + lqs] = std::fmax(sx2x[lqv * nx + lqs], 0.0);
    sx2x[lqv * nx + lqg] = vapor_x_graupel(
        t[oned_vec_index], p[oned_vec_index], rho[oned_vec_index],
        props[oned_vec_index].qg, dvsw, dvsi, dvsw0, dt);
    sx2x[lqg * nx + lqv] = -std::fmin(sx2x[lqv * nx + lqg], 0.0);
    sx2x[lqv * nx + lqg] = std::fmax(sx2x[lqv * nx + lqg], 0.0);
    sx2x[lqs * nx + lqr] =
        snow_to_rain(t[oned_vec_index], p[oned_vec_index], rho[oned_vec_index],
                     dvsw0, props[oned_vec_index].qs);
    sx2x[lqg * nx + lqr] =
        graupel_to_rain(t[oned_vec_index], p[oned_vec_index],
                        rho[oned_vec_index], dvsw0, props[oned_vec_index].qg);
  }

  real_t stot = 0.0;
  real_t sink[nx]; // tendencies
  real_t dqdt[nx]; // tendencies

  // qx_ind = {5, 4, 0, 2, 1, 3}
  // ix = 0, qx_ind[0] = 5
  sink[5] = 0.0;
  if ((is_sig_present) or (5 == lqc) or (5 == lqv) or (5 == lqr)) {
    sink[5] = sx2x[5 * nx + 0] + sx2x[5 * nx + 1] + sx2x[5 * nx + 2] +
              sx2x[5 * nx + 3] + sx2x[5 * nx + 4] + sx2x[5 * nx + 5];
    stot = props[oned_vec_index].qv / dt;
    if ((sink[5] > stot) && (props[oned_vec_index].qv > qmin)) {
#pragma omp simd
      for (size_t i = 0; i < nx; i++) {
        sx2x[5 * nx + i] = sx2x[5 * nx + i] * stot / sink[5];
      }
      sink[5] = sx2x[5 * nx + 0] + sx2x[5 * nx + 1] + sx2x[5 * nx + 2] +
                sx2x[5 * nx + 3] + sx2x[5 * nx + 4] + sx2x[5 * nx + 5];
    }
  }

  // ix = 1, qx_ind[1] = 4
  sink[4] = 0.0;
  if ((is_sig_present) or (4 == lqc) or (4 == lqv) or (4 == lqr)) {
    sink[4] = sx2x[4 * nx + 0] + sx2x[4 * nx + 1] + sx2x[4 * nx + 2] +
              sx2x[4 * nx + 3] + sx2x[4 * nx + 4] + sx2x[4 * nx + 5];
    stot = props[oned_vec_index].qc / dt;
    if ((sink[4] > stot) && (props[oned_vec_index].qc > qmin)) {
#pragma omp simd
      for (size_t i = 0; i < nx; i++) {
        sx2x[4 * nx + i] = sx2x[4 * nx + i] * stot / sink[4];
      }
      sink[4] = sx2x[4 * nx + 0] + sx2x[4 * nx + 1] + sx2x[4 * nx + 2] +
                sx2x[4 * nx + 3] + sx2x[4 * nx + 4] + sx2x[4 * nx + 5];
    }
  }

  // ix = 2, qx_ind[2] = 0
  sink[0] = 0.0;
  if ((is_sig_present) or (0 == lqc) or (0 == lqv) or (0 == lqr)) {
    sink[0] = sx2x[0 * nx + 0] + sx2x[0 * nx + 1] + sx2x[0 * nx + 2] +
              sx2x[0 * nx + 3] + sx2x[0 * nx + 4] + sx2x[0 * nx + 5];
    stot = props[oned_vec_index].qr / dt;
    if ((sink[0] > stot) && (props[oned_vec_index].qr > qmin)) {
#pragma omp simd
      for (size_t i = 0; i < nx; i++) {
        sx2x[0 * nx + i] = sx2x[0 * nx + i] * stot / sink[0];
      }
      sink[0] = sx2x[0 * nx + 0] + sx2x[0 * nx + 1] + sx2x[0 * nx + 2] +
                sx2x[0 * nx + 3] + sx2x[0 * nx + 4] + sx2x[0 * nx + 5];
    }
  }

  // ix = 3, qx_ind[3] = 2
  sink[2] = 0.0;
  if ((is_sig_present) or (2 == lqc) or (2 == lqv) or (2 == lqr)) {
    sink[2] = sx2x[2 * nx + 0] + sx2x[2 * nx + 1] + sx2x[2 * nx + 2] +
              sx2x[2 * nx + 3] + sx2x[2 * nx + 4] + sx2x[2 * nx + 5];
    stot = props[oned_vec_index].qs / dt;
    if ((sink[2] > stot) && (props[oned_vec_index].qs > qmin)) {
#pragma omp simd
      for (size_t i = 0; i < nx; i++) {
        sx2x[2 * nx + i] = sx2x[2 * nx + i] * stot / sink[2];
      }
      sink[2] = sx2x[2 * nx + 0] + sx2x[2 * nx + 1] + sx2x[2 * nx + 2] +
                sx2x[2 * nx + 3] + sx2x[2 * nx + 4] + sx2x[2 * nx + 5];
    }
  }

  // ix = 4, qx_ind[4] = 1
  sink[1] = 0.0;
  if ((is_sig_present) or (1 == lqc) or (1 == lqv) or (1 == lqr)) {
    sink[1] = sx2x[1 * nx + 0] + sx2x[1 * nx + 1] + sx2x[1 * nx + 2] +
              sx2x[1 * nx + 3] + sx2x[1 * nx + 4] + sx2x[1 * nx + 5];
    stot = props[oned_vec_index].qi / dt;
    if ((sink[1] > stot) && (props[oned_vec_index].qi > qmin)) {
#pragma omp simd
      for (size_t i = 0; i < nx; i++) {
        sx2x[1 * nx + i] = sx2x[1 * nx + i] * stot / sink[1];
      }
      sink[1] = sx2x[1 * nx + 0] + sx2x[1 * nx + 1] + sx2x[1 * nx + 2] +
                sx2x[1 * nx + 3] + sx2x[1 * nx + 4] + sx2x[1 * nx + 5];
    }
  }

  // ix = 5, qx_ind[5] = 3
  sink[3] = 0.0;
  if ((is_sig_present) or (3 == lqc) or (3 == lqv) or (3 == lqr)) {
    sink[3] = sx2x[3 * nx + 0] + sx2x[3 * nx + 1] + sx2x[3 * nx + 2] +
              sx2x[3 * nx + 3] + sx2x[3 * nx + 4] + sx2x[3 * nx + 5];
    stot = props[oned_vec_index].qg / dt;
    if ((sink[3] > stot) && (props[oned_vec_index].qg > qmin)) {
#pragma omp simd
      for (size_t i = 0; i < nx; i++) {
        sx2x[3 * nx + i] = sx2x[3 * nx + i] * stot / sink[3];
      }
      sink[3] = sx2x[3 * nx + 0] + sx2x[3 * nx + 1] + sx2x[3 * nx + 2] +
                sx2x[3 * nx + 3] + sx2x[3 * nx + 4] + sx2x[3 * nx + 5];
    }
  }

  // qx_ind = {5, 4, 0, 2, 1, 3}
  // ix = 0, qx_ind[0] = 5
  dqdt[5] = sx2x[0 * nx + 5] + sx2x[1 * nx + 5] + sx2x[2 * nx + 5] +
            sx2x[3 * nx + 5] + sx2x[4 * nx + 5] + sx2x[5 * nx + 5] - sink[5];
  props[oned_vec_index].qv =
      std::fmax(0.0, props[oned_vec_index].qv + dqdt[5] * dt);

  // ix = 1, qx_ind[1] = 4
  dqdt[4] = sx2x[0 * nx + 4] + sx2x[1 * nx + 4] + sx2x[2 * nx + 4] +
            sx2x[3 * nx + 4] + sx2x[4 * nx + 4] + sx2x[5 * nx + 4] - sink[4];
  props[oned_vec_index].qc =
      std::fmax(0.0, props[oned_vec_index].qc + dqdt[4] * dt);

  // ix = 2, qx_ind[2] = 0
  dqdt[0] = sx2x[0 * nx + 0] + sx2x[1 * nx + 0] + sx2x[2 * nx + 0] +
            sx2x[3 * nx + 0] + sx2x[4 * nx + 0] + sx2x[5 * nx + 0] - sink[0];
  props[oned_vec_index].qr =
      std::fmax(0.0, props[oned_vec_index].qr + dqdt[0] * dt);

  // ix = 3, qx_ind[3] = 2
  dqdt[2] = sx2x[0 * nx + 2] + sx2x[1 * nx + 2] + sx2x[2 * nx + 2] +
            sx2x[3 * nx + 2] + sx2x[4 * nx + 2] + sx2x[5 * nx + 2] - sink[2];
  props[oned_vec_index].qs =
      std::fmax(0.0, props[oned_vec_index].qs + dqdt[2] * dt);

  // ix = 4, qx_ind[4] = 1
  dqdt[1] = sx2x[0 * nx + 1] + sx2x[1 * nx + 1] + sx2x[2 * nx + 1] +
            sx2x[3 * nx + 1] + sx2x[4 * nx + 1] + sx2x[5 * nx + 1] - sink[1];
  props[oned_vec_index].qi =
      std::fmax(0.0, props[oned_vec_index].qi + dqdt[1] * dt);

  // ix = 5, qx_ind[5] = 3
  dqdt[3] = sx2x[0 * nx + 3] + sx2x[1 * nx + 3] + sx2x[2 * nx + 3] +
            sx2x[3 * nx + 3] + sx2x[4 * nx + 3] + sx2x[5 * nx + 3] - sink[3];
  props[oned_vec_index].qg =
      std::fmax(0.0, props[oned_vec_index].qg + dqdt[3] * dt);

  real_t qice = props[oned_vec_index].qs + props[oned_vec_index].qi +
                props[oned_vec_index].qg;
  real_t qliq = props[oned_vec_index].qc + props[oned_vec_index].qr;
  real_t qtot = props[oned_vec_index].qv + qice + qliq;
  real_t cv = cvd + (cvv - cvd) * qtot + (clw - cvv) * qliq +
              (ci - cvv) * qice; // qtot? or qv?
  t[oned_vec_index] =
      t[oned_vec_index] +
      dt *
          ((dqdt[lqc] + dqdt[lqr]) * (lvc - (clw - cvv) * t[oned_vec_index]) +
           (dqdt[lqi] + dqdt[lqs] + dqdt[lqg]) *
               (lsc - (ci - cvv) * t[oned_vec_index])) /
          cv;
}

void graupel(size_t nvec, size_t ke, size_t ivstart, size_t ivend,
             size_t kstart, real_t dt, real_t *dz, real_t *t, real_t *rho,
             real_t *p, real_t *qv, real_t *qc, real_t *qi, real_t *qr,
             real_t *qs, real_t *qg, real_t qnc) {
  // std::cout << "sequential graupel" << std::endl;

  std::unique_ptr<Properties[]> props;
  pack_props(nvec, ke, props, ivstart, ivend, kstart, qv, qc, qi, qr, qs, qg);

  // nvec = ncells
  // ivbeg = 0, ivend = ncells
  // kbeg = 0, kend = nlev

  // before:
  // |------ ncells ------|
  // |
  // nlev
  // |
  // |

  // array_1d_t<t_qx_ptr>
  //     q{}; // vector of pointers to point to four hydrometeor inouts
  // array_1d_t<real_t> emptyArray;
  // q.emplace_back(prr_gsp, qr);
  // q.emplace_back(pri_gsp, qi);
  // q.emplace_back(prs_gsp, qs);
  // q.emplace_back(prg_gsp, qg);

  // q.emplace_back(emptyArray, qc);
  // q.emplace_back(emptyArray, qv);

  // |  i  |  q  |    p    |   x  |
  // |  0  | lqr | prr_gsp |  qr  |
  // |  1  | lqi | pri_gsp |  qi  |
  // |  2  | lqs | prs_gsp |  qs  |
  // |  3  | lqg | prg_gsp |  qg  |
  // |  4  | lqc |         |  qc  |
  // |  5  | lqv |         |  qv  |

  // The loop is intentionally i<nlev; since we are using an unsigned integer
  // data type, when i reaches 0, and you try to decrement further, (to -1), it
  // wraps to the maximum value representable by size_t.

  constexpr size_t BLOCK_SIZE = 64 / sizeof(real_t);

  size_t k_end = (lrain) ? ke : kstart - 1;
#pragma omp parallel for
  for (size_t blk = ivstart; blk < ivend; blk += BLOCK_SIZE) {
    size_t blk_end = std::min(ivend, blk + BLOCK_SIZE);

    real_t eflx[BLOCK_SIZE] = {
        0}; // internal energy flux from precipitation (W/m2 )
    real_t vt[BLOCK_SIZE * np] = {0.};   // terminal velocity for different
                                         // hydrometeor categories
    uint8_t kmin_flag[BLOCK_SIZE] = {0}; // flag for kmin

    real_t prg_gsp[BLOCK_SIZE] = {0.};
    real_t prs_gsp[BLOCK_SIZE] = {0.};
    real_t pri_gsp[BLOCK_SIZE] = {0.};
    real_t prr_gsp[BLOCK_SIZE] = {0.};

    for (size_t iv = blk; iv < blk_end; iv++) {
      bool qc_qmin = props[iv].qc > qmin;
      bool qr_qmin = props[iv].qr > qmin;
      bool qs_qmin = props[iv].qs > qmin;
      bool qi_qmin = props[iv].qi > qmin;
      bool qg_qmin = props[iv].qg > qmin;

      kmin_flag[iv - blk] =
          (qr_qmin << 3) | (qi_qmin << 2) | (qs_qmin << 1) | (qg_qmin << 0);

      if ((qc_qmin or qr_qmin or qs_qmin or qi_qmin or qg_qmin) or
          ((t[iv] < tfrz_het2) and
           (props[iv].qv > qsat_ice_rho(t[iv], rho[iv])))) {

        bool is_sig_present =
            qs_qmin or qi_qmin or qg_qmin; // is snow, ice or graupel present?

        solidify(iv, is_sig_present, dt, dz, t, rho, p, props.get(), qnc);
      }
    }

    for (size_t k = kstart; k < k_end; k++) {
      for (size_t iv = blk; iv < blk_end; iv++) {
        bool qr_qmin = false;
        bool qs_qmin = false;
        bool qi_qmin = false;
        bool qg_qmin = false;

        if (k < ke - 1) {
          size_t nexted_vec_index = (k + 1) * ivend + iv;

          bool qc_qmin = props[nexted_vec_index].qc > qmin;
          qr_qmin = props[nexted_vec_index].qr > qmin;
          qs_qmin = props[nexted_vec_index].qs > qmin;
          qi_qmin = props[nexted_vec_index].qi > qmin;
          qg_qmin = props[nexted_vec_index].qg > qmin;

          if ((qc_qmin or qr_qmin or qs_qmin or qi_qmin or qg_qmin) or
              ((t[nexted_vec_index] < tfrz_het2) and
               (props[nexted_vec_index].qv >
                qsat_ice_rho(t[nexted_vec_index], rho[nexted_vec_index])))) {

            bool is_sig_present = qs_qmin or qi_qmin or
                                  qg_qmin; // is snow, ice or graupel present?

            solidify(nexted_vec_index, is_sig_present, dt, dz, t, rho, p,
                     props.get(), qnc);
          }
        }

        size_t oned_vec_index = k * ivend + iv;

        size_t kp1 = std::min(ke - 1, k + 1);
        if (kmin_flag[iv - blk]) {
          real_t qliq = props[oned_vec_index].qc + props[oned_vec_index].qr;
          real_t qice = props[oned_vec_index].qs + props[oned_vec_index].qi +
                        props[oned_vec_index].qg;

          real_t e_int =
              internal_energy(t[oned_vec_index], props[oned_vec_index].qv, qliq,
                              qice, rho[oned_vec_index], dz[oned_vec_index]) +
              eflx[iv - blk];
          real_t zeta = dt / (2.0 * dz[oned_vec_index]);
          real_t xrho = std::sqrt(rho_00 / rho[oned_vec_index]);

          // real_t update[3]; // scratch array with output from precipitation
          // step

          // qp_ind = {0, 1, 2, 3}
          // ix = 0, qp_ind[0] = 0
          if (kmin_flag[iv - blk] & 0b1000) {
            real_t vc =
                vel_scale_factor(0, xrho, rho[oned_vec_index],
                                 t[oned_vec_index], props[oned_vec_index].qr);
            precip<0>(props[oned_vec_index].qr, prr_gsp[iv - blk],
                      vt[(iv - blk) * np], zeta, vc, prr_gsp[iv - blk],
                      vt[(iv - blk) * np], props[oned_vec_index].qr,
                      props[kp1 * ivend + iv].qr, rho[oned_vec_index]);
          }

          // ix = 1, qp_ind[1] = 1
          if (kmin_flag[iv - blk] & 0b0100) {
            real_t vc =
                vel_scale_factor(1, xrho, rho[oned_vec_index],
                                 t[oned_vec_index], props[oned_vec_index].qi);
            precip<1>(props[oned_vec_index].qi, pri_gsp[iv - blk],
                      vt[(iv - blk) * np + 1], zeta, vc, pri_gsp[iv - blk],
                      vt[(iv - blk) * np + 1], props[oned_vec_index].qi,
                      props[kp1 * ivend + iv].qi, rho[oned_vec_index]);
          }

          // ix = 2, qp_ind[2] = 2
          if (kmin_flag[iv - blk] & 0b0010) {
            real_t vc =
                vel_scale_factor(2, xrho, rho[oned_vec_index],
                                 t[oned_vec_index], props[oned_vec_index].qs);
            precip<2>(props[oned_vec_index].qs, prs_gsp[iv - blk],
                      vt[(iv - blk) * np + 2], zeta, vc, prs_gsp[iv - blk],
                      vt[(iv - blk) * np + 2], props[oned_vec_index].qs,
                      props[kp1 * ivend + iv].qs, rho[oned_vec_index]);
          }

          // ix = 3, qp_ind[3] = 3
          if (kmin_flag[iv - blk] & 0b0001) {
            real_t vc =
                vel_scale_factor(3, xrho, rho[oned_vec_index],
                                 t[oned_vec_index], props[oned_vec_index].qg);
            precip<3>(props[oned_vec_index].qg, prg_gsp[iv - blk],
                      vt[(iv - blk) * np + 3], zeta, vc, prg_gsp[iv - blk],
                      vt[(iv - blk) * np + 3], props[oned_vec_index].qg,
                      props[kp1 * ivend + iv].qg, rho[oned_vec_index]);
          }

          real_t pflx =
              prs_gsp[iv - blk] + pri_gsp[iv - blk] + prg_gsp[iv - blk];
          eflx[iv - blk] =
              dt * (prr_gsp[iv - blk] * (clw * t[oned_vec_index] -
                                         cvd * t[kp1 * ivend + iv] - lvc) +
                    pflx * (ci * t[oned_vec_index] - cvd * t[kp1 * ivend + iv] -
                            lsc));
          qliq = props[oned_vec_index].qc + props[oned_vec_index].qr;
          qice = props[oned_vec_index].qs + props[oned_vec_index].qi +
                 props[oned_vec_index].qg;
          e_int = e_int - eflx[iv - blk];
          t[oned_vec_index] = T_from_internal_energy(
              e_int, props[oned_vec_index].qv, qliq, qice, rho[oned_vec_index],
              dz[oned_vec_index]);
        }

        kmin_flag[iv - blk] |=
            (qr_qmin << 3) | (qi_qmin << 2) | (qs_qmin << 1) | (qg_qmin << 0);
      }
    }
  }

  unpack_props(nvec, ke, props, ivstart, ivend, kstart, qv, qc, qi, qr, qs, qg);
}

#endif // MU_ENABLE_OMP
