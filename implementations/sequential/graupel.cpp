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
#include "core/common/graupel.hpp"
#include <algorithm>
#include <iostream>
#include <numeric>

using namespace property;
using namespace thermo;
using namespace transition;
using namespace idx;
using namespace graupel_ct;

struct t_qx_ptr {
  t_qx_ptr(array_1d_t<real_t> &p_, array_1d_t<real_t> &x_) : p(p_), x(x_) {}

  array_1d_t<real_t> &p;
  array_1d_t<real_t> &x;
}; // pointer vector

/**
 * @brief TODO
 *
 * @param precip time step for integration of microphysics  (s)
 * @param params fall speed parameters
 * @param zeta dt/(2dz)
 * @param vc state dependent fall speed correction
 * @param flx flux into cell from above
 * @param vt terminal velocity
 * @param q specific mass of hydrometeor
 * @param q_kp1 specific mass in next lower cell
 * @param rho density
 */
void precip(const real_t (&params)[3], real_t (&precip)[3], real_t zeta,
            real_t vc, real_t flx, real_t vt, real_t q, real_t q_kp1,
            real_t rho) {
  real_t rho_x, flx_eff, flx_partial;
  rho_x = q * rho;
  flx_eff = (rho_x / zeta) + static_cast<real_t>(2.0) * flx;
  flx_partial = rho_x * vc * fall_speed(rho_x, params);
  flx_partial = std::fmin(flx_partial, flx_eff);
  precip[0] = (zeta * (flx_eff - flx_partial)) /
              ((static_cast<real_t>(1.0) + zeta * vt) * rho); // q update
  precip[1] =
      (precip[0] * rho * vt + flx_partial) * static_cast<real_t>(0.5); // flx
  rho_x = (precip[0] + q_kp1) * static_cast<real_t>(0.5) * rho;
  precip[2] = vc * fall_speed(rho_x, params); // vt
}

void graupel(size_t &nvec, size_t &ke, size_t &ivstart, size_t &ivend,
             size_t &kstart, real_t &dt, array_1d_t<real_t> &dz,
             array_1d_t<real_t> &t, array_1d_t<real_t> &rho,
             array_1d_t<real_t> &p, array_1d_t<real_t> &qv,
             array_1d_t<real_t> &qc, array_1d_t<real_t> &qi,
             array_1d_t<real_t> &qr, array_1d_t<real_t> &qs,
             array_1d_t<real_t> &qg, real_t &qnc, array_1d_t<real_t> &prr_gsp,
             array_1d_t<real_t> &pri_gsp, array_1d_t<real_t> &prs_gsp,
             array_1d_t<real_t> &prg_gsp, array_1d_t<real_t> &pflx) {
  // std::cout << "sequential graupel" << std::endl;

  array_1d_t<bool> is_sig_present(nvec *
                                  ke); // is snow, ice or graupel present?

  array_1d_t<size_t> ind_k(nvec * ke),
      ind_i(nvec * ke); // k index of gathered point, iv index of gathered point
  array_2d_t<size_t> kmin(
      nvec, array_1d_t<size_t>(np)); // first level with condensate

  real_t cv, vc, eta, zeta, qvsi, qice, qliq, qtot, dvsw, dvsw0, dvsi, n_ice,
      m_ice, x_ice, n_snow, l_snow, ice_dep, e_int, stot, xrho;

  real_t update[3], // scratch array with output from precipitation step
      sink[nx],     // tendencies
      dqdt[nx];     // tendencies
  array_1d_t<real_t> eflx(
      nvec); // internal energy flux from precipitation (W/m2 )
  array_2d_t<real_t> sx2x(nx, array_1d_t<real_t>(nx, 0.0)), // conversion rates
      vt(nvec,
         array_1d_t<real_t>(
             np)); // terminal velocity for different hydrometeor categories

  array_1d_t<t_qx_ptr>
      q{}; // vector of pointers to point to four hydrometeor inouts
  array_1d_t<real_t> emptyArray;
  q.emplace_back(prr_gsp, qr);
  q.emplace_back(pri_gsp, qi);
  q.emplace_back(prs_gsp, qs);
  q.emplace_back(prg_gsp, qg);

  q.emplace_back(emptyArray, qc);
  q.emplace_back(emptyArray, qv);

  // |  i  |  q  |    p    |   x  |
  // |  0  | lqr | prr_gsp |  qr  |
  // |  1  | lqi | pri_gsp |  qi  |
  // |  2  | lqs | prs_gsp |  qs  |
  // |  3  | lqg | prg_gsp |  qg  |
  // |  4  | lqc |         |  qc  |
  // |  5  | lqv |         |  qv  |

  size_t jmx = 0;
  size_t jmx_ = jmx;

  // The loop is intentionally i<nlev; since we are using an unsigned integer
  // data type, when i reaches 0, and you try to decrement further, (to -1), it
  // wraps to the maximum value representable by size_t.
  size_t oned_vec_index;
  for (size_t i = ke - 1; i < ke; --i) {
    for (size_t j = ivstart; j < ivend; j++) {
      oned_vec_index = i * ivend + j;
      if ((std::max({qc[oned_vec_index], qr[oned_vec_index], qs[oned_vec_index],
                     qi[oned_vec_index], qg[oned_vec_index]}) > qmin) or
          ((t[oned_vec_index] < tfrz_het2) and
           (qv[oned_vec_index] >
            qsat_ice_rho(t[oned_vec_index], rho[oned_vec_index])))) {
        jmx_ = jmx_ + 1;
        ind_k[jmx] = i;
        ind_i[jmx] = j;
        is_sig_present[jmx] = std::max({qs[oned_vec_index], qi[oned_vec_index],
                                        qg[oned_vec_index]}) > qmin;
        jmx = jmx_;
      }

      // qp_ind = {0, 1, 2, 3}
      // ix = 0, qp_ind[0] = 0
      if (i == (ke - 1)) {
        kmin[j][0] = ke + 1;
        prr_gsp[j] = 0.0;
        vt[j][0] = 0.0;
      }

      if (qr[oned_vec_index] > qmin) {
        kmin[j][0] = i;
      }

      // ix = 1, qp_ind[1] = 1
      if (i == (ke - 1)) {
        kmin[j][1] = ke + 1;
        pri_gsp[j] = 0.0;
        vt[j][1] = 0.0;
      }

      if (qi[oned_vec_index] > qmin) {
        kmin[j][1] = i;
      }

      // ix = 2, qp_ind[2] = 2
      if (i == (ke - 1)) {
        kmin[j][2] = ke + 1;
        prs_gsp[j] = 0.0;
        vt[j][2] = 0.0;
      }

      if (qs[oned_vec_index] > qmin) {
        kmin[j][2] = i;
      }

      // ix = 3, qp_ind[3] = 3
      if (i == (ke - 1)) {
        kmin[j][3] = ke + 1;
        prg_gsp[j] = 0.0;
        vt[j][3] = 0.0;
      }

      if (qg[oned_vec_index] > qmin) {
        kmin[j][3] = i;
      }
    }
  }

  size_t k, iv;
  real_t sx2x_sum;
  for (size_t j = 0; j < jmx_; j++) {
    k = ind_k[j];
    iv = ind_i[j];
    oned_vec_index = k * ivend + iv;

    dvsw =
        qv[oned_vec_index] - qsat_rho(t[oned_vec_index], rho[oned_vec_index]);
    qvsi = qsat_ice_rho(t[oned_vec_index], rho[oned_vec_index]);
    dvsi = qv[oned_vec_index] - qvsi;
    n_snow =
        snow_number(t[oned_vec_index], rho[oned_vec_index], qs[oned_vec_index]);
    l_snow = snow_lambda(rho[oned_vec_index], qs[oned_vec_index], n_snow);

    sx2x[lqc][lqr] = cloud_to_rain(t[oned_vec_index], qc[oned_vec_index],
                                   qr[oned_vec_index], qnc);
    sx2x[lqr][lqv] =
        rain_to_vapor(t[oned_vec_index], rho[oned_vec_index],
                      qc[oned_vec_index], qr[oned_vec_index], dvsw, dt);
    sx2x[lqc][lqi] = cloud_x_ice(t[oned_vec_index], qc[oned_vec_index],
                                 qi[oned_vec_index], dt);
    sx2x[lqi][lqc] = -std::fmin(sx2x[lqc][lqi], 0.0);
    sx2x[lqc][lqi] = std::fmax(sx2x[lqc][lqi], 0.0);
    sx2x[lqc][lqs] = cloud_to_snow(t[oned_vec_index], qc[oned_vec_index],
                                   qs[oned_vec_index], n_snow, l_snow);
    sx2x[lqc][lqg] = cloud_to_graupel(t[oned_vec_index], rho[oned_vec_index],
                                      qc[oned_vec_index], qg[oned_vec_index]);

    if (t[oned_vec_index] < tmelt) {
      n_ice = ice_number(t[oned_vec_index], rho[oned_vec_index]);
      m_ice = ice_mass(qi[oned_vec_index], n_ice);
      x_ice = ice_sticking(t[oned_vec_index]);

      if (is_sig_present[j]) {
        eta = deposition_factor(
            t[oned_vec_index],
            qvsi); // neglect cloud depth cor. from gcsp_graupel
        sx2x[lqv][lqi] = vapor_x_ice(qi[oned_vec_index], m_ice, eta, dvsi, dt);
        sx2x[lqi][lqv] = -std::fmin(sx2x[lqv][lqi], 0.0);
        sx2x[lqv][lqi] = std::fmax(sx2x[lqv][lqi], 0.0);
        ice_dep = std::fmin(sx2x[lqv][lqi], dvsi / dt);

        sx2x[lqi][lqs] =
            deposition_auto_conversion(qi[oned_vec_index], m_ice, ice_dep);
        sx2x[lqi][lqs] = sx2x[lqi][lqs] +
                         ice_to_snow(qi[oned_vec_index], n_snow, l_snow, x_ice);
        sx2x[lqi][lqg] =
            ice_to_graupel(rho[oned_vec_index], qr[oned_vec_index],
                           qg[oned_vec_index], qi[oned_vec_index], x_ice);
        sx2x[lqs][lqg] =
            snow_to_graupel(t[oned_vec_index], rho[oned_vec_index],
                            qc[oned_vec_index], qs[oned_vec_index]);
        sx2x[lqr][lqg] = rain_to_graupel(t[oned_vec_index], rho[oned_vec_index],
                                         qc[oned_vec_index], qr[oned_vec_index],
                                         qi[oned_vec_index], qs[oned_vec_index],
                                         m_ice, dvsw, dt);
      }
      sx2x[lqv][lqi] =
          sx2x[lqv][lqi] +
          ice_deposition_nucleation(t[oned_vec_index], qc[oned_vec_index],
                                    qi[oned_vec_index], n_ice, dvsi, dt);
    } else {
      sx2x[lqc][lqr] = sx2x[lqc][lqr] + sx2x[lqc][lqs] + sx2x[lqc][lqg];
      sx2x[lqc][lqs] = 0.0;
      sx2x[lqc][lqg] = 0.0;
      ice_dep = 0.0;
      eta = 0.0;
    }

    if (is_sig_present[j]) {
      dvsw0 = qv[oned_vec_index] - qsat_rho(tmelt, rho[oned_vec_index]);
      sx2x[lqv][lqs] =
          vapor_x_snow(t[oned_vec_index], p[oned_vec_index],
                       rho[oned_vec_index], qs[oned_vec_index], n_snow, l_snow,
                       eta, ice_dep, dvsw, dvsi, dvsw0, dt);
      sx2x[lqs][lqv] = -std::fmin(sx2x[lqv][lqs], 0.0);
      sx2x[lqv][lqs] = std::fmax(sx2x[lqv][lqs], 0.0);
      sx2x[lqv][lqg] = vapor_x_graupel(t[oned_vec_index], p[oned_vec_index],
                                       rho[oned_vec_index], qg[oned_vec_index],
                                       dvsw, dvsi, dvsw0, dt);
      sx2x[lqg][lqv] = -std::fmin(sx2x[lqv][lqg], 0.0);
      sx2x[lqv][lqg] = std::fmax(sx2x[lqv][lqg], 0.0);
      sx2x[lqs][lqr] =
          snow_to_rain(t[oned_vec_index], p[oned_vec_index],
                       rho[oned_vec_index], dvsw0, qs[oned_vec_index]);
      sx2x[lqg][lqr] =
          graupel_to_rain(t[oned_vec_index], p[oned_vec_index],
                          rho[oned_vec_index], dvsw0, qg[oned_vec_index]);
    }

    // qx_ind = {5, 4, 0, 2, 1, 3}
    // ix = 0, qx_ind[0] = 5
    sink[5] = 0.0;
    if ((is_sig_present[j]) or (5 == lqc) or (5 == lqv) or (5 == lqr)) {

      for (size_t i = 0; i < nx; i++) {
        sink[5] = sink[5] + sx2x[5][i];
      }
      stot = qv[oned_vec_index] / dt;

      if ((sink[5] > stot) && (qv[oned_vec_index] > qmin)) {
        real_t nextSink = 0.0;

        for (size_t i = 0; i < nx; i++) {
          sx2x[5][i] = sx2x[5][i] * stot / sink[5];
          nextSink = nextSink + sx2x[5][i];
        }
        sink[5] = nextSink;
      }
    }

    // ix = 1, qx_ind[1] = 4
    sink[4] = 0.0;
    if ((is_sig_present[j]) or (4 == lqc) or (4 == lqv) or (4 == lqr)) {

      for (size_t i = 0; i < nx; i++) {
        sink[4] = sink[4] + sx2x[4][i];
      }
      stot = qc[oned_vec_index] / dt;

      if ((sink[4] > stot) && (qc[oned_vec_index] > qmin)) {
        real_t nextSink = 0.0;

        for (size_t i = 0; i < nx; i++) {
          sx2x[4][i] = sx2x[4][i] * stot / sink[4];
          nextSink = nextSink + sx2x[4][i];
        }
        sink[4] = nextSink;
      }
    }

    // ix = 2, qx_ind[2] = 0
    sink[0] = 0.0;
    if ((is_sig_present[j]) or (0 == lqc) or (0 == lqv) or (0 == lqr)) {

      for (size_t i = 0; i < nx; i++) {
        sink[0] = sink[0] + sx2x[0][i];
      }
      stot = qr[oned_vec_index] / dt;

      if ((sink[0] > stot) && (qr[oned_vec_index] > qmin)) {
        real_t nextSink = 0.0;

        for (size_t i = 0; i < nx; i++) {
          sx2x[0][i] = sx2x[0][i] * stot / sink[0];
          nextSink = nextSink + sx2x[0][i];
        }
        sink[0] = nextSink;
      }
    }

    // ix = 3, qx_ind[3] = 2
    sink[2] = 0.0;
    if ((is_sig_present[j]) or (2 == lqc) or (2 == lqv) or (2 == lqr)) {

      for (size_t i = 0; i < nx; i++) {
        sink[2] = sink[2] + sx2x[2][i];
      }
      stot = qs[oned_vec_index] / dt;

      if ((sink[2] > stot) && (qs[oned_vec_index] > qmin)) {
        real_t nextSink = 0.0;

        for (size_t i = 0; i < nx; i++) {
          sx2x[2][i] = sx2x[2][i] * stot / sink[2];
          nextSink = nextSink + sx2x[2][i];
        }
        sink[2] = nextSink;
      }
    }

    // ix = 4, qx_ind[4] = 1
    sink[1] = 0.0;
    if ((is_sig_present[j]) or (1 == lqc) or (1 == lqv) or (1 == lqr)) {

      for (size_t i = 0; i < nx; i++) {
        sink[1] = sink[1] + sx2x[1][i];
      }
      stot = qi[oned_vec_index] / dt;

      if ((sink[1] > stot) && (qi[oned_vec_index] > qmin)) {
        real_t nextSink = 0.0;

        for (size_t i = 0; i < nx; i++) {
          sx2x[1][i] = sx2x[1][i] * stot / sink[1];
          nextSink = nextSink + sx2x[1][i];
        }
        sink[1] = nextSink;
      }
    }

    // ix = 5, qx_ind[5] = 3
    sink[3] = 0.0;
    if ((is_sig_present[j]) or (3 == lqc) or (3 == lqv) or (3 == lqr)) {

      for (size_t i = 0; i < nx; i++) {
        sink[3] = sink[3] + sx2x[3][i];
      }
      stot = qg[oned_vec_index] / dt;

      if ((sink[3] > stot) && (qg[oned_vec_index] > qmin)) {
        real_t nextSink = 0.0;

        for (size_t i = 0; i < nx; i++) {
          sx2x[3][i] = sx2x[3][i] * stot / sink[3];
          nextSink = nextSink + sx2x[3][i];
        }
        sink[3] = nextSink;
      }
    }

    // qx_ind = {5, 4, 0, 2, 1, 3}
    // ix = 0, qx_ind[0] = 5
    sx2x_sum = 0;
    for (size_t i = 0; i < nx; i++) {
      sx2x_sum = sx2x_sum + sx2x[i][5];
    }
    dqdt[5] = sx2x_sum - sink[5];
    qv[oned_vec_index] = std::fmax(0.0, qv[oned_vec_index] + dqdt[5] * dt);

    // ix = 1, qx_ind[1] = 4
    sx2x_sum = 0;
    for (size_t i = 0; i < nx; i++) {
      sx2x_sum = sx2x_sum + sx2x[i][4];
    }
    dqdt[4] = sx2x_sum - sink[4];
    qc[oned_vec_index] = std::fmax(0.0, qc[oned_vec_index] + dqdt[4] * dt);

    // ix = 2, qx_ind[2] = 0
    sx2x_sum = 0;
    for (size_t i = 0; i < nx; i++) {
      sx2x_sum = sx2x_sum + sx2x[i][0];
    }
    dqdt[0] = sx2x_sum - sink[0];
    qr[oned_vec_index] = std::fmax(0.0, qr[oned_vec_index] + dqdt[0] * dt);

    // ix = 3, qx_ind[3] = 2
    sx2x_sum = 0;
    for (size_t i = 0; i < nx; i++) {
      sx2x_sum = sx2x_sum + sx2x[i][2];
    }
    dqdt[2] = sx2x_sum - sink[2];
    qs[oned_vec_index] = std::fmax(0.0, qs[oned_vec_index] + dqdt[2] * dt);

    // ix = 4, qx_ind[4] = 1
    sx2x_sum = 0;
    for (size_t i = 0; i < nx; i++) {
      sx2x_sum = sx2x_sum + sx2x[i][1];
    }
    dqdt[1] = sx2x_sum - sink[1];
    qi[oned_vec_index] = std::fmax(0.0, qi[oned_vec_index] + dqdt[1] * dt);

    // ix = 5, qx_ind[5] = 3
    sx2x_sum = 0;
    for (size_t i = 0; i < nx; i++) {
      sx2x_sum = sx2x_sum + sx2x[i][3];
    }
    dqdt[3] = sx2x_sum - sink[3];
    qg[oned_vec_index] = std::fmax(0.0, qg[oned_vec_index] + dqdt[3] * dt);

    qice = qs[oned_vec_index] + qi[oned_vec_index] + qg[oned_vec_index];
    qliq = qc[oned_vec_index] + qr[oned_vec_index];
    qtot = qv[oned_vec_index] + qice + qliq;
    cv = cvd + (cvv - cvd) * qtot + (clw - cvv) * qliq +
         (ci - cvv) * qice; // qtot? or qv?
    t[oned_vec_index] =
        t[oned_vec_index] +
        dt *
            ((dqdt[lqc] + dqdt[lqr]) * (lvc - (clw - cvv) * t[oned_vec_index]) +
             (dqdt[lqi] + dqdt[lqs] + dqdt[lqg]) *
                 (lsc - (ci - cvv) * t[oned_vec_index])) /
            cv;

    // reset all values of sx2x to zero
    for (auto &v : sx2x) {
      std::fill(v.begin(), v.end(), 0);
    }
  }

  size_t kp1;
  size_t k_end = (lrain) ? ke : kstart - 1;
  for (size_t k = kstart; k < k_end; k++) {
    for (size_t iv = ivstart; iv < ivend; iv++) {
      oned_vec_index = k * ivend + iv;
      if (k == kstart) {
        eflx[iv] = 0.0;
      }

      kp1 = std::min(ke - 1, k + 1);
      if (k >= *std::min_element(kmin[iv].begin(), kmin[iv].end())) {
        qliq = qc[oned_vec_index] + qr[oned_vec_index];
        qice = qs[oned_vec_index] + qi[oned_vec_index] + qg[oned_vec_index];

        e_int = internal_energy(t[oned_vec_index], qv[oned_vec_index], qliq,
                                qice, rho[oned_vec_index], dz[oned_vec_index]) +
                eflx[iv];
        zeta = dt / (2.0 * dz[oned_vec_index]);
        xrho = std::sqrt(rho_00 / rho[oned_vec_index]);

        // qp_ind = {0, 1, 2, 3}
        // ix = 0, qp_ind[0] = 0
        if (k >= kmin[iv][0]) {
          vc = vel_scale_factor(0, xrho, rho[oned_vec_index], t[oned_vec_index],
                                qr[oned_vec_index]);
          precip(params[0], update, zeta, vc, prr_gsp[iv], vt[iv][0],
                 qr[oned_vec_index], qr[kp1 * ivend + iv], rho[oned_vec_index]);
          qr[oned_vec_index] = update[0];
          prr_gsp[iv] = update[1];
          vt[iv][0] = update[2];
        }

        // ix = 1, qp_ind[1] = 1
        if (k >= kmin[iv][1]) {
          vc = vel_scale_factor(1, xrho, rho[oned_vec_index], t[oned_vec_index],
                                qi[oned_vec_index]);
          precip(params[1], update, zeta, vc, pri_gsp[iv], vt[iv][1],
                 qi[oned_vec_index], qi[kp1 * ivend + iv], rho[oned_vec_index]);
          qi[oned_vec_index] = update[0];
          pri_gsp[iv] = update[1];
          vt[iv][1] = update[2];
        }

        // ix = 2, qp_ind[2] = 2
        if (k >= kmin[iv][2]) {
          vc = vel_scale_factor(2, xrho, rho[oned_vec_index], t[oned_vec_index],
                                qs[oned_vec_index]);
          precip(params[2], update, zeta, vc, prs_gsp[iv], vt[iv][2],
                 qs[oned_vec_index], qs[kp1 * ivend + iv], rho[oned_vec_index]);
          qs[oned_vec_index] = update[0];
          prs_gsp[iv] = update[1];
          vt[iv][2] = update[2];
        }

        // ix = 3, qp_ind[3] = 3
        if (k >= kmin[iv][3]) {
          vc = vel_scale_factor(3, xrho, rho[oned_vec_index], t[oned_vec_index],
                                qg[oned_vec_index]);
          precip(params[3], update, zeta, vc, prg_gsp[iv], vt[iv][3],
                 qg[oned_vec_index], qg[kp1 * ivend + iv], rho[oned_vec_index]);
          qg[oned_vec_index] = update[0];
          prg_gsp[iv] = update[1];
          vt[iv][3] = update[2];
        }

        pflx[oned_vec_index] = prs_gsp[iv] + pri_gsp[iv] + prg_gsp[iv];
        eflx[iv] =
            dt * (prr_gsp[iv] * (clw * t[oned_vec_index] -
                                 cvd * t[kp1 * ivend + iv] - lvc) +
                  pflx[oned_vec_index] * (ci * t[oned_vec_index] -
                                          cvd * t[kp1 * ivend + iv] - lsc));
        pflx[oned_vec_index] = pflx[oned_vec_index] + prr_gsp[iv];
        qliq = qc[oned_vec_index] + qr[oned_vec_index];
        qice = qs[oned_vec_index] + qi[oned_vec_index] + qg[oned_vec_index];
        e_int = e_int - eflx[iv];
        t[oned_vec_index] =
            T_from_internal_energy(e_int, qv[oned_vec_index], qliq, qice,
                                   rho[oned_vec_index], dz[oned_vec_index]);
      }
    }
  }
}
