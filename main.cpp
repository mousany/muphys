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
#include <chrono>
#include <cstdlib>
#include <iostream>

#include "core/common/graupel.hpp"
#include "core/common/types.hpp"
#include "core/common/utils.hpp"
#include "io/io.hpp"
#include <chrono>

int main(int argc, char *argv[]) {
  // Parameters from the command line
  string file;
  string output_file = "output.nc";
  size_t itime;
  real_t dt, qnc, qnc_1;
  io_muphys::parse_args(file, itime, dt, qnc, argc, argv);

  size_t ncells, nlev;
#if defined(MU_ENABLE_SEQ)
  // Parameters from the input file
  array_1d_t<real_t> z, t, p, rho, qv, qc, qi, qr, qs, qg;
  // Pre-calculated parameters
  array_1d_t<real_t> dz;
  // Extra fields required to call graupel
  array_1d_t<real_t> pflx, prr_gsp, pri_gsp, prs_gsp, prg_gsp;
#elif defined(MU_ENABLE_OMP)
  // Parameters from the input file
  buffer_1d_t<real_t> z, t, p, rho, qv, qc, qi, qr, qs, qg;
  // Pre-calculated parameters
  buffer_1d_t<real_t> dz;
  // Extra fields required to call graupel
  buffer_1d_t<real_t> pflx, prr_gsp, pri_gsp, prs_gsp, prg_gsp;
#else
#error No implementation was selected. Options are: seq, omp
#endif

  // start-end indices
  size_t kend, kbeg, ivend, ivbeg, nvec;

  const string input_file = file;
  io_muphys::read_fields(input_file, itime, ncells, nlev, z, t, p, rho, qv, qc,
                         qi, qr, qs, qg);
  utils_muphys::calc_dz(z, dz, ncells, nlev);

#if defined(MU_ENABLE_SEQ)
  prr_gsp.resize(ncells, 0.0);
  pri_gsp.resize(ncells, 0.0);
  prs_gsp.resize(ncells, 0.0);
  prg_gsp.resize(ncells, 0.0);
  pflx.resize(ncells * nlev, 0.0);
#elif defined(MU_ENABLE_OMP)
  prr_gsp = new real_t[ncells]();
  pri_gsp = new real_t[ncells]();
  prs_gsp = new real_t[ncells]();
  prg_gsp = new real_t[ncells]();
  pflx = new real_t[ncells * nlev]();
#else
#error No implementation was selected. Options are: seq, omp
#endif

  kbeg = 0;
  kend = nlev;
  ivbeg = 0;
  ivend = ncells;
  nvec = ncells;
  qnc_1 = qnc;

  auto start_time = std::chrono::steady_clock::now();

  graupel(nvec, kend, ivbeg, ivend, kbeg, dt, dz, t, rho, p, qv, qc, qi, qr, qs,
          qg, qnc_1, prr_gsp, pri_gsp, prs_gsp, prg_gsp, pflx);

  auto end_time = std::chrono::steady_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      end_time - start_time);

  io_muphys::write_fields(output_file, ncells, nlev, t, qv, qc, qi, qr, qs, qg);
  std::cout << "time taken : " << duration.count() << " milliseconds"
            << std::endl;

  if (std::getenv("MU_LOG_TIME")) {
    io_muphys::log_time(duration.count());
  }

#if defined(MU_ENABLE_SEQ)

#elif defined(MU_ENABLE_OMP)
  delete[] prr_gsp;
  delete[] pri_gsp;
  delete[] prs_gsp;
  delete[] prg_gsp;
  delete[] pflx;

  delete[] dz;

  delete[] z;
  delete[] t;
  delete[] p;
  delete[] rho;
  delete[] qv;
  delete[] qc;
  delete[] qi;
  delete[] qr;
  delete[] qs;
  delete[] qg;
#else
#error No implementation was selected. Options are: seq, omp
#endif

  return 0;
}
