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

#if defined(MU_ENABLE_GPU)
#include <cuda_runtime.h>
#endif

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
  std::unique_ptr<real_t[]> z, t, p, rho, qv, qc, qi, qr, qs, qg;
  // Pre-calculated parameters
  std::unique_ptr<real_t[]> dz;
#else
  // Parameters from the input file
  real_t *z, *t, *p, *rho, *qv, *qc, *qi, *qr, *qs, *qg;
  // Pre-calculated parameters
  real_t *dz;
#endif

  // start-end indices
  size_t kend, kbeg, ivend, ivbeg, nvec;

  const string input_file = file;
  io_muphys::read_fields(input_file, itime, ncells, nlev, z, t, p, rho, qv, qc,
                         qi, qr, qs, qg);
#if defined(MU_ENABLE_SEQ)
  utils_muphys::calc_dz(z, dz, ncells, nlev);
#elif defined(MU_ENABLE_OMP)
  utils_muphys::calc_dz(z.get(), dz, ncells, nlev);
#else
  utils_muphys::calc_dz(z, dz, ncells, nlev);
#endif

#if defined(MU_ENABLE_SEQ)
  prr_gsp.resize(ncells, 0.0);
  pri_gsp.resize(ncells, 0.0);
  prs_gsp.resize(ncells, 0.0);
  prg_gsp.resize(ncells, 0.0);
  pflx.resize(ncells * nlev, 0.0);
#endif

  kbeg = 0;
  kend = nlev;
  ivbeg = 0;
  ivend = ncells;
  nvec = ncells;
  qnc_1 = qnc;

  auto start_time = std::chrono::steady_clock::now();

#if defined(MU_ENABLE_SEQ)
  graupel(nvec, kend, ivbeg, ivend, kbeg, dt, dz, t, rho, p, qv, qc, qi, qr, qs,
          qg, qnc_1, prr_gsp, pri_gsp, prs_gsp, prg_gsp, pflx);
#elif defined(MU_ENABLE_OMP)
  graupel(nvec, kend, ivbeg, ivend, kbeg, dt, dz.get(), t.get(), rho.get(),
          p.get(), qv.get(), qc.get(), qi.get(), qr.get(), qs.get(), qg.get(),
          qnc_1);
#else
  graupel(nvec, kend, ivbeg, ivend, kbeg, dt, dz, t, rho, p, qv, qc, qi, qr, qs,
          qg, qnc_1);
#endif

  auto end_time = std::chrono::steady_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      end_time - start_time);

#if defined(MU_ENABLE_SEQ)
  io_muphys::write_fields(output_file, ncells, nlev, t, qv, qc, qi, qr, qs, qg);
#elif defined(MU_ENABLE_OMP)
  io_muphys::write_fields(output_file, ncells, nlev, t.get(), qv.get(),
                          qc.get(), qi.get(), qr.get(), qs.get(), qg.get());
#else
  io_muphys::write_fields(output_file, ncells, nlev, t, qv, qc, qi, qr, qs, qg);
#endif

  std::cout << "time taken : " << duration.count() << " milliseconds"
            << std::endl;

  if (std::getenv("MU_LOG_TIME")) {
    io_muphys::log_time(duration.count());
  }
#ifdef MU_ENABLE_GPU
  cudaFreeHost(dz);
  cudaFreeHost(t);
  cudaFreeHost(rho);
  cudaFreeHost(p);
  cudaFreeHost(qv);
  cudaFreeHost(qc);
  cudaFreeHost(qi);
  cudaFreeHost(qr);
  cudaFreeHost(qs);
  cudaFreeHost(qg);
#endif

  return 0;
}
