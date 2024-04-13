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
#include "../core/common/types.hpp"
#include <fstream>
#include <iostream>
#include <netcdf>
using namespace netCDF;
using namespace netCDF::exceptions;

namespace io_muphys {

#ifdef __SINGLE_PRECISION
using NCreal_t = NcFloat;
#else
using NCreal_t = NcDouble;
#endif

void parse_args(std::string &file, size_t &itime, real_t &dt, real_t &qnc,
                int argc, char **argv);

#if defined(MU_ENABLE_SEQ)
void input_vector(NcFile &datafile, array_1d_t<real_t> &v,
                  const std::string input, size_t &ncells, size_t &nlev,
                  size_t itime);
void input_vector(NcFile &datafile, array_1d_t<real_t> &v,
                  const std::string input, size_t &ncells, size_t &nlev);
#elif defined(MU_ENABLE_OMP)
void input_vector(NcFile &datafile, std::unique_ptr<real_t[]> &v,
                  const std::string &input, size_t ncells, size_t nlev,
                  size_t itime);
void input_vector(NcFile &datafile, std::unique_ptr<real_t[]> &v,
                  const std::string &input, size_t ncells, size_t nlev);
#else
void input_vector(NcFile &datafile, real_t *&v, const std::string &input,
                  size_t ncells, size_t nlev, size_t itime);
void input_vector(NcFile &datafile, real_t *&v, const std::string &input,
                  size_t ncells, size_t nlev);

#endif

#if defined(MU_ENABLE_SEQ)
void output_vector(NcFile &datafile, array_1d_t<NcDim> &dims,
                   const std::string output, array_1d_t<real_t> &v,
                   size_t &ncells, size_t &nlev, int &deflate_level);
void output_vector(NcFile &datafile, array_1d_t<NcDim> &dims,
                   const std::string output, std::map<std::string, NcVarAtt>,
                   array_1d_t<real_t> &v, size_t &ncells, size_t &nlev,
                   int &deflate_level);
#else
void output_vector(NcFile &datafile, std::vector<NcDim> &dims,
                   const std::string &output, real_t *v, size_t ncells,
                   size_t nlev, int deflate_level);

#endif

#if defined(MU_ENABLE_SEQ)
void read_fields(const std::string input_file, size_t &itime, size_t &ncells,
                 size_t &nlev, array_1d_t<real_t> &z, array_1d_t<real_t> &t,
                 array_1d_t<real_t> &p, array_1d_t<real_t> &rho,
                 array_1d_t<real_t> &qv, array_1d_t<real_t> &qc,
                 array_1d_t<real_t> &qi, array_1d_t<real_t> &qr,
                 array_1d_t<real_t> &qs, array_1d_t<real_t> &qg);
#elif defined(MU_ENABLE_OMP)
void read_fields(const std::string &input_file, size_t &itime, size_t &ncells,
                 size_t &nlev, std::unique_ptr<real_t[]> &z,
                 std::unique_ptr<real_t[]> &t, std::unique_ptr<real_t[]> &p,
                 std::unique_ptr<real_t[]> &rho, std::unique_ptr<real_t[]> &qv,
                 std::unique_ptr<real_t[]> &qc, std::unique_ptr<real_t[]> &qi,
                 std::unique_ptr<real_t[]> &qr, std::unique_ptr<real_t[]> &qs,
                 std::unique_ptr<real_t[]> &qg);
#else
void read_fields(const std::string &input_file, size_t &itime, size_t &ncells,
                 size_t &nlev, real_t *&z, real_t *&t, real_t *&p, real_t *&rho,
                 real_t *&qv, real_t *&qc, real_t *&qi, real_t *&qr,
                 real_t *&qs, real_t *&qg);
#endif

#if defined(MU_ENABLE_SEQ)
void write_fields(const string output_file, size_t &ncells, size_t &nlev,
                  array_1d_t<real_t> &t, array_1d_t<real_t> &qv,
                  array_1d_t<real_t> &qc, array_1d_t<real_t> &qi,
                  array_1d_t<real_t> &qr, array_1d_t<real_t> &qs,
                  array_1d_t<real_t> &qg);
void write_fields(const string output_file, const string input_file,
                  size_t &ncells, size_t &nlev, array_1d_t<real_t> &t,
                  array_1d_t<real_t> &qv, array_1d_t<real_t> &qc,
                  array_1d_t<real_t> &qi, array_1d_t<real_t> &qr,
                  array_1d_t<real_t> &qs, array_1d_t<real_t> &qg);
#else
void write_fields(const string &output_file, size_t ncells, size_t nlev,
                  real_t *t, real_t *qv, real_t *qc, real_t *qi, real_t *qr,
                  real_t *qs, real_t *qg);

#endif

void log_time(int64_t value);
} // namespace io_muphys
