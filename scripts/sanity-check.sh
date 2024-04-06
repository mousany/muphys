#!/bin/bash

if [ -z "$1" ]; then
  MU_IMPL="seq"
else
  MU_IMPL=$1
fi

# if MU_IMPL is gpu, use nvc++
if [ $MU_IMPL == "gpu" ]; then
  BUILD_CXX_COMPILER=nvc++
else
  BUILD_CXX_COMPILER=g++
fi

BUILD_PREFIX=build_sanity

# rm -rf $BUILD_PREFIX 
mkdir -p $BUILD_PREFIX

cmake  -DCMAKE_CXX_COMPILER=$BUILD_CXX_COMPILER \
       -DCMAKE_CXX_FLAGS="-O0" \
       -DMU_IMPL=$MU_IMPL \
       -DMU_ENABLE_SINGLE=ON \
       -B $BUILD_PREFIX/build_single

cmake --build $BUILD_PREFIX/build_single

$BUILD_PREFIX/build_single/bin/graupel tasks/dbg.nc
cdo diffn output.nc reference_results/dbg_single.nc

$BUILD_PREFIX/build_single/bin/graupel tasks/input.nc
cdo diffn output.nc reference_results/sequential_single_output.nc

$BUILD_PREFIX/build_single/bin/graupel tasks/20k.nc
cdo diffn output.nc reference_results/sequential_single_20k.nc


# cmake  -DCMAKE_CXX_COMPILER=$BUILD_CXX_COMPILER \
#        -DCMAKE_CXX_FLAGS="-O0" \
#        -DMU_IMPL=$MU_IMPL \
#        -DMU_ENABLE_SINGLE=OFF \
#        -B $BUILD_PREFIX/build_double

# cmake --build $BUILD_PREFIX/build_double

# $BUILD_PREFIX/build_double/bin/graupel tasks/dbg.nc
# cdo diffn output.nc reference_results/dbg_double.nc

# $BUILD_PREFIX/build_double/bin/graupel tasks/input.nc
# cdo diffn output.nc reference_results/sequential_double_output.nc

# $BUILD_PREFIX/build_double/bin/graupel tasks/20k.nc
# cdo diffn output.nc reference_results/sequential_double_20k.nc
