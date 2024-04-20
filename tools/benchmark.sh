#!/bin/bash

set -e

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

if [ -z "$2" ] || [ $2 != "--double" ]; then
  MU_ENABLE_SINGLE=ON
  build_precision="single"
else
  MU_ENABLE_SINGLE=OFF
  build_precision="double"
fi

BUILD_PREFIX=build_benchmark

# rm -rf $BUILD_PREFIX 
mkdir -p $BUILD_PREFIX

cmake  -DCMAKE_CXX_COMPILER=$BUILD_CXX_COMPILER \
       -DCMAKE_CXX_FLAGS="-Ofast" \
       -DMU_IMPL=$MU_IMPL \
       -DMU_ENABLE_SINGLE=$MU_ENABLE_SINGLE \
       -B $BUILD_PREFIX/build_$build_precision

cmake --build $BUILD_PREFIX/build_$build_precision --parallel

$BUILD_PREFIX/build_$build_precision/bin/graupel tasks/1500k.nc

