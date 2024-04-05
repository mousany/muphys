#!/bin/bash

if [ -z "$1" ]; then
  MU_IMPL="seq"
else
  MU_IMPL=$1
fi

BUILD_PREFIX=build_benchmark

# rm -rf $BUILD_PREFIX 
mkdir -p $BUILD_PREFIX

cmake  -DCMAKE_CXX_COMPILER=g++ \
       -DCMAKE_CXX_FLAGS="-Ofast" \
       -DMU_IMPL=$MU_IMPL \
       -DMU_ENABLE_SINGLE=ON \
       -B $BUILD_PREFIX/build_single

cmake --build $BUILD_PREFIX/build_single

$BUILD_PREFIX/build_single/bin/graupel tasks/1500k.nc

