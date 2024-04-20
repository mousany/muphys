#!/bin/bash

rm -rf build_single build_double

cmake -DCMAKE_CXX_COMPILER=nvc++ \
      -DCMAKE_CXX_FLAGS="-Ofast" \
      -DMU_IMPL=gpu \
      -DMU_ENABLE_SINGLE=ON \
      -B build_single -S . \
      && cmake --build build_single

cmake -DCMAKE_CXX_COMPILER=nvc++ \
        -DCMAKE_CXX_FLAGS="-Ofast" \
        -DMU_IMPL=gpu \
        -DMU_ENABLE_SINGLE=OFF \
        -B build_double -S . \
        && cmake --build build_double
