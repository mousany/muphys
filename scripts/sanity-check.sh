#!/bin/bash

function diffn-gpu {
  threshold=1e-11

  awk -v thr="$threshold" '{
      all_rows=all_rows $0 "\n"; # Store all rows
      if (NF < 10) next; # Skip non-data lines

      min = $8; mean = $9; max = $10;

      if ((min < -thr || min > thr) || (mean < -thr || mean > thr) || (max < -thr || max > thr)) {
          count_differs++;
      }
      total++;
  }

  END {
      if (count_differs > 0) {
          print all_rows;
          printf "%d out of %d records differ\n", count_differs, total;
      }
  }'
}

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

function check {
  if [ $MU_IMPL == "gpu" ]; then
    cdo infon -sub output.nc $1 | diffn-gpu
  else
    cdo diffn output.nc $1
  fi
}

BUILD_PREFIX=build_sanity

# rm -rf $BUILD_PREFIX 
mkdir -p $BUILD_PREFIX

cmake  -DCMAKE_CXX_COMPILER=$BUILD_CXX_COMPILER \
       -DCMAKE_CXX_FLAGS="-O0" \
       -DMU_IMPL=$MU_IMPL \
       -DMU_ENABLE_SINGLE=ON \
       -B $BUILD_PREFIX/build_single && \
cmake --build $BUILD_PREFIX/build_single

if [ $? -ne 0 ]; then
  exit 1
fi

$BUILD_PREFIX/build_single/bin/graupel tasks/dbg.nc && \
check reference_results/dbg_single.nc

$BUILD_PREFIX/build_single/bin/graupel tasks/input.nc && \
check reference_results/sequential_single_output.nc

$BUILD_PREFIX/build_single/bin/graupel tasks/20k.nc && \
check reference_results/sequential_single_20k.nc


# cmake  -DCMAKE_CXX_COMPILER=$BUILD_CXX_COMPILER \
#        -DCMAKE_CXX_FLAGS="-O0" \
#        -DMU_IMPL=$MU_IMPL \
#        -DMU_ENABLE_SINGLE=OFF \
#        -B $BUILD_PREFIX/build_double && \
# cmake --build $BUILD_PREFIX/build_double

# if [ $? -ne 0 ]; then
#   exit 1
# fi

# $BUILD_PREFIX/build_double/bin/graupel tasks/dbg.nc && \
# check reference_results/dbg_double.nc

# $BUILD_PREFIX/build_double/bin/graupel tasks/input.nc && \
# check reference_results/sequential_double_output.nc

# $BUILD_PREFIX/build_double/bin/graupel tasks/20k.nc && \
# check reference_results/sequential_double_20k.nc
