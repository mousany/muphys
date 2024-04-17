#!/bin/bash

set -e

function diffn-gpu {
  threshold=10e-11

  awk -v thr="$threshold" '{
      all_rows=all_rows $0 "\n"; # Store all rows
      if (NF < 10) next; # Skip non-data lines

      min = $8+0; mean = $9+0; max = $10+0; # Convert to numerical values

      if ((min < -thr || min > thr) || (mean < -thr || mean > thr) || (max < -thr || max > thr)) {
          count_differs++;
      }
      total++;
  }

  END {
      if (count_differs > 0) {
          print all_rows;
          printf "%d out of %d records differ\n", count_differs, total;
          exit 1;
      }
  }'
}

if [ -z "$1" ]; then
  MU_IMPL="seq"
else
  MU_IMPL=$1
fi

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
       -DMU_ENABLE_SINGLE=$MU_ENABLE_SINGLE \
       -B $BUILD_PREFIX/build_$build_precision

cmake --build $BUILD_PREFIX/build_$build_precision --parallel

$BUILD_PREFIX/build_$build_precision/bin/graupel tasks/dbg.nc
check reference_results/dbg_$build_precision.nc

$BUILD_PREFIX/build_$build_precision/bin/graupel tasks/input.nc
check reference_results/sequential_${build_precision}_output.nc

$BUILD_PREFIX/build_$build_precision/bin/graupel tasks/20k.nc
check reference_results/sequential_${build_precision}_20k.nc
