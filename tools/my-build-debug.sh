#!/bin/bash
#SBATCH --account=ka1273
#SBATCH --job-name=scc
#SBATCH --partition=gpu
#SBATCH --gpus=1
#SBATCH --output=%x.%j.log
#SBATCH --exclusive
#SBATCH --time=00:10:00

ulimit -s unlimited
ulimit -c 0

. scripts/levante-setup.sh nvidia gpu
CC=nvc CXX=nvc++ cmake -DMU_ENABLE_SINGLE=on -DCMAKE_BUILD_TYPE=Debug -DMU_ENABLE_TESTS=false -DMU_IMPL=seq -DMU_GPU_VERSION=cc80 -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -B build-debug -S .
cmake --build build-debug

