#!/bin/bash
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

export LD_LIBRARY_PATH=/sw/spack-levante/nvhpc-23.7-xasprs/Linux_x86_64/23.7/profilers/Nsight_Systems/host-linux-x64/:/sw/spack-levante/nvhpc-23.7-xasprs/Linux_x86_64/23.7/cuda/lib64

ncu --set full --metric all --export ncu-report build/bin/graupel $1
