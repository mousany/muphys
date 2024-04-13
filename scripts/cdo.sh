#!/bin/bash
module load cdo
LD_LIBRARY_PATH=/sw/spack-levante/cdo-2.2.2-4z4icb/lib cdo infon -sub output.nc $1
