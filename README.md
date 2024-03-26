# muphys-cpp

## Description
C++ prototypes of [muphys](https://gitlab.dkrz.de/icon-libraries/muphys) project using heterogeneous libraries and C extensions.

## Documentation
Online documentation is generated automatically using `doxygen` .

## Installation

### Dependencies
* [NetCDF for CXX](https://github.com/Unidata/netcdf-cxx4)
  * for Levante: `spack load netcdf-cxx4@4.3.1`

Other dependency like [googletest](https://github.com/google/googletest) is built in-tree from github archives. 

### Available compile options 
* _Implementation_ - The sequential implementation is selected by default. The user can choose of the following options:
  * MU_IMPL=seq - C++ serial implementation
 * _Precision_ (default is `double`)
  * MU_ENABLE_SINGLE - to switch to `float` 
* _Unit-test_ - compile tests together with the main executable (default is `true`)
  * MU_ENABLE_TESTS

### Compile the project (with default flags and Seq frontend)

`cmake -DMU_IMPL=seq -B build -S .`

`cmake --build build`

## Usage

`./build/bin/graupel input_file.nc`

### Automated tests

- Run tests manually:
`cd build && ctest` 

## License

muphys-cpp is available under a BSD 3-clause license. See [LICENSES/](./LICENSES) for license information and [AUTHORS.TXT](./AUTHORS.TXT) for a list of authors.

## Coding Challenge

- Task: Creat an optimized parallel implementation of muphys
    - Requirements:
        - 4 Versions: Correctness/Performance x CPU/GPU
        - Directive Programming: Use OpenMP (Offload to GPU) !No CUDA Kernels            
    - Parallelisation Task:
        - Parallelize the sequential implementation
        - Can rewrite the algorithm entirely
    - Optimization Task:
        - Look reordering/collapse/rewrte
        - Minimize critical regions
        - Data structure to optimize cache utilization
        - Optimal CPU-GPU memory transfer

- Layout:
    - code: 
        - /core: physics interactions
        - /io: read/write/ I/O files
        - /implementations/sequential
    - external dependencies:
        - /extern
    - validation:
        - /test: automated unit tests
        - /tasks: input files
        - /reference_results: output fiels
        - /scripts: Levante batch files

- Validation:
    - CPU-CPU: `cdo diffv` bit-identical
    - CPU-GPU: `cdo sub` within floating point tolerance

- Submission:
    - Path to the implementation
    - Script to build the 4 versions
    - Slurm logs to confirm the results
    - Summary list of optimizations performed
    - Plots to confirm the performance results
    - (Optional) Profiler analysis output & interpretation
    - (Optional) Experience report for using OpenMP

- Grading:
    - Correctness(50%): results are correct on both CPU and GPU 
    - Performance(50%): the code is faster than the sequential version
        - How to determine? 
        - Will the teams be ranked?
            - No, but the fastest may get a bonus
    - Bonus(20%): readability, portability, extreme performance, optional documentation ...        
