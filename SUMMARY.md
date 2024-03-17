# Summary 

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