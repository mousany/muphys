# Summary list of optimization performed

0. Baseline
40107af545e2b5d47a334d2557968c3a46933001
1. Parallelize the second and third for loop
7d3006105f5e9cf5e4b2ec9663c5c87685f91955
2. Merge the first and second for loop to remove critical section
213bfe7eb4b711fe2fa76664bdacacf36c7aa49d
3. Eliminate std::vector, std::min and std::max and remove unnecessary temporary variables
cc6bcc96a259b8d1a9a83a1da616f2018d81db94
    - Remove std::min and std::max to help the compiler optimize out redundant loads
    - Remove unnecessary temporary variables which are shared between threads to avoid false sharing
    - Remove std::vector that is used temporarily in each iteration to avoid memory allocation and deallocation
4. Cache blocking
035dd1c332900ccefd7fc6a862ecd25b18f82997
5. Use -Ofast for faster pow function
36a0300926eb251184c22df4f9d6d23fa28b7882
5.5 Merge the remaining two for loops
d15f88b1d54b923051e3b7b70a3353a13747cbd4
6. Vectorize short loops and align with cache line
4f26e3fcc204cb453aa92cf746856ee9c56e810b
7. Move to GPU
b9aec28649097344fc825444624581ebe8bc5d98
8. Use pinned memory for faster data transfer
3587520ada8fea1c5d81148038f9a85554afcefd
9. Pipeline the memory transfer
549f67bc492b5a30d6dadfbe580a82959d61d577
