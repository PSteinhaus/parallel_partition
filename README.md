# p_partition

This project holds the source code for a parallelized `partition` function, using OpenMP, as well as `quicksort` and `nth_element`.

The code can be found in the header-only library file [`library.h`](my_lib/include/library.h).

[`catch_tests.cpp`](my_lib/catch_tests.cpp) and [`main.cpp`](main.cpp) hold unit tests and a simple benchmark respectively. As the benchmark compares the performance of quicksort with `__gnu_parallel::sort` it depends upon using the GCC, or an alternative way to access this function.
