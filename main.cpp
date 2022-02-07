#include <iostream>
#include "library.h"
#include <vector>
#include <parallel/algorithm>
#include <omp.h>

int main(int argc, char *argv[])
{
    std::vector<float> v = {9,8,7,6,5,4,3,2,1,0};
    int blocklength = 5;
    float pivot = 4.5;


    // test parallel partition (phase one)  (block size is currently 2)
    int leftNeutralized, rightNeutralized;
    p_partition::partition(v.begin(),
                           v.end(),
                           [pivot](const auto &em) { return em < pivot; },
                           1);
    for (int i=0; i<v.size(); i++)
        std::cout << v[i] << " ";
    std::cout << std::endl;


    std::cout << leftNeutralized << " " << rightNeutralized << std::endl;

    //p_partition::testFunc();


    // BENCHMARK:
    //srand (time(nullptr)); // randomly seed rand

    const int ITERATIONS = 100;
    const int VECTOR_SIZE = 10000;
    auto comparator = [](const auto& em1, const auto& em2){ return em1 < em2; };

    double start_gnu = omp_get_wtime();
    for(int i=0; i<ITERATIONS; ++i) {
        std::vector<int> vec(VECTOR_SIZE);
        std::generate(vec.begin(), vec.end(), rand);
        __gnu_parallel::sort(vec.begin(), vec.end(), comparator);
    }
    double end_gnu = omp_get_wtime();
    double total_gnu = end_gnu - start_gnu;

    std::cout << "GNU implementation: " << total_gnu << "\n";

    auto num_threads = omp_get_max_threads();
    double start_ours = omp_get_wtime();
    for(int i=0; i<ITERATIONS; ++i) {
        std::vector<int> vec(VECTOR_SIZE);
        std::generate(vec.begin(), vec.end(), rand);
        p_partition::quicksort(vec.begin(), vec.end(), comparator, num_threads);

    }
    double end_ours = omp_get_wtime();
    double total_ours = end_ours - start_ours;

    std::cout << "our implementation: " << total_ours << "\n";
    std::cout << "num_threads: " << num_threads << "\n";
}


