#include <iostream>
#include "library.h"
#include <vector>
#include <parallel/algorithm>
#include <omp.h>

int main(int argc, char *argv[])
{
    std::vector<float> v = {9,8,7,6,5,4,3,2,1,0};
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


    // BENCHMARK:
    srand (time(nullptr)); // randomly seed rand

    const unsigned long ITERATIONS = 1;
    const unsigned long VECTOR_SIZE = 10000000;
    auto comparator = [](const auto& em1, const auto& em2){ return em1 < em2; };

    for (int config=0; config<6; ++config) {
        unsigned long iterations =  ITERATIONS * (unsigned long)pow(10, config);
        unsigned long vector_size = VECTOR_SIZE / (unsigned long)pow(10, config);

        // gnu sort
        double start_gnu = omp_get_wtime();
        for(int i=0; i<iterations; ++i) {
            std::vector<int> vec(vector_size);
            std::generate(vec.begin(), vec.end(), rand);
            __gnu_parallel::sort(vec.begin(), vec.end(), comparator);
        }
        double end_gnu = omp_get_wtime();
        double total_gnu = end_gnu - start_gnu;

        // std::sort
        double start_std = omp_get_wtime();
        for(int i=0; i<iterations; ++i) {
            std::vector<int> vec(vector_size);
            std::generate(vec.begin(), vec.end(), rand);
            std::sort(vec.begin(), vec.end(), comparator);
        }
        double end_std = omp_get_wtime();
        double total_std = end_std - start_std;

        //std::cout << "config " << config << " (iterations: " << iterations << ", size: " << vector_size << "): \t" << "GNU: " << total_gnu << ",\t" << "std: " << total_std << "\n";

        for (int blockBytes=200; blockBytes<=60000; blockBytes *= 2) {
            p_partition::BLOCK_BYTES = blockBytes;

            unsigned long blockSize = blockBytes / sizeof(int);
            // blockSizes larger than the vector size don't need to be tested, as they behave all the same
            if (blockSize > vector_size) continue;

            for (int bFactor = 1; bFactor <= 64 * 64; bFactor *= 2) {

                unsigned long breakoff = blockSize * bFactor / 64;
                // breakoffs larger than the vector size don't need to be tested, as they behave all the same
                if (breakoff > vector_size) continue;

                p_partition::breakoffFactor = bFactor;
                // our sort
                double start_ours = omp_get_wtime();
                for (int i = 0; i < iterations; ++i) {
                    std::vector<int> vec(vector_size);
                    std::generate(vec.begin(), vec.end(), rand);
                    p_partition::quicksort(vec.begin(), vec.end(), comparator);
                }
                double end_ours = omp_get_wtime();
                double total_ours = end_ours - start_ours;

                //std::cout << "blockBytes: " << blockBytes << "    \t, blockSize: " << blockSize << "    \t, breakoff: " << breakoff << "    \t, bFactor: " << bFactor << "    \t ours: " << total_ours << std::endl;
                // csv output version:
                std::cout << config << "," << blockBytes << "," << blockSize << "," << breakoff << "," << bFactor << "," << total_ours << std::endl;
            }
        }
    }


}


