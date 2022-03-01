#include <iostream>
#include "library.h"
#include <vector>
#include <parallel/algorithm>
#include <omp.h>


void benchmark(bool csv, bool v0, bool others_only=false) {

    unsigned long VECTOR_SIZE = v0 ? 10000000 : 100000000;
    auto comparator = [](const auto& em1, const auto& em2){ return em1 < em2; };

    for (int config=0; config < (v0 ? 6 : 3); ++config) {
        unsigned long iterations = (unsigned long)pow(10, config);
        unsigned long vector_size = VECTOR_SIZE / (unsigned long)pow(10, config);

        if (!csv || others_only) {
            // gnu sort
            double total_gnu = 0.;
            for (int i = 0; i < iterations; ++i) {
                std::vector<int> vec(vector_size);
                std::generate(vec.begin(), vec.end(), rand);

                double start_gnu = omp_get_wtime();
                __gnu_parallel::sort(vec.begin(), vec.end(), comparator);
                double end_gnu = omp_get_wtime();

                total_gnu += end_gnu - start_gnu;
            }
            total_gnu /= (double) iterations;

            // std::sort
            double total_std = 0.;
            for (int i = 0; i < iterations; ++i) {
                std::vector<int> vec(vector_size);
                std::generate(vec.begin(), vec.end(), rand);

                double start_std = omp_get_wtime();
                std::sort(vec.begin(), vec.end(), comparator);
                double end_std = omp_get_wtime();

                total_std += end_std - start_std;
            }
            total_std /= (double) iterations;

            std::cout << "config " << config << " (iterations: " << iterations << ", size: " << vector_size << "): \t"
                      << "GNU: " << total_gnu << ",\t" << "std: " << total_std << "\n";
        }

        if (!others_only) {
            for (int blockBytes = v0 ? 200 : 20000; blockBytes <= (v0 ? 60000 : 8000000); blockBytes *= 2) {
                p_partition::BLOCK_BYTES = blockBytes;

                unsigned long blockSize = blockBytes / sizeof(int);
                // blockSizes larger than the vector size don't need to be tested, as they behave all the same
                if (blockSize > vector_size) continue;

                for (int bFactor = 1; bFactor <= 64 * 64 * (v0 ? 1 : 8); bFactor *= 2) {

                    unsigned long breakoff = blockSize * bFactor / 64;
                    // breakoffs larger than the vector size don't need to be tested, as they behave all the same
                    if (breakoff > vector_size) continue;

                    p_partition::breakoffFactor = bFactor;
                    // our sort
                    double total_ours = 0.;
                    for (int i = 0; i < iterations; ++i) {
                        std::vector<int> vec(vector_size);
                        std::generate(vec.begin(), vec.end(), rand);

                        double start_ours = omp_get_wtime();
                        p_partition::quicksort(vec.begin(), vec.end(), comparator);
                        double end_ours = omp_get_wtime();

                        total_ours += end_ours - start_ours;
                    }
                    total_ours /= (double) iterations;

                    if (!csv)
                        std::cout << "blockBytes: " << blockBytes << "    \t, blockSize: " << blockSize
                                  << "    \t, breakoff: "
                                  << breakoff << "    \t, bFactor: " << bFactor << "    \t ours: " << total_ours
                                  << std::endl;
                    else // csv output version:
                        std::cout << config << "," << blockBytes << "," << blockSize << "," << breakoff << ","
                                  << bFactor << "," << total_ours << std::endl;
                }
            }
        }
    }
}

int main(int argc, char *argv[])
{
    // BENCHMARK:
    srand (time(nullptr)); // randomly seed rand

    benchmark(true, true, false);
}