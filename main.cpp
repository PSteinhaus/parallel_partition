#include <iostream>
#include "library.h"
#include <vector>



int main(int argc, char *argv[])
{
    std::vector<float> v = {9,8,7,6,5,4,3,2,1,0};
    int blocklength = 5;
    float pivot = 4.5;


    // test parallel partition (phase one)  (block size is currently 2)
    int leftNeutralized, rightNeutralized;
    auto remainingBlocks = p_partition::parallel_partition_phase(v.begin(),
                                          v.end(),
                                          [pivot](const auto& em){ return em < pivot; },
                                          &leftNeutralized,
                                          &rightNeutralized,
                                          1,1);
    for (int i=0; i<v.size(); i++)
        std::cout << v[i] << " ";
    std::cout << std::endl;

    for (int i=0; i<remainingBlocks.size(); i++) {
        auto iter = remainingBlocks[i];
        std::cout << "remaining: ";
        for (int j = 0; j < p_partition::BLOCK_SIZE; j++) {
            std::cout << iter[j] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << leftNeutralized << " " << rightNeutralized << std::endl;

    p_partition::testFunc();
}


