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
    p_partition::partition(v.begin(),
                           v.end(),
                           [pivot](const auto &em) { return em < pivot; },
                           1);
    for (int i=0; i<v.size(); i++)
        std::cout << v[i] << " ";
    std::cout << std::endl;


    std::cout << leftNeutralized << " " << rightNeutralized << std::endl;

    p_partition::testFunc();
}


