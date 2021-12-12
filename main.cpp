#include <iostream>
#include "library.h"
#include <vector>



int main(int argc, char *argv[])
{
    std::vector<float> v = {9,8,7,6,5,4,3,2,1,0};
    int blocklength = 5;
    float pivot = 4.5;

    // test neutralize
    int ret = neutralize(v.begin(),v.end()-5, blocklength, [pivot](const auto& em){ return em < pivot; });

    std::cout << "ret = " << ret << std::endl;

    for (int i=0; i<v.size(); i++)
        std::cout << v[i] << " ";

}


