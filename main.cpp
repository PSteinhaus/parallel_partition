#include <iostream>
#include "library.h"
#include <vector>



int main(int argc, char *argv[])
{
    std::vector<float> v = {9,8,7,6,5,4,3,2,1,0};
    int blocklength = 5;
    float pivot = 4.5;

    //GetMax <int> (x,y);
    std::cout << "test" << std::endl;
    hello();
    int ret = neutralize(v.begin(),v.end()-1,blocklength, pivot);

    std::cout << ret << std::endl;

    for (int i=0; i<v.size(); i++)
        std::cout << v[i] << " ";

}


