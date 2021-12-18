#include "catch.h"
#include "library.h"
#include <vector>
#include <random>

TEST_CASE("test neutralize", "[correctness]"){
    std::vector<float> v = {4,8,7,6,5,4,3,2,1,0};
    std::vector<float> expected = {4,4,3,2,1,8,7,6,5,0};
    int blocklength = 5;
    float pivot = 4.5;

    // test neutralize
    p_partition::neutralize(v.begin(),v.begin()+5, blocklength, [pivot](const auto& em){ return em < pivot; });
    for (int i=0; i<v.size(); i++){
        REQUIRE(v[i] == expected[i]);
    }

}

TEST_CASE("test parallel_partition", "[correctness]"){
    int numberOfValues = 50;

    std::vector<float> v(numberOfValues);
    std::generate(v.begin(), v.end(), []() {
        return rand() % 100;
    });

    float pivot = 50;

    // test neutralize
    int leftNeutralized, rightNeutralized;
    p_partition::parallel_partition(v.begin(),
                                    v.end(),
                                    [pivot](const auto& em){ return em < pivot; },
                                    12);

    for (auto i: v)
        std::cout << i << ' ';
    std::cout<< std::endl;

    for (int i=0; i<v.size(); i++){
        REQUIRE(v[i] == 1);
    }

}