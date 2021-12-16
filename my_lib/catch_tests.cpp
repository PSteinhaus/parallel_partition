#include "catch.h"
#include "library.h"
#include <vector>

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

TEST_CASE("test parallel_partition Phase1 ", "[correctness]"){
    std::vector<float> v = {9,8,7,6,5,4,3,2,1,0};
    std::vector<float> expected = {0,1,2,3,4,5,6,7,8,9};
    float pivot = 4.5;

    // test neutralize
    int leftNeutralized, rightNeutralized;
    p_partition::parallel_partition_phase(v.begin(),
                                          v.end(),
                                          [pivot](const auto& em){ return em < pivot; },
                                          &leftNeutralized,
                                          &rightNeutralized,
                                          1,1);
    for (int i=0; i<v.size(); i++){
        REQUIRE(v[i] == expected[i]);
    }

}