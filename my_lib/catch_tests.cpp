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
    float pivot = 50;

    std::vector<float> v(numberOfValues);
    std::generate(v.begin(), v.end(), []() {
        return rand() % 100;
    });
    std::vector<float> expected = v;
    std::partition(expected.begin(), expected.end(),[pivot](const auto& em){ return em < pivot;});


    // test neutralize
    int leftNeutralized, rightNeutralized;
    p_partition::parallel_partition(v.begin(),
                                    v.end(),
                                    [pivot](const auto& em){ return em < pivot; },
                                    12);


    for (int i=0; i<v.size(); i++){
        REQUIRE((v[i] < pivot) == (expected[i] < pivot));
    }
}

TEST_CASE("test parallel_partition_phase1", "[correctness]"){
    int numberOfValues = 50;
    int numThreads = 12;
    int blockSize = 3;
    float pivot = 50;
    int ln, rn;


    std::vector<int> remainingBlocks;

    std::vector<float> v(numberOfValues);
    std::generate(v.begin(), v.end(), []() {
        return rand() % 100;
    });

    remainingBlocks = p_partition::parallel_partition_phase_one(v.begin(),v.end(), [pivot](const auto& em){ return em < pivot; },
                                              numThreads, numberOfValues, blockSize, &ln,&rn);

    for (int i=0; i<ln; i+= blockSize){
        bool inRemaining =  std::find(remainingBlocks.begin(), remainingBlocks.end(), i) != remainingBlocks.end();
        //not in remaining => values < pivot
        if(!inRemaining){
            for (int j=0; j<blockSize; j++) {
                REQUIRE(v[i + j] < pivot);
            }
        }
        //in remainig => at least one value > pivot
        else{
            bool temp = false;
            for (int j=0; j<blockSize; j++) {
                temp = temp || (v[i+j] > pivot);
            }
            REQUIRE(temp);
        }
    }

    for (int i=0; i>rn; i+= blockSize){
        int pos = numberOfValues-i-blockSize;
        bool inRemaining =  std::find(remainingBlocks.begin(), remainingBlocks.end(), pos) != remainingBlocks.end();
        //not in remaining => values > pivot
        if(!inRemaining){
            for (int j=0; j<blockSize; j++) {
                REQUIRE(v[pos + j] > pivot);
            }
        }
            //in remainig => at least one value < pivot
        else{
            bool temp = false;
            for (int j=0; j<blockSize; j++) {
                temp = temp || (v[pos+j] < pivot);
            }
            REQUIRE(temp);
        }
    }
}