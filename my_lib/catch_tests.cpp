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

TEST_CASE("test partition", "[correctness]"){
    int numberOfValues = 50;
    float pivot = 50;

    std::vector<float> v(numberOfValues);
    std::generate(v.begin(), v.end(), []() {
        return rand() % 100;
    });
    std::vector<float> expected = v;
    std::partition(expected.begin(), expected.end(),[pivot](const auto& em){ return em < pivot;});

    p_partition::partition(v.begin(),
                           v.end(),
                           [pivot](const auto &em) { return em < pivot; },
                           12);

    for (int i=0; i<v.size(); i++){
        REQUIRE((v[i] < pivot) == (expected[i] < pivot));
    }
}

TEST_CASE("test parallel_partition_phase_1", "[correctness]"){
    int numberOfValues = 50;
    int numThreads = 12;
    int blockSize = 6;
    float pivot = 50;
    int ln, rn;


    std::vector<int> remainingBlocks;

    std::vector<float> v(numberOfValues);
    std::generate(v.begin(), v.end(), []() {
        return rand() % 100;
    });

    //print array
    std::cout << "all elements in Input: ";
    auto temp = v.begin();
    while(temp != v.end()){
        std::cout << *temp << ' ';
        temp = std::next(temp);
    }
    std::cout << std::endl;

    remainingBlocks = p_partition::parallel_partition_phase_one(v.begin(),v.end(), [pivot](const auto& em){ return em < pivot; },
                                              numThreads, numberOfValues, blockSize, &ln,&rn);

    //print some info
    std::cout << std::endl;
    std::cout << "ln: "<< ln << " rn: " << numberOfValues- rn<< std::endl;
    std::cout << "all remaining Blocks: ";
    for (auto i: remainingBlocks)
        std::cout << i << ' ';
    std::cout<< std::endl;

    //print array
    std::cout << "all elements in Output: ";
    temp = v.begin();
    while(temp != v.end()){
        std::cout << *temp << ' ';
        temp = std::next(temp);
    }
    std::cout << std::endl;

    for (int i=0; i<ln; i+= blockSize){
        bool inRemaining =  std::find(remainingBlocks.begin(), remainingBlocks.end(), i) != remainingBlocks.end();
        //not in remaining => values < pivot
        if(!inRemaining){
            for (int j=0; j<blockSize; j++) {

                REQUIRE(v[i + j] < pivot);
            }
        }
        else{
            ln += blockSize;
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
        else{
            rn += blockSize;
        }
    }
}

TEST_CASE("test parallel_partition_phase_2", "[correctness]"){
    int numberOfValues = 50;
    int blockSize = 3;
    float pivot = 50;
    int ln = 24, rn = 18;

    std::vector<float> v = {41, 11, 34, 0, 41, 24, 21, 16, 18, 11, 5, 45, 35, 27, 3, 12, 38, 42, 27, 36,
                            47, 4, 2, 95, 92, 82, 78, 58, 62, 95, 91, 26, 71, 53, 69, 91, 67, 99, 81, 94,
                            61, 64, 22, 33, 73, 64, 69, 67, 53, 68};

    std::vector<int> remainingBlocks = {41, 29, 32};
    p_partition::parallel_partition_phase_two(v.begin(), v.end(), [pivot](const auto& em){ return em < pivot;},
                                 numberOfValues, blockSize, ln, rn, remainingBlocks);

    std::vector<float> expected = v;
    std::partition(expected.begin(), expected.end(),[pivot](const auto& em){ return em < pivot;});

    for (auto i: v)
        std::cout << i << ' ';
    std::cout<< std::endl;

    for (int i=0; i<v.size(); i++){
        if(expected[i] < pivot){
            REQUIRE(v[i] < pivot);
        }
        else{
            REQUIRE(v[i] >= pivot);
        }
    }
}

TEST_CASE("test quicksort", "[correctness]"){
    int numberOfValues = 500;
    int numThreads = 4;

    std::vector<float> v(numberOfValues);
    std::generate(v.begin(), v.end(), []() {
        return rand() % 100;
    });

    std::vector<float> vCopy(v);

    //print array
    std::cout << "all elements in Input: ";
    auto temp = v.begin();
    while(temp != v.end()){
        std::cout << *temp << ' ';
        temp = std::next(temp);
    }
    std::cout << std::endl;

    auto comparator = [](const auto& em1, const auto& em2){ return em1 < em2; };

    p_partition::quicksort(v.begin(),v.end(), comparator, numThreads);
    // for comparison sort the copy with std::sort
    std::sort(vCopy.begin(), vCopy.end(), comparator);

    //print output
    std::cout << "all elements in Output: ";
    temp = v.begin();
    while(temp != v.end()){
        std::cout << *temp << ' ';
        temp = std::next(temp);
    }
    std::cout << std::endl;

    for (int i=0; i<numberOfValues; i++){
        REQUIRE(v[i] == vCopy[i]);
    }
}

TEST_CASE("test nth_element", "[correctness]"){
    int numberOfValues = 500;
    int numThreads = 4;
    int n = rand() % numberOfValues;

    std::vector<float> v(numberOfValues);
    std::generate(v.begin(), v.end(), []() {
        return rand() % 100;
    });

    std::vector<float> vCopy(v);

    //print array
    std::cout << "all elements in Input: ";
    auto temp = v.begin();
    while(temp != v.end()){
        std::cout << *temp << ' ';
        temp = std::next(temp);
    }
    std::cout << std::endl;

    auto comparator = [](const auto& em1, const auto& em2){ return em1 < em2; };

    p_partition::nth_element(v.begin(), std::next(v.begin(), n), v.end(), comparator, numThreads);
    // for comparison use std::nth_element on the copy
    std::nth_element(vCopy.begin(), std::next(vCopy.begin(), n), vCopy.end(), comparator);

    //print output
    std::cout << "all elements in Output: ";
    temp = v.begin();
    while(temp != v.end()){
        std::cout << *temp << ' ';
        temp = std::next(temp);
    }
    std::cout << std::endl;

    // check that all smaller indexes hold smaller values
    for (int i=0; i<n; i++){
        REQUIRE(v[i] <= v[n]);
    }

    // check that n_th is actually the n_th element of the sorted list
    REQUIRE(v[n] == vCopy[n]);
}