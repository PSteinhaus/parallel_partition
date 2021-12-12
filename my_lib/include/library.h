#ifndef PARALLEL_PARTITION_LIBRARY_H
#define PARALLEL_PARTITION_LIBRARY_H

#include <vector>
#include <omp.h>

namespace p_partition  {
    const int BOTH = 0;
    const int LEFT = 1;
    const int RIGHT = 2;

    //returns 0 if Both neutralized, 1 if left ist neutralized, 2 if right is neutralized
    template <typename ForwardIt, typename UnaryPredicate>
    int neutralize(ForwardIt left, ForwardIt right, int blocksize, UnaryPredicate p){
        int i=0 , j=0;
        do {
            for (; i < blocksize; i++){
                if (!p(*left))
                    break;
                left = std::next(left);
            }
            for (; j < blocksize; j++) {
                if (p(*right))
                    break;
                right = std::next(right);
            }
            if ((i == blocksize) || (j == blocksize)){
                break;
            }

            std::swap(*left, *right);
            i++; left = std::next(left);
            j++; right = std::next(right);
        } while ( i < blocksize && j < blocksize );

        std::cout << i << " " << j << std::endl;
        if (i == blocksize) {
            if (j == blocksize)
                return BOTH;
            return LEFT;
        }
        return RIGHT;
    }

    const int BLOCK_SIZE = 2;

    template <typename It1, typename It2>
    bool block_available(It1 target, It2 counterpart) {
        return std::distance(target, counterpart) >= BLOCK_SIZE-1;  // TODO: does this work? The reference says that negative distances are only possible for random access iterators (such as pointers)...
    }

    template <typename It>
    It get_block(It *iter) {
        It block = *iter;
        *iter = std::next(*iter, BLOCK_SIZE);
        return block;
    }

    // returns the remaining blocks
    template <typename ForwardIt, typename BackwardIt, typename UnaryPredicate>
    std::vector<ForwardIt> parallel_partition_phase(ForwardIt left, BackwardIt right, UnaryPredicate p, int *leftNeutralized, int *rightNeutralized) {
        //auto pivot = *std::next(left, std::distance(left,right)/2); // naive pivot
        std::vector<ForwardIt> remainingBlocks;
        int ln, rn;
#pragma omp parallel reduction(+: ln, rn)
        {
#pragma omp single
            {
                remainingBlocks = new std::vector<ForwardIt>(omp_get_num_threads());
            }
            ForwardIt leftBlock, rightBlock;
            bool gotLeftBlock, gotRightBlock;
            int leftCounter = 0, rightCounter = 0;
#pragma omp critical
            {   // get your first left block
                gotLeftBlock = block_available(left, right);
                if (gotLeftBlock) {
                    leftBlock = get_block(&left);
                }
                // get your first right block
                gotRightBlock = block_available(right, left);
                if (gotRightBlock) {
                    rightBlock = get_block(&right);
                }
            }
            while (gotLeftBlock && gotRightBlock) {
                int side = neutralize(leftBlock, rightBlock, p);
#pragma omp critical
                {
                    // try to get new blocks, depending on which side was completed in the neutralize call
                    switch (side) {
                        case BOTH:
                            ++leftCounter;
                            gotLeftBlock = block_available(left, right);
                            if (gotLeftBlock) {
                                leftBlock = get_block(&left);
                            }
                            // don't break, just continue with the RIGHT case to get a right block as well
                        case RIGHT:
                            ++rightCounter;
                            gotRightBlock = block_available(right, left);
                            if (gotRightBlock) {
                                rightBlock = get_block(&right);
                            }
                            break;
                        case LEFT:
                            ++leftCounter;
                            gotLeftBlock = block_available(left, right);
                            if (gotLeftBlock) {
                                leftBlock = get_block(&left);
                            }
                            break;
                    }
                }
            }
#pragma omp critical
            {
                if (gotLeftBlock) remainingBlocks.push_back(leftBlock);
                if (gotRightBlock) remainingBlocks.push_back(rightBlock);
            }
            ln += leftCounter * BLOCK_SIZE;
            rn += rightCounter * BLOCK_SIZE;
        }
        *leftNeutralized = ln;
        *rightNeutralized = ln;
        return remainingBlocks;
    }
}

#endif //PARALLEL_PARTITION_LIBRARY_H
