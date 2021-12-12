#ifndef PARALLEL_PARTITION_LIBRARY_H
#define PARALLEL_PARTITION_LIBRARY_H

#include <vector>
#include <omp.h>
#include <atomic>

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

    inline
    bool block_available(int leftTaken, int rightTaken, int size) {
        return leftTaken + rightTaken < size;
    }

    template <typename ForwardIt>
    ForwardIt get_left_block(ForwardIt leftBegin, int *leftTaken) {
        ForwardIt block = std::next(leftBegin, *leftTaken);
        *leftTaken += BLOCK_SIZE;
        return block;
    }

    template <typename ForwardIt>
    ForwardIt get_right_block(ForwardIt leftBegin, int *rightTaken, int size) {
        ForwardIt block = std::next(leftBegin, size - 1 - *rightTaken);
        *rightTaken += BLOCK_SIZE;
        return block;
    }

    // returns the remaining blocks
    template <typename ForwardIt, typename UnaryPredicate>
    std::vector<ForwardIt> parallel_partition_phase(ForwardIt left, ForwardIt afterLast, UnaryPredicate p, int *leftNeutralized, int *rightNeutralized) {
        int leftTaken = 0, rightTaken = 0;
        int size = std::distance(left, afterLast);
        int ln = 0, rn = 0;
        std::vector<ForwardIt> remainingBlocks;

#pragma omp parallel reduction(+: ln, rn)
        {
            ForwardIt leftBlock, rightBlock;
            bool gotLeftBlock, gotRightBlock;
            int leftCounter = 0, rightCounter = 0;
#pragma omp critical
            {   // get your first left block
                gotLeftBlock = block_available(leftTaken, rightTaken, size);
                if (gotLeftBlock) {
                    leftBlock = get_left_block(left, &leftTaken);
                }
                // get your first right block
                gotRightBlock = block_available(leftTaken, rightTaken, size);
                if (gotRightBlock) {
                    rightBlock = get_right_block(left, &rightTaken, size);
                }
            }
            while (gotLeftBlock && gotRightBlock) {
                int side = neutralize(leftBlock, rightBlock, BLOCK_SIZE, p);
#pragma omp critical
                {
                    // try to get new blocks, depending on which side was completed in the neutralize call
                    switch (side) {
                        case BOTH:
                            ++leftCounter;
                            gotLeftBlock = block_available(leftTaken, rightTaken, size);
                            if (gotLeftBlock) {
                                leftBlock = get_left_block(left, &leftTaken);
                            }
                            // don't break, just continue with the RIGHT case to get a right block as well
                        case RIGHT:
                            ++rightCounter;
                            gotRightBlock = block_available(leftTaken, rightTaken, size);
                            if (gotRightBlock) {
                                rightBlock = get_right_block(left, &rightTaken, size);
                            }
                            break;
                        case LEFT:
                            ++leftCounter;
                            gotLeftBlock = block_available(leftTaken, rightTaken, size);
                            if (gotLeftBlock) {
                                leftBlock = get_left_block(left, &leftTaken);
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
        *rightNeutralized = rn;
        return remainingBlocks;
    }
}

#endif //PARALLEL_PARTITION_LIBRARY_H
