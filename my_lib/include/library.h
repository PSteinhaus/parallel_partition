#ifndef PARALLEL_PARTITION_LIBRARY_H
#define PARALLEL_PARTITION_LIBRARY_H

#include <vector>
#include <omp.h>
#include <atomic>
#include <mutex>

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

        //std::cout << i << " " << j << std::endl;
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
        ForwardIt block = std::next(leftBegin, size - BLOCK_SIZE - *rightTaken);
        *rightTaken += BLOCK_SIZE;
        return block;
    }

    // returns the remaining blocks
    template <typename ForwardIt, typename UnaryPredicate>
    std::vector<ForwardIt> parallel_partition_phase(ForwardIt left, ForwardIt afterLast, UnaryPredicate p, int *leftNeutralized, int *rightNeutralized) {
        int leftTaken = 0, rightTaken = 0;  // Share two variables to count how many blocks have already been taken in by working threads from each side.
                                            // This works as long as we take care that only one thread ever uses them at once.
        int size = std::distance(left, afterLast);
        int ln = 0, rn = 0;
        std::vector<ForwardIt> remainingBlocks;
        std::mutex taken_mtx;   // This mutex protects the functions where leftTaken and rightTaken are read or modified.

#pragma omp parallel reduction(+: ln, rn)
        {
            ForwardIt leftBlock, rightBlock;
            bool gotLeftBlock, gotRightBlock;
            int leftCounter = 0, rightCounter = 0;
            {
                taken_mtx.lock();
                // get your first left block
                gotLeftBlock = block_available(leftTaken, rightTaken, size);
                if (gotLeftBlock) {
                    leftBlock = get_left_block(left, &leftTaken);
                }
                // get your first right block
                gotRightBlock = block_available(leftTaken, rightTaken, size);
                if (gotRightBlock) {
                    rightBlock = get_right_block(left, &rightTaken, size);
                }
                taken_mtx.unlock();
            }
            while (gotLeftBlock && gotRightBlock) {
                int side = neutralize(leftBlock, rightBlock, BLOCK_SIZE, p);
                // try to get new blocks, depending on which side was completed in the neutralize call
                switch (side) {
                    case BOTH:
                        ++leftCounter;
                        taken_mtx.lock();
                        gotLeftBlock = block_available(leftTaken, rightTaken, size);
                        if (gotLeftBlock) {
                            leftBlock = get_left_block(left, &leftTaken);
                        }
                        taken_mtx.unlock();
                        // don't break, just continue with the RIGHT case to get a right block as well
                    case RIGHT:
                        ++rightCounter;
                        taken_mtx.lock();
                        gotRightBlock = block_available(leftTaken, rightTaken, size);
                        if (gotRightBlock) {
                            rightBlock = get_right_block(left, &rightTaken, size);
                        }
                        taken_mtx.unlock();
                        break;
                    case LEFT:
                        ++leftCounter;
                        taken_mtx.lock();
                        gotLeftBlock = block_available(leftTaken, rightTaken, size);
                        if (gotLeftBlock) {
                            leftBlock = get_left_block(left, &leftTaken);
                        }
                        taken_mtx.unlock();
                        break;
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
