#ifndef PARALLEL_PARTITION_LIBRARY_H
#define PARALLEL_PARTITION_LIBRARY_H

#include <vector>
#include <omp.h>
#include <atomic>
#include <mutex>
#include <cstdlib>
#include <iostream>


namespace p_partition  {
    const int BOTH = 0;
    const int LEFT = 1;
    const int RIGHT = 2;
    const int NONE = 3;

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

    template <typename ForwardIt>
    void swapBlocks(ForwardIt left, ForwardIt right, int blockSize){
        for(int i = 0; i< blockSize; i++){
            std::swap(*left, *right);
        }
    }



    inline
    bool block_available(int leftTaken, int rightTaken, int size) {
        return leftTaken + rightTaken < size;
    }

    template <typename ForwardIt>
    ForwardIt get_left_block(ForwardIt leftBegin, int *leftTaken, int blockSize) {
        ForwardIt block = std::next(leftBegin, *leftTaken);
        *leftTaken += blockSize;
        return block;
    }

    //TODO maybe use prev instead of next for faster runtime
    template <typename ForwardIt>
    ForwardIt get_right_block(ForwardIt leftBegin, int *rightTaken, int size, int blockSize) {
        ForwardIt block = std::next(leftBegin, size - blockSize - *rightTaken);
        *rightTaken += blockSize;
        return block;
    }


    int get_blocks_from_remaining(std::vector<int> *remainingBlocks, int ln, int rn, int *left, int *right){
        int blockPos;
        int foundBlocks = NONE;
        bool foundLeft = false, foundRight = false;
        std::vector<int>::iterator remBegin = (*remainingBlocks).begin(), remEnd = (*remainingBlocks).end();

        //go through remainingBlocks and try to find a block for left and right
        while (foundBlocks != LEFT) {
            if (*remBegin < ln && foundBlocks != LEFT) {
                *left = *remBegin;
                //std::cout << "found Left: " << *left << std::endl;
                (*remainingBlocks).erase(remBegin);
                foundBlocks = LEFT;
            }

            //check if List ended
            else if (std::distance(remBegin, remEnd) > 1) {
                remBegin = std::next(remBegin);
            }
            else{break;}
        }
        remBegin = (*remainingBlocks).begin();
        remEnd = (*remainingBlocks).end();
        while((foundBlocks != RIGHT) && (foundBlocks != BOTH)){
            if (*remBegin > rn && foundBlocks != RIGHT) {
                *right = *remBegin;
                //std::cout << "found Right: "<< *right << std::endl;
                (*remainingBlocks).erase(remBegin);
                if(foundBlocks == LEFT){foundBlocks = BOTH;}
                else{foundBlocks = RIGHT;}
            }
            else if(std::distance(remBegin, remEnd) > 1) {
                remBegin = std::next(remBegin);
            }
            else{break;}
        }
        return foundBlocks;
    }


    // returns the remaining blocks
    template <typename ForwardIt, typename UnaryPredicate>
    void parallel_partition(ForwardIt left, ForwardIt afterLast, UnaryPredicate p, int numThreads) {
        int leftTaken = 0, rightTaken = 0;  // Share two variables to count how many blocks have already been taken in by working threads from each side.
                                            // This works as long as we take care that only one thread ever uses them at once.
        int size = std::distance(left, afterLast);
        int ln = 0, rn = 0;
        int blockSize = 3; // TODO find a way to calculate a good blookSize

        //TODO maybe use long to make user there is no overflow on big arrays
        std::vector<int> remainingBlocks; //TODO maybe use std::list ?
        std::mutex taken_mtx;   // This mutex protects the functions where leftTaken and rightTaken are read or modified.

        //Phase 1

#pragma omp parallel reduction(+: ln, rn) num_threads(numThreads)
        {
            if(omp_get_thread_num() == 1){
                std::cout << "number of threads: " << omp_get_num_threads() << std::endl;
            }

            ForwardIt leftBlock, rightBlock;
            bool gotLeftBlock, gotRightBlock;
            int leftCounter = 0, rightCounter = 0;
            int posLeftBlock, posRightBlock;
            {
                taken_mtx.lock();
                // get your first left block
                gotLeftBlock = block_available(leftTaken, rightTaken, size);
                if (gotLeftBlock) {
                    posLeftBlock = leftTaken;
                    leftBlock = get_left_block(left, &leftTaken, blockSize);
                }
                // get your first right block
                gotRightBlock = block_available(leftTaken, rightTaken, size);
                if (gotRightBlock) {
                    posRightBlock = rightTaken;
                    rightBlock = get_right_block(left, &rightTaken, size, blockSize);
                }
                taken_mtx.unlock();
            }
            while (gotLeftBlock && gotRightBlock) {
                int side = neutralize(leftBlock, rightBlock, blockSize, p);
                // try to get new blocks, depending on which side was completed in the neutralize call
                switch (side) {
                    case BOTH:
                        ++leftCounter;
                        taken_mtx.lock();
                        gotLeftBlock = block_available(leftTaken, rightTaken, size);
                        if (gotLeftBlock) {
                            posLeftBlock = leftTaken;
                            leftBlock = get_left_block(left, &leftTaken, blockSize);
                        }
                        taken_mtx.unlock();
                        // don't break, just continue with the RIGHT case to get a right block as well
                    case RIGHT:
                        ++rightCounter;
                        taken_mtx.lock();
                        gotRightBlock = block_available(leftTaken, rightTaken, size);
                        if (gotRightBlock) {
                            posRightBlock = rightTaken;
                            rightBlock = get_right_block(left, &rightTaken, size, blockSize);
                        }
                        taken_mtx.unlock();
                        break;
                    case LEFT:
                        ++leftCounter;
                        taken_mtx.lock();
                        gotLeftBlock = block_available(leftTaken, rightTaken, size);
                        if (gotLeftBlock) {
                            posLeftBlock = leftTaken;
                            leftBlock = get_left_block(left, &leftTaken, blockSize);
                        }
                        taken_mtx.unlock();
                        break;
                }
            }
#pragma omp critical
            {
                if (gotLeftBlock) remainingBlocks.push_back(posLeftBlock);
                if (gotRightBlock) remainingBlocks.push_back(size-posRightBlock);
            }
            ln += leftCounter * blockSize;
            rn += rightCounter * blockSize;
        }

        //print some info
        std::cout << "size remBlock: " << remainingBlocks.size() << std::endl;
        std::cout << "ln: "<< ln << " rn: " << size- rn<< std::endl;
        std::cout << "all remaining Blocks: ";
        for (auto i: remainingBlocks)
            std::cout << i << ' ';
        std::cout<< std::endl;

        //start Phase 2
        int gotBlocks = BOTH;
        int leftFromRemaining, rightFromRemaining;
        int usedBlocksLeft = ln, usedBlocksRight = size-rn-blockSize;
        ForwardIt leftBlock, rightBlock;


        while(gotBlocks != NONE) {
            gotBlocks = NONE;
            if (!remainingBlocks.empty()) {
                gotBlocks = get_blocks_from_remaining(&remainingBlocks, ln, size - rn, &leftFromRemaining,
                                                &rightFromRemaining);
            }
            //get the actual blocks
            if (gotBlocks == LEFT || gotBlocks == BOTH) {
                int blockPos = leftFromRemaining;
                std::cout << "leftBlockPos: " <<blockPos << std::endl;
                leftBlock = get_left_block(left, &blockPos, blockSize);
            }
            if (gotBlocks == RIGHT || gotBlocks == BOTH) {
                int blockPos = rightFromRemaining - blockSize;
                std::cout << "rightBlockPos: " <<blockPos << std::endl;
                rightBlock = get_left_block(left, &blockPos, blockSize);
            }

            //TODO possible optimise to only get the missing Block after neutralize;
            if (gotBlocks == BOTH) {
                int side = neutralize(leftBlock, rightBlock, blockSize, p);
                switch (side) {
                    case BOTH:
                        ln += blockSize;
                        rn += blockSize;
                        break;
                    case RIGHT:
                        rn += blockSize;
                        remainingBlocks.push_back(leftFromRemaining);
                        break;
                    case LEFT:
                        ln += blockSize;
                        remainingBlocks.push_back(rightFromRemaining);
                        break;
                }
            }
            if(gotBlocks == LEFT){
                int pos = usedBlocksLeft;
                while(std::find(remainingBlocks.begin(), remainingBlocks.end(), pos) != remainingBlocks.end()){
                    pos += blockSize;
                }
                std::cout << "left swap with pos: " << pos << std::endl;
                rightBlock = get_left_block(left, &pos, blockSize);
                std::cout << "leftBlock: " << *leftBlock << " rightBlock: " << *rightBlock << std::endl;
                //swapBlocks(leftBlock, rightBlock,blockSize);
                usedBlocksLeft += blockSize;
            }
            if(gotBlocks == RIGHT){
                int pos = usedBlocksRight;
                std::cout << "start search: "<< pos << std::endl;
                while(std::find(remainingBlocks.begin(), remainingBlocks.end(), pos) != remainingBlocks.end()){
                    pos -= blockSize;
                }
                std::cout << "right swap with pos: " << pos << std::endl;
                rightBlock = get_left_block(left, &pos, blockSize);
                std::cout << "leftBlock: " << *leftBlock << " rightBlock: " << *rightBlock << std::endl;
                //swapBlocks(leftBlock, rightBlock,blockSize);
                usedBlocksRight -= blockSize;
            }

        }

        std::cout << "all remaining Blocks: ";
        for (auto i: remainingBlocks)
            std::cout << i << ' ';
        std::cout<< std::endl;

        std::cout << "ln: " << ln << " rn: " << size - rn << std::endl;

        //sort remaining values

        ForwardIt leftIt = std::next(left,ln-1);
        ForwardIt rightIt = std::next(left, size-rn-1);

        std::partition(leftIt, rightIt, p);



    }




    // test func for nested parallelism
    void testFunc(){
        omp_set_nested(2);

#pragma omp parallel for num_threads(5)
        for(int i = 0; i< 5 ; i++){
#pragma omp parallel for num_threads(2)
            for(int j = 0; j<2; j++) {
                int tid = omp_get_thread_num();
                int num = omp_get_num_threads();
                printf("Hello world from omp thread %d %d\n", tid, num);
                _sleep(3000);
            }
        }
    }


}

#endif //PARALLEL_PARTITION_LIBRARY_H
