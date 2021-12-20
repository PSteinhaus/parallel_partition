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

    template <typename ForwardIt>
    ForwardIt get_right_block(ForwardIt leftBegin, int *rightTaken, int size, int blockSize) {
        ForwardIt block = std::next(leftBegin, size - blockSize - *rightTaken);
        *rightTaken += blockSize;
        return block;
    }

    bool get_left_block_from_remaining(std::vector<int> *remainingBlocks, int ln, int *left){
        std::vector<int>::iterator remBegin = (*remainingBlocks).begin(), remEnd = (*remainingBlocks).end();
        if(!(*remainingBlocks).empty()) {
            while (true) {
                if (*remBegin < ln) {
                    *left = *remBegin;
                    (*remainingBlocks).erase(remBegin);
                    return true;
                }
                //check if List ended
                if (std::distance(remBegin, remEnd) > 1) {
                    remBegin = std::next(remBegin);
                } else { return false; }
            }
        }
        return false;
    }

    bool get_right_block_from_remaining(std::vector<int> *remainingBlocks, int rn, int *right){
        std::vector<int>::iterator remBegin = (*remainingBlocks).begin(), remEnd = (*remainingBlocks).end();
        if(!(*remainingBlocks).empty()) {
            while (true) {
                if (*remBegin >= rn) {
                    *right = *remBegin;
                    (*remainingBlocks).erase(remBegin);
                    return true;
                } else if (std::distance(remBegin, remEnd) > 1) {
                    remBegin = std::next(remBegin);
                } else { return false; }
            }
        }
        return false;
    }

    int get_blocks_from_remaining(std::vector<int> *remainingBlocks, int ln, int rn, int *left, int *right){
        bool foundBlock = false;

        foundBlock = get_left_block_from_remaining(remainingBlocks, ln, left);
        if(foundBlock == false){return false;}
        foundBlock = get_right_block_from_remaining(remainingBlocks, rn, right);
        if(foundBlock == false){
            (*remainingBlocks).push_back(*left);
            return false;
        }
        return true;
    }


    // returns the remaining blocks
    template <typename ForwardIt, typename UnaryPredicate>
    std::vector<int> parallel_partition_phase_one(ForwardIt left, ForwardIt afterLast, UnaryPredicate p, int numThreads, int size, int blockSize, int* leftNeutralized, int *rightNeutralized) {
        int leftTaken = 0, rightTaken = 0;  // Share two variables to count how many blocks have already been taken in by working threads from each side.
        // This works as long as we take care that only one thread ever uses them at once.
        int ln = 0, rn = 0;

        std::vector<int> remainingBlocks;
        std::mutex taken_mtx;   // This mutex protects the functions where leftTaken and rightTaken are read or modified.


#pragma omp parallel reduction(+: ln, rn) num_threads(numThreads)
        {
            if (omp_get_thread_num() == 1) {
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
                if (gotRightBlock) remainingBlocks.push_back(size - posRightBlock-blockSize);
            }
            ln += leftCounter * blockSize;
            rn += rightCounter * blockSize;
        }
        *leftNeutralized = ln;
        *rightNeutralized = rn;

        return remainingBlocks;
    }

    template <typename ForwardIt, typename UnaryPredicate>
    void parallel_partition_phase_two(ForwardIt left, ForwardIt afterLast, UnaryPredicate p, int size, int blockSize, int ln, int rn, std::vector<int> remainingBlocks){

        bool gotBlocks;
        int leftFromRemaining, rightFromRemaining;
        ForwardIt leftBlock, rightBlock;
        int blockPos;

        //try to get left and right block and neutralize both
        gotBlocks = get_blocks_from_remaining(&remainingBlocks, ln, size-rn, &leftFromRemaining,
                                              &rightFromRemaining);
        while(gotBlocks == true) {
            blockPos = leftFromRemaining;
            std::cout << "leftBlockPos: " << blockPos << std::endl;
            leftBlock = get_left_block(left, &blockPos, blockSize);

            blockPos = rightFromRemaining;
            std::cout << "rightBlockPos: " << blockPos << std::endl;
            rightBlock = get_left_block(left, &blockPos, blockSize);

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
            gotBlocks = get_blocks_from_remaining(&remainingBlocks, ln, size-rn, &leftFromRemaining,
                                                  &rightFromRemaining);
        }

        int usedBlocksLeft = ln, usedBlocksRight = size-rn-blockSize;
        std::cout << "ln: " << ln << " rn: " << rn << std::endl;

        //try to get left block for swap
        gotBlocks = get_left_block_from_remaining(&remainingBlocks, ln, &leftFromRemaining);
        while(gotBlocks == true){


            //get left Block
            blockPos = leftFromRemaining;
            std::cout << "leftBlockPos: " << blockPos << std::endl;
            leftBlock = get_left_block(left, &blockPos, blockSize);

            //get right Block
            int pos = usedBlocksLeft;
            int shift = 0;
            while(std::find(remainingBlocks.begin(), remainingBlocks.end(), pos+shift) != remainingBlocks.end()){
                shift += blockSize;
            }
            pos = pos + shift;
            std::cout << "left swap with pos: " << pos << std::endl;
            rightBlock = get_left_block(left, &pos, blockSize);

            //swap
            for(int i = 0; i< blockSize; i++){
                std::cout << "leftBlock: " << *leftBlock << " rightBlock: " << *rightBlock << std::endl;
                std::swap(*leftBlock, *rightBlock);
                leftBlock = std::next(leftBlock);
                rightBlock = std::next(rightBlock);
            }
            usedBlocksLeft += shift+blockSize;
            gotBlocks = get_left_block_from_remaining(&remainingBlocks, ln, &leftFromRemaining);
        }

        //try to get right block for swap
        gotBlocks = get_right_block_from_remaining(&remainingBlocks, size-rn, &rightFromRemaining);
        while(gotBlocks == true) {
            //get right Block
            blockPos = rightFromRemaining;
            std::cout << "rightBlockPos: " << blockPos << std::endl;
            rightBlock = get_left_block(left, &blockPos, blockSize);

            //get left Block
            int pos = usedBlocksRight;
            int shift = 0;
            while(std::find(remainingBlocks.begin(), remainingBlocks.end(), pos-shift) != remainingBlocks.end()){
                shift = blockSize;
            }
            pos = pos-shift;
            std::cout << "right swap with pos: " << pos << std::endl;
            leftBlock = get_left_block(left, &pos, blockSize);

            //swap
            for(int i = 0; i< blockSize; i++){
                std::cout << "leftBlock: " << *leftBlock << " rightBlock: " << *rightBlock << std::endl;
                std::swap(*leftBlock, *rightBlock);
                leftBlock = std::next(leftBlock);
                rightBlock = std::next(rightBlock);
            }
            usedBlocksRight -= shift+blockSize;
            gotBlocks = get_right_block_from_remaining(&remainingBlocks, size-rn, &rightFromRemaining);
        }


        std::cout << "all remaining Blocks: ";
        for (auto i: remainingBlocks)
            std::cout << i << ' ';
        std::cout<< std::endl;




        //sort remaining values
        ForwardIt leftIt = std::next(left,ln-1);
        ForwardIt rightIt = std::next(left, size-rn);

        std::partition(leftIt, rightIt, p);



    }

    template <typename ForwardIt, typename UnaryPredicate>
    void parallel_partition(ForwardIt left, ForwardIt afterLast, UnaryPredicate p, int numThreads){
        //TODO maybe use long to make user there is no overflow on big arrays
        int size = std::distance(left, afterLast);
        int ln, rn;

        // TODO find a way to calculate a good blookSize
        int blockSize = 3;


        //TODO maybe use list ?
        //TODO maybe only give pointers to functions
        std::vector<int> remainingBlocks;

        //print array
        std::cout << "all elements in Input: ";
        ForwardIt temp = left;
        while(temp != afterLast){
            std::cout << *temp << ' ';
            temp = std::next(temp);
        }
        std::cout << std::endl;


        remainingBlocks = parallel_partition_phase_one(left, afterLast, p, numThreads, size, blockSize, &ln, &rn);

        //print some info
        std::cout << std::endl;
        std::cout << "ln: "<< ln << " rn: " << rn << std::endl;
        std::cout << "all remaining Blocks: ";
        for (auto i: remainingBlocks)
            std::cout << i << ' ';
        std::cout<< std::endl;

        //print array
        std::cout << "all elements in Input: ";
        temp = left;
        while(temp != afterLast){
            std::cout << *temp << ' ';
            temp = std::next(temp);
        }
        std::cout << std::endl;

        parallel_partition_phase_two(left, afterLast, p, size, blockSize, ln, rn, remainingBlocks);

        //print array
        std::cout << std::endl;
        std::cout << "all elements in Input: ";
        temp = left;
        while(temp != afterLast){
            std::cout << *temp << ' ';
            temp = std::next(temp);
        }
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
