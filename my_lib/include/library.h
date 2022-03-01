#ifndef PARALLEL_PARTITION_LIBRARY_H
#define PARALLEL_PARTITION_LIBRARY_H

#include <vector>
#include <omp.h>
#include <atomic>
#include <mutex>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <stack>

namespace p_partition  {
    const int BOTH = 0;
    const int LEFT = 1;
    const int RIGHT = 2;
    static size_t BLOCK_BYTES = 1280000;
    static size_t breakoffFactor = 1;

    //returns 0 if Both neutralized, 1 if left ist neutralized, 2 if right is neutralized
    template <typename ForwardIt, typename UnaryPredicate>
    int neutralize(ForwardIt left, ForwardIt right, size_t blocksize, UnaryPredicate p){
        size_t i=0 , j=0;
        do {
            //std::cout << "1" << std::endl;
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

    inline
    bool block_available(size_t leftTaken, size_t rightTaken, size_t size, size_t blockSize) {
        return size - (leftTaken + rightTaken) >= blockSize;
    }

    template <typename ForwardIt>
    ForwardIt get_left_block(ForwardIt leftBegin, size_t *leftTaken, size_t blockSize) {
        ForwardIt block = std::next(leftBegin, *leftTaken);
        *leftTaken += blockSize;
        return block;
    }

    template <typename ForwardIt>
    ForwardIt get_right_block(ForwardIt leftBegin, size_t *rightTaken, size_t size, size_t blockSize) {
        ForwardIt block = std::next(leftBegin, size - blockSize - *rightTaken);
        *rightTaken += blockSize;
        return block;
    }

    bool get_left_block_from_remaining(std::vector<size_t> *remainingBlocks, size_t ln, size_t *left){
        std::vector<size_t>::iterator remBegin = (*remainingBlocks).begin(), remEnd = (*remainingBlocks).end();
        while (remBegin != remEnd) {
            //std::cout << "2" << std::endl;
            if (*remBegin < ln) {
                *left = *remBegin;
                (*remainingBlocks).erase(remBegin);
                return true;
            }
            remBegin = std::next(remBegin);
        }
        return false;
    }

    bool get_right_block_from_remaining(std::vector<size_t> *remainingBlocks, size_t rn, size_t *right){
        std::vector<size_t>::iterator remBegin = (*remainingBlocks).begin(), remEnd = (*remainingBlocks).end();
        while (remBegin != remEnd) {
            //std::cout << "3" << std::endl;
            if (*remBegin >= rn) {
                *right = *remBegin;
                (*remainingBlocks).erase(remBegin);
                return true;
            }
            remBegin = std::next(remBegin);
        }
        return false;
    }

    int get_blocks_from_remaining(std::vector<size_t> *remainingBlocks, size_t ln, size_t rn, size_t *left, size_t *right){
        bool foundBlock = get_left_block_from_remaining(remainingBlocks, ln, left);
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
    std::vector<size_t> parallel_partition_phase_one(ForwardIt left, UnaryPredicate p, int numThreads, size_t size, size_t blockSize, size_t* leftNeutralized, size_t *rightNeutralized) {
        size_t leftTaken = 0, rightTaken = 0;  // Share two variables to count how many blocks have already been taken in by working threads from each side.
        // This works as long as we take care that only one thread ever uses them at once.
        size_t ln = 0, rn = 0;

        std::vector<size_t> remainingBlocks;
        std::mutex taken_mtx;   // This mutex protects the functions where leftTaken and rightTaken are read or modified.

#pragma omp parallel reduction(+: ln, rn) num_threads(numThreads)
        {
            /*
            if (omp_get_thread_num() == 1) {
                std::cout << "number of threads: " << omp_get_num_threads() << std::endl;
            }
            */
            ForwardIt leftBlock, rightBlock;
            bool gotLeftBlock, gotRightBlock;
            size_t leftCounter = 0, rightCounter = 0;
            size_t posLeftBlock, posRightBlock;
            // DEBUG
            //int t_num = omp_get_thread_num();
            {
                taken_mtx.lock();
                // get your first left block
                gotLeftBlock = block_available(leftTaken, rightTaken, size, blockSize);
                if (gotLeftBlock) {
                    posLeftBlock = leftTaken;
                    leftBlock = get_left_block(left, &leftTaken, blockSize);
                    // DEBUG
                    //std::cout << "thread "<< t_num << " took left block: " << posLeftBlock << std::endl;
                }
                // get your first right block
                gotRightBlock = block_available(leftTaken, rightTaken, size, blockSize);
                if (gotRightBlock) {
                    posRightBlock = rightTaken;
                    rightBlock = get_right_block(left, &rightTaken, size, blockSize);
                    // DEBUG
                    //std::cout << "thread "<< t_num << " took right block: " <<  size - posRightBlock - blockSize << std::endl;
                }
                taken_mtx.unlock();
            }
            while (gotLeftBlock && gotRightBlock) {
                //std::cout << "4" << std::endl;
                int side = neutralize(leftBlock, rightBlock, blockSize, p);
                // try to get new blocks, depending on which side was completed in the neutralize call
                switch (side) {
                    case BOTH:
                        ++leftCounter;
                        ++rightCounter;
                        taken_mtx.lock();
                        // DEBUG
                        //std::cout << "thread "<< t_num << " finished left block: " << posLeftBlock << std::endl;
                        gotLeftBlock = block_available(leftTaken, rightTaken, size, blockSize);
                        if (gotLeftBlock) {
                            posLeftBlock = leftTaken;
                            leftBlock = get_left_block(left, &leftTaken, blockSize);
                            // DEBUG
                            //std::cout << "thread "<< t_num << " took left block: " << posLeftBlock << std::endl;
                        }
                        // don't break, just continue with the RIGHT case to get a right block as well
                        // jump to where the lock is already taken to avoid having to un- and relock it
                        goto right_locked;
                    case RIGHT:
                        ++rightCounter;
                        taken_mtx.lock();
                        right_locked:
                        // DEBUG
                        //std::cout << "thread "<< t_num << " finished right block: " << size - posRightBlock - blockSize << std::endl;
                        gotRightBlock = block_available(leftTaken, rightTaken, size, blockSize);
                        if (gotRightBlock) {
                            posRightBlock = rightTaken;
                            rightBlock = get_right_block(left, &rightTaken, size, blockSize);
                            // DEBUG
                            //std::cout << "thread "<< t_num << " took right block: " <<  size - posRightBlock - blockSize << std::endl;
                        }
                        taken_mtx.unlock();
                        break;
                    case LEFT:
                        ++leftCounter;
                        taken_mtx.lock();
                        // DEBUG
                        //std::cout << "thread "<< t_num << " finished left block: " << posLeftBlock << std::endl;
                        gotLeftBlock = block_available(leftTaken, rightTaken, size, blockSize);
                        if (gotLeftBlock) {
                            posLeftBlock = leftTaken;
                            leftBlock = get_left_block(left, &leftTaken, blockSize);
                            // DEBUG
                            //std::cout << "thread "<< t_num << " took left block: " << posLeftBlock << std::endl;
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
    ForwardIt parallel_partition_phase_two(ForwardIt left, UnaryPredicate p, size_t size, size_t blockSize, size_t ln, size_t rn, std::vector<size_t> &remainingBlocks){

        bool gotBlocks;
        size_t leftFromRemaining, rightFromRemaining;
        ForwardIt leftBlock, rightBlock;

        // try to get left and right blocks and neutralize both
        gotBlocks = get_blocks_from_remaining(&remainingBlocks, ln, size-rn, &leftFromRemaining,
                                              &rightFromRemaining);
        while(gotBlocks) {
            //std::cout << "5" << std::endl;
            //std::cout << "leftBlockPos: " << blockPos << std::endl;
            leftBlock = std::next(left, leftFromRemaining);

            //std::cout << "rightBlockPos: " << blockPos << std::endl;
            rightBlock = std::next(left, rightFromRemaining);

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

        // now, neutralized but not well-placed blocks have to be swapped, if there are any

        size_t rightOfLN = ln, leftOfRN = size-rn-blockSize;
        //std::cout << "ln: " << ln << " rn: " << rn << std::endl;

        //try to get left block for swap
        gotBlocks = get_left_block_from_remaining(&remainingBlocks, ln, &leftFromRemaining);
        while(gotBlocks){
            //std::cout << "6" << std::endl;

            //get left block iterator
            //std::cout << "leftBlockPos: " << blockPos << std::endl;
            leftBlock = std::next(left, leftFromRemaining);

            //get neutralized block between LN and RN to swap with
            size_t pos = rightOfLN;
            // skip over all blocks that are still non-neutralized
            size_t shift = 0;  // store the number of non-neutralized elements right of LN in this shift, to use it later
            while(std::find(remainingBlocks.begin(), remainingBlocks.end(), pos+shift) != remainingBlocks.end()){
                //std::cout << "7" << std::endl;
                shift += blockSize;
            }
            pos = pos + shift;
            //std::cout << "left swap with pos: " << pos << std::endl;
            rightBlock = std::next(left, pos);

            //swap
            for(size_t i = 0; i< blockSize; i++){
                //std::cout << "leftBlock: " << *leftBlock << " rightBlock: " << *rightBlock << std::endl;
                std::swap(*leftBlock, *rightBlock);
                leftBlock = std::next(leftBlock);
                rightBlock = std::next(rightBlock);
            }
            // remember that you've already looked through everything in the shift
            // and the neutralized block found has been swapped, so it can be skipped in the next search too
            rightOfLN += shift+blockSize;
            gotBlocks = get_left_block_from_remaining(&remainingBlocks, ln, &leftFromRemaining);
        }

        //try to get right block for swap
        gotBlocks = get_right_block_from_remaining(&remainingBlocks, size-rn, &rightFromRemaining);
        while(gotBlocks) {
            //std::cout << "8" << std::endl;
            //get right block iterator
            //std::cout << "rightBlockPos: " << blockPos << std::endl;
            rightBlock = std::next(left, rightFromRemaining);

            //get neutralized block between LN and RN to swap with
            size_t pos = leftOfRN;
            // skip over all blocks that are still non-neutralized
            size_t shift = 0;
            while(std::find(remainingBlocks.begin(), remainingBlocks.end(), pos-shift) != remainingBlocks.end()){
                //std::cout << "9" << std::endl;
                shift += blockSize;
            }
            pos = pos-shift;
            //std::cout << "right swap with pos: " << pos << std::endl;
            leftBlock = std::next(left, pos);

            //swap
            for(size_t i = 0; i< blockSize; i++){
                //std::cout << "leftBlock: " << *leftBlock << " rightBlock: " << *rightBlock << std::endl;
                std::swap(*leftBlock, *rightBlock);
                leftBlock = std::next(leftBlock);
                rightBlock = std::next(rightBlock);
            }
            // remember that you've already looked through everything in the shift
            // and the neutralized block found has been swapped, so it can be skipped in the next search too
            leftOfRN -= shift+blockSize;
            gotBlocks = get_right_block_from_remaining(&remainingBlocks, size-rn, &rightFromRemaining);
        }

        //sort remaining values
        ForwardIt leftIt = std::next(left,ln);
        ForwardIt rightIt = std::next(left, size-rn);

        return std::partition(leftIt, rightIt, p);
    }

    template <typename ForwardIt, typename UnaryPredicate>
    ForwardIt partition(ForwardIt left, ForwardIt afterLast, UnaryPredicate p, int numThreads){
        //TODO maybe use long to make user there is no overflow on big arrays
        size_t size = std::distance(left, afterLast);
        size_t ln, rn;

        size_t blockSize = BLOCK_BYTES / sizeof(typename ForwardIt::value_type);

        //TODO maybe use list ?
        //TODO maybe only give pointers to functions
        std::vector<size_t> remainingBlocks;

        //print array
        std::cout << "all elements in Input: ";
        ForwardIt temp = left;
        while(temp != afterLast){
            std::cout << *temp << ' ';
            temp = std::next(temp);
        }
        std::cout << std::endl;


        remainingBlocks = parallel_partition_phase_one(left, p, numThreads, size, blockSize, &ln, &rn);

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

        auto split = parallel_partition_phase_two(left, p, size, blockSize, ln, rn, remainingBlocks);

        //print array
        std::cout << std::endl;
        std::cout << "all elements in Input: ";
        temp = left;
        while(temp != afterLast){
            std::cout << *temp << ' ';
            temp = std::next(temp);
        }

        return split;
    }

    template <typename ForwardIt, typename Compare>
    inline
    auto choose_pivot(ForwardIt begin, ForwardIt end, Compare c) {
        auto first = *begin;
        auto size = std::distance(begin, end);
        ForwardIt midIter = std::next(begin, size/2);
        auto middle = *midIter;
        auto last = *std::next(midIter, size - size/2 - 1);

        // DEBUG: these outputs illustrate how quickly a state of near sorted-ness is reached
        //std::cout << "FIRST : " << first << std::endl;
        //std::cout << "MIDDLE: " << middle << std::endl;
        //std::cout << "LAST  : " << last << std::endl;

        auto pivot = (std::min({first, middle, last}, c) + std::max({first, middle, last}, c)) / 2;

        //std::cout << "PIVOT  : " << pivot << std::endl;

        return pivot;
    }

    // source: https://stackoverflow.com/questions/18453945/c-generic-insertion-sort
    template<typename Iter, typename Compare>
    inline
    void insertion_sort(Iter first, Iter last, Compare c)
    {
        for (Iter it = first; it != last; ++it)
            std::rotate(std::upper_bound(first, it, *it, c), it, std::next(it));
    }

    template <typename ForwardIt, typename Compare>
    void helpingSort(ForwardIt left, ForwardIt afterLast, Compare c, unsigned long blockSize, std::stack<std::pair<ForwardIt, ForwardIt>> *sortStack, std::mutex *stackMutex, std::atomic<size_t> *completedElements, size_t totalElements) {
        size_t breakoff = blockSize * breakoffFactor / 64;
        start:
        auto size = std::distance(left, afterLast);
        if (size <= breakoff) {    // empirically found factor
            // stop recursion and simply sort sequentially
            seq:
            *completedElements += size;
            //insertion_sort(left, afterLast, c);
            std::sort(left, afterLast, c);
        } else {
            // recursively part the interval in two
            auto pivot = choose_pivot(left, afterLast, c);
            auto predicate = [pivot, c](const auto& em){ return c(em, pivot); };
            auto split = std::partition(left, afterLast, predicate);

            // if the following condition holds true, then the problem won't be fixed through recursion
            // and the remaining values have to be sorted sequentially (though they probably really are already)
            if (split == left || split == afterLast) {
                //std::cout << "test" << std::endl;
                goto seq;
            }


            // left isn't small so push it to the work stack
            //int leftSize = std::distance(left, split);
            //std::cout << "leftSize: " << leftSize << std::endl;
            stackMutex->lock();
            sortStack->push(std::pair<ForwardIt, ForwardIt>(left, split));
            stackMutex->unlock();

            // right isn't small so push it to the work stack
            left = split;
            goto start;


            /*  // AN OLDER ATTEMPT
            // push the first half onto the stack and recursively sort the latter
            // but only if you can; if the lock is taken start sorting first
            if (stackMutex->try_lock()) {
                sortStack->push(std::pair<ForwardIt, ForwardIt>(left, split));
                stackMutex->unlock();
                helpingSort(split, afterLast, c, blockSize, sortStack, stackMutex, completedThreads);
            } else {
                // lock is taken, sort first, then push
                helpingSort(split, afterLast, c, blockSize, sortStack, stackMutex, completedThreads);
                stackMutex->lock();
                sortStack->push(std::pair<ForwardIt, ForwardIt>(left, split));
                stackMutex->unlock();
            }
            */

            /*  // A SLIGHTLY YOUNGER ATTEMPT TO IMPROVE/FIX THE OLDER ATTEMPT
            if (!leftSmall) {
                if (stackMutex->try_lock()) {
                    sortStack->push(std::pair<ForwardIt, ForwardIt>(left, split));
                    stackMutex->unlock();
                    if (!rightSmall) {
                        helpingSort(split, afterLast, c, blockSize, sortStack, stackMutex, completedThreads);
                    } else {
                        // right is small, so sort sequentially
                        std::sort(split, afterLast, c);
                    }
                } else {
                    // lock is taken, sort first, then push
                    if (!rightSmall) {
                        helpingSort(split, afterLast, c, blockSize, sortStack, stackMutex, completedThreads);
                    } else {
                        // right is small, so sort sequentially
                        std::sort(split, afterLast, c);
                    }
                    stackMutex->lock();
                    sortStack->push(std::pair<ForwardIt, ForwardIt>(left, split));
                    stackMutex->unlock();
                }
            } else {
                // left is small, so sort sequentially
                std::sort(left, split, c);
                if (!rightSmall) {
                    helpingSort(split, afterLast, c, blockSize, sortStack, stackMutex, completedThreads);
                } else {
                    // right is small, so sort sequentially
                    std::sort(split, afterLast, c);
                }
            }
            */
        }
        // start helping by popping sorting jobs
        while(true) {
            //std::cout <<"in while"<<std::endl;
            // try to get a job
            if (stackMutex->try_lock()) {
                if (!sortStack->empty()) {
                    std::pair<ForwardIt, ForwardIt> sortJob = sortStack->top();
                    sortStack->pop();
                    //int sSize = sortStack->size();
                    //std::cout << "STACK SIZE: " << sSize << std::endl;
                    stackMutex->unlock();
                    left = sortJob.first;
                    afterLast = sortJob.second;
                    goto start;
                    // DEBUG: recursion depth might become a problem here (for some reason...)
                    //helpingSort(sortJob.first, sortJob.second, c, blockSize, sortStack, stackMutex, completedThreads);
                } else {
                    stackMutex->unlock();
                }
            }
            // check for total completion
            size_t completed = *completedElements;
            //std::cout << "COMPLETED: " << completed << std::endl;
            //std::cout << "TOTAL: " << totalElements << std::endl;
            if (completed >= totalElements) {
                break;
            }
        }
    }

    template <typename ForwardIt, typename Compare>
    void _quicksort(ForwardIt left, ForwardIt afterLast, Compare c, int numThreads, size_t blockSize, std::stack<std::pair<ForwardIt, ForwardIt>> *sortStack, std::mutex *stackMutex, std::atomic<size_t> *completedElements, size_t totalElements) {

        if (numThreads == 1) {
            // classic: just sort sequentially
            //std::sort(left, afterLast, c);
            // advanced: sort sequentially with helping scheme
            helpingSort(left, afterLast, c, blockSize, sortStack, stackMutex, completedElements, totalElements);
            return;
        }

        size_t size = std::distance(left, afterLast);
        size_t ln, rn;
        auto pivot = choose_pivot(left, afterLast, c);
        auto predicate = [pivot, c](const auto& em){ return c(em, pivot); };

        auto remaining = parallel_partition_phase_one(left, predicate, numThreads, size, blockSize, &ln, &rn);

        if (ln == 0) {
            // no blocks could be neutralized, default back to std::sort
            //std::sort(left, afterLast, c);
            helpingSort(left, afterLast, c, blockSize, sortStack, stackMutex, completedElements, totalElements);
            return;
        }

        ForwardIt split = parallel_partition_phase_two(left, predicate, size, blockSize, ln, rn, remaining);

        // split the threads based on the size of the new partitions
        auto s = std::distance(left, split);
        int processorSplit = (size_t)numThreads * s / size;
        if (processorSplit == 0)
            processorSplit = 1;


#pragma omp parallel num_threads(2)
        {
            if(omp_get_thread_num() == 0){
                _quicksort(left, split, c, processorSplit, blockSize, sortStack, stackMutex, completedElements, totalElements);
            }
            else{
                _quicksort(split, afterLast, c, numThreads - processorSplit, blockSize, sortStack, stackMutex, completedElements, totalElements);
            }
        }
    }

    template <typename ForwardIt, typename Compare>
    void quicksort(ForwardIt left, ForwardIt afterLast, Compare c) {
        omp_set_nested(256);
        int numThreads = omp_get_max_threads();
        size_t totalElements = std::distance(left, afterLast);
        size_t blockSize = BLOCK_BYTES / sizeof(typename ForwardIt::value_type);

        auto sortStack = std::stack<std::pair<ForwardIt, ForwardIt>>();
        std::mutex stackMutex;
        std::atomic<size_t> completedElements{0};

        _quicksort(left, afterLast, c, numThreads, blockSize, &sortStack, &stackMutex, &completedElements, totalElements);
    }

    template <typename ForwardIt, typename Compare>
    void nth_element(ForwardIt begin, ForwardIt n_th, ForwardIt end, Compare c, int numThreads) {
        auto pivot = choose_pivot(begin, end, c);
        auto predicate = [pivot, c](const auto& em){ return c(em, pivot); };
        auto size = std::distance(begin, end);
        size_t blockSize = BLOCK_BYTES / sizeof(typename ForwardIt::value_type);
        size_t ln, rn;

        const size_t BREAKOFF = 4;
        if (size < BREAKOFF) {
            // for small sizes just use std::nth_element
            std::nth_element(begin, n_th, end, c);
            return;
        }

        auto remaining = parallel_partition_phase_one(begin, predicate, numThreads, size, blockSize, &ln, &rn);

        ForwardIt split = parallel_partition_phase_two(begin, predicate, size, blockSize, ln, rn, remaining);

        // check whether you have to search left or right of the split
        // TODO: find out whether there's a better way to calculate this without making any assumptions about the iterators,
        //       or whether we should try to somehow check if the iterators are random access iterators so we can just do `n_th < split`
        auto n = std::distance(begin, n_th);
        auto s = std::distance(begin, split);
        if (n < s) {
            // search left
            nth_element(begin, n_th, split, c, numThreads);
        } else {
            // search right
            nth_element(split, n_th, end, c, numThreads);
        }
    }

}

#endif //PARALLEL_PARTITION_LIBRARY_H
