#ifndef PARALLEL_PARTITION_LIBRARY_H
#define PARALLEL_PARTITION_LIBRARY_H

#include <vector>

//returns 0 if Both neutralized, 1 if left ist neutralized, 2 if right is neutralized
template <typename ForwardIt, typename UnaryPredicate>
int neutralize(ForwardIt left, ForwardIt right, int blocksize, UnaryPredicate p){
    int i , j ;
    do {
        for (i = 0; i < blocksize; i++){
            left = std::next(left);
            if (!p(*left))
                break;
        }
        for (j = 0; j < blocksize; j++) {
            right = std::next(right);
            if (p(*right))
                break;
        }
        if ((i == blocksize) || (j == blocksize)){
            break;
        }

        std::swap(*left, *right);
        i++;
        j++;
    } while ( i < blocksize && j < blocksize );

    std::cout << i << " " << j << std::endl;
    if (i == blocksize) {
        if (j == blocksize)
            return 0;
        return 1;
    }
    return 2;
}



#endif //PARALLEL_PARTITION_LIBRARY_H
