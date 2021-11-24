#include "library.h"

#include <iostream>
#include <vector>

void hello() {
    std::cout << "Hello, World!" << std::endl;
}

int neutralize(std::vector<float>::iterator left,std::vector<float>::iterator right, int blocksize, float pivot){


    int i , j ;
    do {
        for (i = 0; i < blocksize; i++){
            if (*(left + i) > pivot) {
                break;
            }
        }
        for (j = 0; j < blocksize; j++) {
            if (*(right - j) < pivot) {
                break;
            }
        }
        if ((i == blocksize) || (j == blocksize)){
            break;
        }

        //maybe use standard function for swap
        float temp = *(left+i);
        *(left+i) = *(right-j);
        *(right-j) = temp;

        i++;
        j++;
    }while ( i < blocksize && j < blocksize );
    std::cout << i << " " << j << std::endl;
    if ( i == blocksize && j == blocksize )
        return 0;
    if ( i == blocksize)
        return 1;
    return 2;



}
