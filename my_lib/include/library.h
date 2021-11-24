#ifndef PARALLEL_PARTITION_LIBRARY_H
#define PARALLEL_PARTITION_LIBRARY_H

#include <vector>

void hello();

//returns 0 if Both neutralized, 1 if left ist neutralized, 2 if right is neutralized
int neutralize(std::vector<float>::iterator left,std::vector<float>::iterator right,int blocklength, float pivot);

#endif //PARALLEL_PARTITION_LIBRARY_H
