cmake_minimum_required(VERSION 3.20)
project(parallel_partition)

set(CMAKE_CXX_STANDARD 14)

add_library(parallel_partition STATIC include/library.h src/library.cpp)
target_include_directories(parallel_partition PUBLIC include PRIVATE src)

