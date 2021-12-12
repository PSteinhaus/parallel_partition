cmake_minimum_required(VERSION 3.20)
project(parallel_partition LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

find_package(OpenMP REQUIRED)

add_library(parallel_partition STATIC include/library.h src/library.cpp)

target_include_directories(parallel_partition PUBLIC include PRIVATE src)
target_link_libraries(parallel_partition PUBLIC OpenMP::OpenMP_CXX)