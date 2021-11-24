cmake_minimum_required(VERSION 3.20)
project(parallel_partition)

set(CMAKE_CXX_STANDARD 14)

add_subdirectory(my_lib)
add_executable(main main.cpp)
target_link_libraries(main PUBLIC parallel_partition)


