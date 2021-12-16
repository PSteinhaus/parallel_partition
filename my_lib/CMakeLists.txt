cmake_minimum_required(VERSION 3.20)
project(parallel_partition LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
find_package(OpenMP REQUIRED)


#compile library
add_library(parallel_partition STATIC include/library.h src/library.cpp)

#opemp for library
target_include_directories(parallel_partition PUBLIC include PRIVATE src)
target_link_libraries(parallel_partition PUBLIC OpenMP::OpenMP_CXX)



#tests with catch
add_subdirectory(catch)

add_executable(library_tests catch_main.cpp catch_tests.cpp)

target_link_libraries(library_tests PRIVATE parallel_partition PRIVATE catch)