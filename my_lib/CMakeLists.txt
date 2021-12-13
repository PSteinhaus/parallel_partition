cmake_minimum_required(VERSION 3.20)
project(parallel_partition LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

find_package(OpenMP REQUIRED)

add_subdirectory(catch)

add_library(parallel_partition STATIC include/library.h src/library.cpp)

add_executable(library_tests catch_main.cpp catch_tests.cpp )

target_include_directories(parallel_partition PUBLIC include PRIVATE src)
target_link_libraries(parallel_partition PUBLIC OpenMP::OpenMP_CXX)

#catch to seems to need at least O1 or higher
#(Error: too many sections (32835))
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")

target_link_libraries(library_tests PRIVATE catch)