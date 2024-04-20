cmake_minimum_required(VERSION 3.22)
set(CMAKE_CXX_COMPILER "$ENV{CXX}")
message( "$ENV{CXX}")

project(Dinics)

set(MODE "-g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fopenmp ${MODE}")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_library(DinicsBase dinics-seq.cpp mgraph.cpp)
add_executable(main main.cpp)
target_link_libraries(main PUBLIC DinicsBase)

