cmake_minimum_required(VERSION 3.27)

project(Dinics)

set(MODE "-g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 ${MODE}")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_executable(main dinics-seq.cpp)
