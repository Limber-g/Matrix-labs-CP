cmake_minimum_required(VERSION 3.20)
project(Labs_1_4_second_module)

set(CMAKE_CXX_STANDARD 17)

add_library(Matrix SHARED Matrix.h Matrix.cpp)
add_executable(Labs_1_4_second_module main.cpp)
target_link_libraries(Labs_1_4_second_module Matrix)
