cmake_minimum_required(VERSION 3.10)
project(voxelsolve C)

set(CMAKE_C_STANDARD 11)

add_executable(voxelsolve main.c)

target_link_libraries(voxelsolve m)
