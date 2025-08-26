cmake_minimum_required(VERSION 3.10)
project(voxelsolve)

set(CMAKE_C_STANDARD 11)

#sdl 
#add_subdirectory(external/sdl)
#include_directories(external/sdl/include)


add_executable(run main.c)

target_link_libraries(voxelsolve m)
