cmake_minimum_required(VERSION 3.27)
project(ND_SpinMC)

set(CMAKE_CXX_STANDARD 17)

include_directories(include)

add_executable(ND_SpinMC
        src/main.cpp
        src/arguments.cpp
        src/metropolis.cpp
        src/stats.cpp)
