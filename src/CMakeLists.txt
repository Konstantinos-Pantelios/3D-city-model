cmake_minimum_required(VERSION 3.17)
project(3DCityModel)

set(CMAKE_CXX_STANDARD 17)

include_directories(.)

add_executable(src
        main.cpp
        Point.h
        json.hpp
        kdtree.h
)

