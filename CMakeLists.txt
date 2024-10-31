cmake_minimum_required(VERSION 3.27)
project(NM_3)

set(CMAKE_CXX_STANDARD 17)

add_executable(NM_3 main.cpp
        spline/smoothing_spline.cpp
        spline/point/point.cpp
)
