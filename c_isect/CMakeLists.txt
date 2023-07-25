cmake_minimum_required(VERSION 3.25)
project(c_isect)

set(CMAKE_CXX_STANDARD 23)
add_executable(c_isect main.cpp csv.h geometry.cpp geometry.h)
