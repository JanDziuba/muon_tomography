cmake_minimum_required(VERSION 3.25)
project(c_poca)

set(CMAKE_CXX_STANDARD 23)

add_executable(c_poca main.cpp csv.h Geometry.cpp Geometry.h)
