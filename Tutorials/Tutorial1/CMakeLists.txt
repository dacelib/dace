# DACE Examples
#
# This file shows how to use the DACE in a cmake project
# Simply import the DACE cmake package. It automatically provides the correct
# locations for headers and libraries.
# If the DACE is installed in a non-default location, tell cmake where to look:
#   set(dace_DIR "/my/custom/root/usr/local/lib/cmake/dace")
# It exposes two targets:
#   dace::dace   (the dynamic DACE library)
#   dace::dace_s (the static DACE library)
# Make sure to link your executable with one of those and include the correct
# header (<dace/dace.h> or <dace/dace_s.h>) in the code.

cmake_minimum_required (VERSION 2.8.4)

project(Examples1 CXX)

find_package(dace 2.0.0 REQUIRED)

add_executable(Example1 Example1.cpp)
target_link_libraries(Example1 PUBLIC dace::dace)

add_executable(Example2 Example2.cpp)
target_link_libraries(Example2 PUBLIC dace::dace)

add_executable(Example3 Example3.cpp)
target_link_libraries(Example3 PUBLIC dace::dace)

add_executable(Example4 Example4.cpp)
target_link_libraries(Example4 PUBLIC dace::dace)

add_executable(Example5 Example5.cpp)
target_link_libraries(Example5 PUBLIC dace::dace)

add_executable(Example6 Example6.cpp)
target_link_libraries(Example6 PUBLIC dace::dace)

add_executable(Example7 Example7.cpp)
target_link_libraries(Example7 PUBLIC dace::dace)

add_executable(Example8 Example8.cpp)
target_link_libraries(Example8 PUBLIC dace::dace)

add_executable(Example9 Example9.cpp)
target_link_libraries(Example9 PUBLIC dace::dace)

add_executable(Example10 Example10.cpp)
target_link_libraries(Example10 PUBLIC dace::dace)

add_executable(Example11Pointwise Example11Pointwise.cpp)
target_link_libraries(Example11Pointwise PUBLIC dace::dace)

add_executable(Example11DA Example11DA.cpp)
target_link_libraries(Example11DA PUBLIC dace::dace)

add_executable(Example12 Example12.cpp)
target_link_libraries(Example12 PUBLIC dace::dace)
