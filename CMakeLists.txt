cmake_minimum_required(VERSION 3.12.0)
project(navier_stokes LANGUAGES CXX C)

include(cmake-common.cmake)

add_executable(navier_stokes3D src/main3D.cpp src/NavierStokes3D.cpp)
deal_ii_setup_target(navier_stokes3D)

add_executable(navier_stokes2D src/main2D.cpp src/NavierStokes2D.cpp)
deal_ii_setup_target(navier_stokes2D)
