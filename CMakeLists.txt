cmake_minimum_required(VERSION 3.12.0)
project(navier_stokes LANGUAGES CXX C)

include(cmake-common.cmake)

add_executable(navier_stokes src/main.cpp src/NavierStokes.cpp)
deal_ii_setup_target(navier_stokes)
