


cmake_minimum_required(VERSION 3.10.0)

project(SimpleSundialSolver VERSION 0.1.0 LANGUAGES C CXX)

find_package(SUNDIALS REQUIRED)

add_subdirectory(src)
