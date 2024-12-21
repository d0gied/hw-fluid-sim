cmake_minimum_required(VERSION 3.28)

# Set the project name
project(Fluid)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add the executable
add_executable(main main.cpp fixed.h simulation.h)
