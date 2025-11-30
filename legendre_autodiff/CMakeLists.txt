cmake_minimum_required(VERSION 3.10)
project(LegendreAutoDiff LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Create the executable
add_executable(legendre_autodiff main.cpp)

# (Optional) Make it faster
set(CMAKE_CXX_FLAGS_RELEASE "-O3")