cmake_minimum_required(VERSION 3.20)
project(ASC_ode)

set (CMAKE_CXX_STANDARD 20)

add_compile_options(-fpermissive)

include_directories(src nanoblas/src)

add_subdirectory (src)
add_subdirectory (nanoblas)
add_subdirectory (Exercises)
add_subdirectory (mechsystem)

file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/Exercises/data")

add_executable (test_ode demos/test_ode.cpp)
target_link_libraries (test_ode PUBLIC nanoblas)

add_executable (demo_autodiff demos/demo_autodiff.cpp)
target_link_libraries (demo_autodiff PUBLIC nanoblas)

