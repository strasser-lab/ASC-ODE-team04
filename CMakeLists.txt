cmake_minimum_required(VERSION 3.20)
project(ASC_ode)

set (CMAKE_CXX_STANDARD 20)

add_compile_options(-fpermissive)

include_directories(src nanoblas/src)

add_subdirectory (src)
add_subdirectory (nanoblas)

add_subdirectory (mechsystem)

add_executable (test_ode demos/test_ode.cpp)
target_link_libraries (test_ode PUBLIC nanoblas)

add_executable (demo_autodiff demos/demo_autodiff.cpp)
target_link_libraries (demo_autodiff PUBLIC nanoblas)

add_executable (Exercise_17_4_1_mass_spring demos/Exercise_17_4_1_mass_spring.cpp)
target_link_libraries (Exercise_17_4_1_mass_spring PUBLIC nanoblas)

add_executable (Exercise_17_4_1_RCCircuit demos/Exercise_17_4_1_RCCircuit.cpp)
target_link_libraries (Exercise_17_4_1_RCCircuit PUBLIC nanoblas)
