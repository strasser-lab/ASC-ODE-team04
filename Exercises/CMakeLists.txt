cmake_minimum_required(VERSION 3.10)
project(ASC_ode)

set(CMAKE_CXX_STANDARD 20)

add_executable (Exercise_17_4_1_mass_spring Exercise_17_4_1_mass_spring.cpp)
target_link_libraries (Exercise_17_4_1_mass_spring PUBLIC nanoblas)

add_executable (Exercise_17_4_1_RCCircuit Exercise_17_4_1_RCCircuit.cpp)
target_link_libraries (Exercise_17_4_1_RCCircuit PUBLIC nanoblas)

add_executable (PendulumAD_18_5 PendulumAD_18_5.cpp)
target_link_libraries (PendulumAD_18_5 PUBLIC nanoblas)