cmake_minimum_required(VERSION 3.30)

project(cpp_numerical_methods_finance LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

MESSAGE(STATUS "the source directory is: ${CMAKE_SOURCE_DIR}")

find_package(GSL REQUIRED)
find_package(Boost REQUIRED COMPONENTS date_time)

add_subdirectory(numerical_methods)

add_subdirectory(binomial_pricer)
add_subdirectory(binomial_pricer_revisited)
add_subdirectory(american_options)
add_subdirectory(non_linear_solvers)

add_subdirectory(test)

