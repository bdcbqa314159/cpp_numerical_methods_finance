cmake_minimum_required(VERSION 3.30)

project(cpp_numerical_methods_finance LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(GSL REQUIRED)
find_package(Boost REQUIRED COMPONENTS date_time)

add_subdirectory(numerical_methods)

add_subdirectory(chapter_01)

