project(testing_the_lib)

set(exe1 test)

add_executable(${exe1} test.cpp)
target_include_directories(${exe1} PUBLIC ${CMAKE_SOURCE_DIR}/numerical_methods/)
target_link_libraries(${exe1} numerical_methods)