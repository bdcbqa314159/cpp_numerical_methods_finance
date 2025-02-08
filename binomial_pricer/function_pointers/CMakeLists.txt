project(function_pointers)

set(exe1 main07)

set(sources1 bin_model01.cpp options03.cpp ${exe1}.cpp)
set(includes1 ./)

add_executable(${exe1} ${sources1})
target_include_directories(${exe1} PUBLIC ${includes1})