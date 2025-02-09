project(virtual_inheritance)

set(exe1 main12)

set(sources1 bin_model02.cpp options08.cpp ${exe1}.cpp)
set(includes1 ./)

add_executable(${exe1} ${sources1})
target_include_directories(${exe1} PUBLIC ${includes1})