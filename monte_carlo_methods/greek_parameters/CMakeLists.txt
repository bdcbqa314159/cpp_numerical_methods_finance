project(greek_parameters)

set(exe1 main21)
set(sources1 bs_model01.cpp path_dep_option03.cpp ${exe1}.cpp)
set(includes1 ./)

add_executable(${exe1} ${sources1})
target_include_directories(${exe1} PUBLIC ${includes1})
