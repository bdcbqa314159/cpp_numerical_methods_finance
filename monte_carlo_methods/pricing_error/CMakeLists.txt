project(pricing_error)

set(exe1 main20)
set(sources1 bs_model01.cpp path_dep_option02.cpp ${exe1}.cpp)
set(includes1 ./)

add_executable(${exe1} ${sources1})
target_include_directories(${exe1} PUBLIC ${includes1})
