project(path_dependent_basket_options)

set(exe1 main23)
set(sources1 matrix.cpp bs_model02.cpp path_dep_option05.cpp ${exe1}.cpp)
set(includes1 ./)

add_executable(${exe1} ${sources1})
target_include_directories(${exe1} PUBLIC ${includes1})
