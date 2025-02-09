project(summing_up)

set(exe1 main11)

set(sources1 bin_model02.cpp options06.cpp double_digital_option.cpp ${exe1}.cpp)
set(includes1 ./)

add_executable(${exe1} ${sources1})
target_include_directories(${exe1} PUBLIC ${includes1})