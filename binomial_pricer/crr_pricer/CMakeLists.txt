project(crr_pricer)

set(exe1 main05)

set(sources1 bin_model01.cpp options01.cpp ${exe1}.cpp)
set(includes1 ./)

add_executable(${exe1} ${sources1})
target_include_directories(${exe1} PUBLIC ${includes1})