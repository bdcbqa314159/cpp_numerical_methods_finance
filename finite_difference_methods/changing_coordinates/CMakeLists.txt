project(changing_coordinates)

set(exe1 main26)
set(sources1 bs_eq.cpp
            heat_eq.cpp
            bs_model01.cpp 
            implicit_scheme.cpp 
            fd_method.cpp
            matrix.cpp
            "option.cpp" 
            ${exe1}.cpp)
set(includes1 ./)

add_executable(${exe1} ${sources1})
target_include_directories(${exe1} PUBLIC ${includes1})