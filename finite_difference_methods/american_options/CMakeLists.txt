project(american_option)

set(exe1 main27)
set(sources1 bs_eq.cpp
            bs_model01.cpp
            explicit_lcp.cpp
            explicit_method.cpp 
            fd_method.cpp
            "option.cpp" 
            ${exe1}.cpp)
set(includes1 ./)

add_executable(${exe1} ${sources1})
target_include_directories(${exe1} PUBLIC ${includes1})