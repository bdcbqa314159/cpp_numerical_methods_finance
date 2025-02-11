project(computing_implied_volatility)

set(exe1 main18)

set(sources1 eur_call.cpp ${exe1}.cpp)
set(includes1 ./)

add_executable(${exe1} ${sources1})
target_include_directories(${exe1} PUBLIC ${includes1})
