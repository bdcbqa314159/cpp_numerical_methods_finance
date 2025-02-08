project(numerical_methods)

set(includes includes/)
set(sources src/message.cpp
            src/binomial_pricing.cpp)

add_library(${PROJECT_NAME} ${sources})
target_include_directories(${PROJECT_NAME} PUBLIC ${includes})
target_link_libraries(${PROJECT_NAME} PUBLIC GSL::gsl GSL::gslcblas Boost::boost Boost::date_time)
