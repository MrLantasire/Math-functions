cmake_minimum_required(VERSION 3.24)

project(plane)

include_directories("../../" "../../algebra/" "../../algebra/linear/")

set(SOURCE polygon.c)
set(HEADER polygon.h)
add_library(polygon ${SOURCE} ${HEADER})
target_include_directories(polygon PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCE plane_curve_approximation.c)
set(HEADER plane_curve_approximation.h)
add_library(plane_curve_approximation ${SOURCE} ${HEADER})
target_link_libraries(plane_curve_approximation complex)
target_link_libraries(plane_curve_approximation equation)
target_include_directories(plane_curve_approximation PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
