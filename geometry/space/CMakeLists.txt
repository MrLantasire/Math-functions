cmake_minimum_required(VERSION 3.24)

project(space)

include_directories("../../" "../../algebra/" "../../algebra/linear/")

set(SOURCE space_curve_approximation.c)
set(HEADER space_curve_approximation.h)
add_library(space_curve_approximation ${SOURCE} ${HEADER})

target_link_libraries(plane_curve_approximation complex)
target_link_libraries(plane_curve_approximation equation)

target_include_directories(space_curve_approximation PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
