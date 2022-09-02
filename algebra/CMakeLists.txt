cmake_minimum_required(VERSION 3.24)

project(algebra)

include_directories("../")

set(SOURCE complex.c)
set(HEADER complex.h)
add_library(complex ${SOURCE} ${HEADER})

set(SOURCE equation.c)
set(HEADER equation.h)
add_library(equation ${SOURCE} ${HEADER})
target_link_libraries(equation complex)

set(SOURCE quaternion.c)
set(HEADER quaternion.h)
add_library(quaternion ${SOURCE} ${HEADER})
target_link_libraries(quaternion vector)

target_include_directories(complex PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(equation PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(quaternion PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
