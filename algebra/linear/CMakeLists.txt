cmake_minimum_required(VERSION 3.24)

project(linear_algebra)

include_directories("../../")

set(SOURCE matrix.c)
set(HEADER matrix.h)

add_library(matrix ${SOURCE} ${HEADER})

target_include_directories(matrix PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCE vector.c)
set(HEADER vector.h)

add_library(vector ${SOURCE} ${HEADER})

target_include_directories(vector PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})