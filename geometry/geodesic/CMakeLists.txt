cmake_minimum_required(VERSION 3.24)

project(geodesic)

include_directories("../../")

set(SOURCE gis.c)
set(HEADER gis.h)

add_library(gis ${SOURCE} ${HEADER})

target_include_directories(gis PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
