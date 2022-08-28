cmake_minimum_required(VERSION 3.24)

project(test)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")

include_directories("algebra/linear/")
include_directories("geometry/geodesic/")
include_directories("geometry/plane/")
include_directories("geometry/space/")

add_executable(${PROJECT_NAME} examples.c)

add_subdirectory(algebra/)
add_subdirectory(algebra/linear/)
add_subdirectory(geometry/geodesic/)
add_subdirectory(geometry/plane/)
add_subdirectory(geometry/space/)

target_link_libraries(${PROJECT_NAME} matrix)
target_link_libraries(${PROJECT_NAME} vector)
target_link_libraries(${PROJECT_NAME} complex)
target_link_libraries(${PROJECT_NAME} equation)
target_link_libraries(${PROJECT_NAME} gis)
target_link_libraries(${PROJECT_NAME} polygon)
target_link_libraries(${PROJECT_NAME} plane_curve_approximation)
target_link_libraries(${PROJECT_NAME} space_curve_approximation)