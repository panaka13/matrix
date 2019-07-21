cmake_minimum_required (VERSION 2.6)
project (matrix)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread --std=c++17")
add_subdirectory("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_SOURCE_DIR}/include")

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/exec")
add_executable(MatrixTest "${PROJECT_SOURCE_DIR}/src/test.cpp") 
