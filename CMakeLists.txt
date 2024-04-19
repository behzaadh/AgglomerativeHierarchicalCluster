cmake_minimum_required(VERSION 3.2)

project(AgglomerativeHierarchicalCluster LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Add the external project directory
SET(EIGEN3_INCLUDE_DIR "YOUR EIGEN PATH") # Here you must add the path of EIGEN3_INCLUDE_DIR link
SET(CMAKE_MODULE_PATH "${EIGEN3_INCLUDE_DIR}/cmake")
SET(EIGEN3_DIR "${EIGEN3_INCLUDE_DIR}/cmake")

# include_directories is needed for the compiler to know where looking for Eigen3 header files to be included
include_directories( ${EIGEN3_INCLUDE_DIR} )

add_executable(AgglomerativeHierarchicalCluster main.cpp
    Measures.h
    AgglomerativeHierarchicalCluster.h)
