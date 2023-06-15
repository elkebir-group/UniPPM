cmake_minimum_required(VERSION 3.21)
project(unippm)

set(CMAKE_CXX_STANDARD 14)
cmake_minimum_required(VERSION 3.21)
project(unippm)

set(CMAKE_CXX_STANDARD 20)
set(CXX 1)

set (CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/")

find_package(GUROBI REQUIRED)
MESSAGE(STATUS "Include dir of Gurobi: ${GUROBI_INCLUDE_DIRS}")
MESSAGE(STATUS "C library for Gurobi: ${GUROBI_LIBRARY}")
MESSAGE(STATUS "CXX library for Gurobi: ${GUROBI_CXX_DEBUG_LIBRARY}")
find_package(Boost)

include_directories(${BOOST_INCLUDE_DIRS} ${GUROBI_INCLUDE_DIRS})

set( CommonLibs m ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARY} )

link_libraries(${CommonLibs} Boost::boost )

add_executable(unippm src/unippm.cpp src/Input.h src/Input.cpp src/AncestryGraph.h src/AncestryGraph.cpp src/Hashing.cpp src/Hashing.h src/ILP.cpp src/ILP.h)