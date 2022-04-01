cmake_minimum_required(VERSION 3.21)
project(unippm)

set(CMAKE_CXX_STANDARD 11)

set (CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/")

#include(FetchContent)

find_package(GUROBI REQUIRED)
MESSAGE(STATUS "Include dir of Gurobi: ${GUROBI_INCLUDE_DIRS}")
find_package(cryptominisat5 REQUIRED)
MESSAGE(STATUS "Include dir of cryptominisat5: ${CRYPTOMINISAT5_INCLUDE_DIRS}")
find_package(approxmc REQUIRED)
MESSAGE(STATUS "Include dir of approxmc: ${APPROXMC_INCLUDE_DIRS}")
find_package(unigen REQUIRED)
MESSAGE(STATUS "Include dir of unigen: ${UNIGEN_INCLUDE_DIRS}")
find_package(Boost)

#set(BUILD_DEPS ON)
#FetchContent_Declare(
#        or-tools
#        GIT_REPOSITORY https://github.com/google/or-tools.git
#        GIT_TAG        v9.3
#)
#
#FetchContent_MakeAvailable(or-tools)
#find_package(GUROBI REQUIRED)

#include_directories(${CRYPTOMINISAT5_INCLUDE_DIRS} ${APPROXMC_INCLUDE_DIRS} ${UNIGEN_INCLUDE_DIRS} ${BOOST_INCLUDE_DIRS} ${GUROBI_INCLUDE_DIRS})
include_directories(${CRYPTOMINISAT5_INCLUDE_DIRS} ${APPROXMC_INCLUDE_DIRS} ${UNIGEN_INCLUDE_DIRS} ${BOOST_INCLUDE_DIRS} ${Clp_INCLUDE_DIRS})

set( CommonLibs m ${CRYPTOMINISAT5_LIBRARIES} ${APPROXMC_LIBRARIES} ${UNIGEN_LIBRARIES} ${GUROBI_LIBRARY} )

link_libraries(${CommonLibs} Boost::boost )

add_executable(unippm-reads main.cpp Input.cpp Input.h AncestryGraph.cpp AncestryGraph.h Solver.cpp Solver.h CNF.cpp CNF.h Likelihood.cpp Likelihood.h)# MultithreadSampler.cpp MultithreadSampler.h) #Parallel.cpp Parallel.h)

add_executable(unippm-freqs main-freqs.cpp Input.cpp Input.h AncestryGraph.cpp AncestryGraph.h Solver.cpp Solver.h CNF.cpp CNF.h) #Parallel.cpp Parallel.h)

add_executable(enumerate-likelihood main-enumllh.cpp Input.cpp Input.h AncestryGraph.cpp AncestryGraph.h Solver.cpp Solver.h CNF.cpp CNF.h Likelihood.cpp Likelihood.h)

