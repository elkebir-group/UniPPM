cmake_minimum_required(VERSION 3.21 FATAL_ERROR)
project(unippm)

set(CMAKE_CXX_STANDARD 11)

find_package(cryptominisat5 REQUIRED)
MESSAGE(STATUS "Include dir of cryptominisat5: ${CRYPTOMINISAT5_INCLUDE_DIRS}")
find_package(approxmc REQUIRED)
find_package(unigen REQUIRED)
find_package(Boost)

include_directories(${CRYPTOMINISAT5_INCLUDE_DIRS} ${APPROXMC_INCLUDE_DIRS} ${UNIGEN_INCLUDE_DIRS} ${BOOST_INCLUDE_DIRS})
set( CommonLibs m ${CRYPTOMINISAT5_LIBRARIES} ${APPROXMC_LIBRARIES} ${UNIGEN_LIBRARIES} )

#add_executable(unippm main.cpp Input.cpp Input.h AncestryGraph.cpp AncestryGraph.h Solver.cpp Solver.h CNF.cpp CNF.h Likelihood.cpp Likelihood.h)
add_executable(test test.cpp Input.cpp Input.h AncestryGraph.cpp AncestryGraph.h Solver.cpp Solver.h CNF.cpp CNF.h)

target_link_libraries(test ${CommonLibs} Boost::boost)


