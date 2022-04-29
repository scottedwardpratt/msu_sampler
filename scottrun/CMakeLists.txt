cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
set(GITHOME_MSU $ENV{GITHOME_MSU})

project(samplertest CXX)
project(samplertest_hydro CXX)
project(rhocount CXX)
project(sheartest CXX)
project(bulktest CXX)

add_custom_target(extern_msu_sampler COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/msu_sampler/software)
add_custom_target(extern_msu_commonutils COMMAND make WORKING_DIRECTORY ${GITHOME_MSU}/commonutils/software)

#include(${CMAKE_CURRENT_SOURCE_DIR}/../software/CMakeLists.txt)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Wold-style-cast -Werror=switch -std=c++11 -mfpmath=sse")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL 1.15 REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${GITHOME_MSU}/commonutils/software/include ${GITHOME_MSU}/msu_sampler/software/include)

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../software/lib ${GITHOME_MSU}/commonutils/software/lib)

set(SOURCE_FILES ./samplertest.cc)
SET_SOURCE_FILES_PROPERTIES(
      samplertest.cc PROPERTIES OBJECT_DEPENDS
    ../software/lib/libmsu_sampler.a
 )
add_executable(samplertest ${SOURCE_FILES})
 
set(SOURCE_FILES ./samplertest_hydro.cc)
SET_SOURCE_FILES_PROPERTIES(
      samplertest_hydro.cc PROPERTIES OBJECT_DEPENDS
    ../software/lib/libmsu_sampler.a
 )
add_executable(samplertest_hydro ${SOURCE_FILES})

set(SOURCE_FILES ./rhocount.cc)
SET_SOURCE_FILES_PROPERTIES(
      rhocount.cc PROPERTIES OBJECT_DEPENDS
    ../software/lib/libmsu_sampler.a
 )
add_executable(rhocount ${SOURCE_FILES})

set(SOURCE_FILES ./bulktest.cc)
 SET_SOURCE_FILES_PROPERTIES(
	bulktest.cc PROPERTIES OBJECT_DEPENDS ../software/lib/libmsu_sampler.a
)
add_executable(bulktest ${SOURCE_FILES})

set(SOURCE_FILES ./sheartest.cc)
SET_SOURCE_FILES_PROPERTIES(
	sheartest.cc PROPERTIES OBJECT_DEPENDS ../software/lib/libmsu_sampler.a
)
add_executable(sheartest ${SOURCE_FILES})



add_dependencies(samplertest extern_msu_sampler)
add_dependencies(samplertest_hydro extern_msu_sampler)
add_dependencies(rhocount extern_msu_sampler)
add_dependencies(sheartest extern_msu_sampler)
add_dependencies(bulktest extern_msu_sampler)
find_library(${GITHOME_MSU}/msu_sampler/software/lib msu_sampler)
find_library(${GITHOME_MSU}/commonutils/software/lib msu_commonutils)
target_link_libraries(samplertest msu_sampler msu_commonutils ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})
target_link_libraries(samplertest_hydro msu_sampler msu_commonutils ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})
target_link_libraries(rhocount msu_sampler ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})
target_link_libraries(sheartest msu_sampler ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})
target_link_libraries(bulktest msu_sampler ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${EIGEN3_LIBRARY})
