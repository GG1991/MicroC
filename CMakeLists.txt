
cmake_minimum_required(VERSION 3.0)

project(MicroC VERSION 0.1 LANGUAGES C)

find_library(
	PETSC_LIB petsc
	PATHS ENV PETSC_DIR
	PATH_SUFFIXES "$ENV{PETSC_ARCH}/lib"
	)

find_path(
	PETSC_INC_A petsc.h
	PATHS ENV PETSC_DIR
	PATH_SUFFIXES "include"
	)

find_path(
	PETSC_INC_B petscconf.h
	PATHS ENV PETSC_DIR
	PATH_SUFFIXES "$ENV{PETSC_ARCH}/include"
	)

set(PETSC_INC ${PETSC_INC_A} ${PETSC_INC_B})

option(TIMER "Enable time benchmarks instrumentation" OFF)
option(CGDEBUG "Enable CG algorithm debug" OFF)
option(FILTER "Enable Filtering on Average Stress and Ctan" OFF)

if (TIMER)
	add_definitions(-DTIMER)

	find_path(INSTRUMENT_INC instrument.h PATHS ~/GIT/Instrument PATH_SUFFIXES "include")
	find_library(INSTRUMENT_LIB instrument PATHS ~/GIT/Instrument PATH_SUFFIXES "build")
endif()

if (CGDEBUG)
	add_definitions(-DCGDEBUG)
endif()

if (FILTER)
	add_definitions(-DFILTER)
endif()

# Include Directories (for all targets)
include_directories(include ${CMAKE_BINARY_DIR} ${MICROC_INC} ${PETSC_INC} ${INSTRUMENT_INC})

# wildcard all the sources in src
file(GLOB SOURCESLIB src/*.c)

add_library(microc ${SOURCESLIB})
target_link_libraries(microc ${PETSC_LIB} ${INSTRUMENT_LIB} m)
set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -fopenmp")

enable_testing()
add_subdirectory(test)
