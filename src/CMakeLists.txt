CMAKE_MINIMUM_REQUIRED ( VERSION 2.8.12 )

SET ( LMP_PATH /path/to/lammps )
SET ( RANDSPG_PATH /path/to/randspg )
PROJECT ( "SPINNER" )
SET ( PROJECT_VERSION_MAJOR 1 )
SET ( PROJECT_VERSION_MINOR 1 )

find_package( MPI REQUIRED )
if ( ${MPI_CXX_VERSION_MAJOR} LESS 3 )
    message( FATAL_ERROR "MPI-3 is required." )
endif()

SET ( CMAKE_VERBOSE_MAKEFILE true )
SET ( SRC_FILES 
        main.cpp
        calculator.cpp
        config.cpp
        input.cpp
        utils.cpp )
SET ( CMAKE_CXX_COMPILER "mpicxx" )

ADD_COMPILE_OPTIONS ( -O3 -g -Wall -std=c++11 )
LINK_LIBRARIES( lammps RandSpgLib m )
LINK_DIRECTORIES( ${LMP_PATH}/src ${RANDSPG_PATH}/build )
INCLUDE_DIRECTORIES ( ${LMP_PATH}/src ${RANDSPG_PATH}/include )

ADD_EXECUTABLE ( ${CMAKE_PROJECT_NAME} ${SRC_FILES} )
