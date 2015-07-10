#  CMakeLists.txt
#  Created on: April 13, 2013
#      Author: cdremy@umich.edu
#
###############################################################
# BASIC CONFIGURATION
CMAKE_MINIMUM_REQUIRED (VERSION 2.6)  # Check the version type

PROJECT( "RAMone" CXX )       # Define the project name

MESSAGE("MESSAGE: CMakeList File for Series and Parallel Elastic Actuation")
MESSAGE(${PROJECT_SOURCE_DIR})

LIST( APPEND CMAKE_MODULE_PATH  #Appends the two paths to the (empty) list of module paths
	${PROJECT_SOURCE_DIR}/CMake
	${PROJECT_SOURCE_DIR}/CMake/Modules
)
# These paths are used in any include and find_package statement
MESSAGE(${CMAKE_MODULE_PATH})   #There are only 2 entries

MESSAGE("Before include")
INCLUDE( SimOptDefaults )       # Process the listfile given in "SimOptDefaults"
# It is located in ${PROJECT_SOURCE_DIR}/CMake/Modules
MESSAGE("After include")

###############################################################
# GET USED PACKAGES
FIND_PACKAGE( MUSCOD REQUIRED )  # find and include "MUSCODConfig.cmake"
IF( MUSCOD_FOUND )
	INCLUDE( ${MUSCOD_USE_FILE} )
ENDIF( MUSCOD_FOUND )


###############################################################
# INCLUDES FOR THE LIBRARY
# Our own stuff:
INCLUDE_DIRECTORIES(SRC)
# Eigen for Matrix computations:
SET( EIGEN3_INCLUDE_DIR /usr/include/eigen3)
#FIND_PACKAGE(Eigen3 REQUIRED)
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})

###############################################################
# BUILD LIBRARIES OF THE OPTIMAL CONTROL PROBLEM
# NOTE: Here you need to add all source files that are required for your
#       optimal control problem.

# Toy Problem
ADD_LIBRARY ( RAMone_Walking
    AB_Walking.cpp
	SRC/AB_Constraints.cpp
	SRC/AB_Dynamics.cpp
	)

###############################################################
# BUILD THE EXECUTABLES FOR THE VELOCITY STUDIES
#ADD_EXECUTABLE ( PB_DDH_Processing
#	SRC/PB_DDH_Processing.cpp
#	)


###############################################################
#  SELF-TESTING
INCLUDE( TESTDartConfig ) # NOTE: line must be called BEFORE INCLUDE( Dart )!
INCLUDE( Dart )       # Makes CTest and Dart available, including call of ENABLE_TESTING()
INCLUDE( TestSuite )   # Defines the actual tests.
