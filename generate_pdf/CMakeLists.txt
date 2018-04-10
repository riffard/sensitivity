###############  CMakeLists for MimacTools build   #####################
#  Written by Quentin Riffard (riffard@lpsc.in2p3.fr)                  #
#  To create an Xcode project:                                         #
#         cmake -G Xcode .                                             #

###############             cmake version          #####################
cmake_minimum_required(VERSION 2.6)

###############              Project name          #####################
project(SoftAnalysis)

###############        Root include/lib dir        #####################
find_package(ROOT COMPONENTS RooStats MathCore)
IF(NOT ${ROOT_FOUND})
	message("Root not found... Let's try an other way ;)")
	include($ENV{ROOTSYS}/etc/root/cmake/FindROOT.cmake)
ENDIF()

include_directories(${ROOT_INCLUDE_DIRS})

link_directories( ${ROOT_LIBRARY_DIR})
########################################################################

###############        NEST include/lib dir        #####################
find_package(NEST REQUIRED)
include_directories(${NEST_INCLUDE_DIRS}/include)
########################################################################


###############     dm_calc include/lib dir        #####################
find_package(dm_calc REQUIRED)
include_directories("${dm_calc_INCLUDE_DIRS}")

########################################################################


########################################################################
if(NOT DEFINED CMAKE_MACOSX_RPATH)
  set(CMAKE_MACOSX_RPATH 0)
endif()
########################################################################

###############            Minimal flags           #####################
IF(${UNIX})
  set (CMAKE_CXX_FLAGS "-std=c++11 -g3 -Wall -O0")
ENDIF(${UNIX})

IF(${APPLE})
  set (CMAKE_CXX_FLAGS "-std=c++0x  -g3 -Wall -O0")
ENDIF(${APPLE})
########################################################################

#######   Check the compiler and set the compile and link flags  #######
set(CMAKE_BUILD_TYPE Debug)
########################################################################

###############          Output directory          #####################
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}")
########################################################################


###############             Root flags            #####################
if (CMAKE_SYSTEM_NAME MATCHES Darwin)
   EXEC_PROGRAM("sw_vers -productVersion | cut -d . -f 1-2" OUTPUT_VARIABLE MAC_OS_VERSION)
   MESSAGE("--- Found a Mac OS X System ${MAC_OS_VERSION}")
   if (CMAKE_COMPILER_IS_GNUCXX)
      MESSAGE("--- Found GNU compiler collection")

      STRING(COMPARE EQUAL "10.5" "${MAC_OS_VERSION}" MAC_OS_10_5)
      IF(MAC_OS_10_5)
        SET(CMAKE_CXX_FLAGS "-m64")
        SET(CMAKE_Fortran_FLAGS "-m64")
      ENDIF(MAC_OS_10_5)

      SET(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS} -flat_namespace -single_module -undefined dynamic_lookup")
      SET(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -flat_namespace -single_module -undefined dynamic_lookup")

      # Select flags.
      set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
      set(CMAKE_CXX_FLAGS_RELEASE        "-O2")
      set(CMAKE_CXX_FLAGS_DEBUG          "-g -O2 -fno-reorder-blocks -fno-schedule-insns -fno-inline")
      set(CMAKE_CXX_FLAGS_DEBUGFULL      "-g3 -fno-inline -Wnon-virtual-dtor -Wno-long-long -ansi -Wundef -Wcast-align -Wchar-subscripts -Wall -W -Wpointer-arith -Wformat-security -fno-exceptions -fno-check-new -fno-common")
      set(CMAKE_CXX_FLAGS_PROFILE        "-g3 -fno-inline -ftest-coverage -fprofile-arcs")
      set(CMAKE_C_FLAGS_RELWITHDEBINFO   "-O2 -g")
      set(CMAKE_C_FLAGS_RELEASE          "-O2")
      set(CMAKE_C_FLAGS_DEBUG            "-g -O2 -fno-reorder-blocks -fno-schedule-insns -fno-inline")
      set(CMAKE_C_FLAGS_DEBUGFULL        "-g3 -fno-inline -Wno-long-long -std=iso9899:1990 -Wundef -Wcast-align -Werror-implicit-function-declaration -Wchar-subscripts -Wall -W -Wpointer-arith -Wwrite-strings -Wformat-security -Wmissing-format-attribute -fno-common")
      set(CMAKE_C_FLAGS_PROFILE          "-g3 -fno-inline -ftest-coverage -fprofile-arcs")

      endif (CMAKE_COMPILER_IS_GNUCXX)

endif (CMAKE_SYSTEM_NAME MATCHES Darwin)

########################################################################

####################### Application sources ############################
file( GLOB  APP_SOURCES RELATIVE "${CMAKE_BINARY_DIR}" "*.cc" )
########################################################################


########################################################################
include_directories (${CMAKE_BINARY_DIR}/include)
include_directories(${CMAKE_BINARY_DIR}/ressources)
include_directories(${CMAKE_BINARY_DIR})
########################################################################

########################################################################
file(GLOB sources ${CMAKE_BINARY_DIR}/src/*.cc)
file(GLOB headers ${CMAKE_BINARY_DIR}/include/*.hh)
file(GLOB ressources ${CMAKE_BINARY_DIR}/ressources/*.hh)
########################################################################

########################################################################
foreach( testsourcefile ${APP_SOURCES} )

    string( REPLACE ".cc" ".bin" testname ${testsourcefile} )

    add_executable( ${testname} ${testsourcefile}  ${sources} ${headers})

    target_link_libraries( ${testname} ${ROOT_LIBRARIES} ${dm_calc_LIBRARIES} ${NEST_LIBRARIES}  )
    
endforeach( testsourcefile ${APP_SOURCES} )
########################################################################
