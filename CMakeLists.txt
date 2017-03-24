cmake_minimum_required(VERSION 2.8)
project(gmx_hole C)

# ============== Eclipse CDT Error Handling======
# IF(CMAKE_COMPILER_IS_GNUCC)
#  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmessage-length=0")
# ENDIF(CMAKE_COMPILER_IS_GNUCC)
# IF(CMAKE_COMPILER_IS_GNUCXX)
#  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0")
# ENDIF(CMAKE_COMPILER_IS_GNUCXX)
# set(CMAKE_VERBOSE_MAKEFILE ON)
# ===============================

set(gmx_hole_MAJOR_VERSION 2)
set(gmx_hole_MINOR_VERSION 0)
set(gmx_hole_PATCH_VERSION 0)
set(gmx_hole_VERSION
  ${gmx_hole_MAJOR_VERSION}.${gmx_hole_MINOR_VERSION}.${gmx_hole_PATCH_VERSION})

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")
if(WIN32 AND NOT CYGWIN)
	set(DEF_INSTALL_CMAKE_DIR CMake)
else()
set(DEF_INSTALL_CMAKE_DIR lib/CMake/gmx_hole)
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH "Installation directory for CMake files")

# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE CMAKE)
	set(var INSTALL_${p}_DIR)
	if(NOT IS_ABSOLUTE "${${var}}")
		set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
	endif()
endforeach()

# Cmake modules/macros are in a subdirectory to keep this file cleaner
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# set up include-directories
include_directories(
	"${gmx_hole_SOURCE_DIR}"   # to find src/*.h
	"${gmx_hole_BINARY_DIR}")  # to find src/config.h

# Add sub-directories
add_subdirectory(src)


# ============== Eclipse CDT Error Handling======
#IF(CMAKE_COMPILER_IS_GNUCC)
#  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmessage-length=0")
#ENDIF(CMAKE_COMPILER_IS_GNUCC)
#IF(CMAKE_COMPILER_IS_GNUCXX)
#  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0")
#ENDIF(CMAKE_COMPILER_IS_GNUCXX)
#set(CMAKE_VERBOSE_MAKEFILE ON)

# ===============================
