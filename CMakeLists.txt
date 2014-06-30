cmake_minimum_required(VERSION 2.8)
project(gmx_hole C)

set(gmx_hole_MAJOR_VERSION 1)
set(gmx_hole_MINOR_VERSION 0)
set(gmx_hole_PATCH_VERSION 0)
set(gmx_hole_VERSION
  ${FOOBAR_MAJOR_VERSION}.${FOOBAR_MINOR_VERSION}.${FOOBAR_PATCH_VERSION})

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

# set up include-directories
include_directories(
	"${gmx_hole_SOURCE_DIR}"   # to find src/*.h
	"${gmx_hole_BINARY_DIR}")  # to find src/config.h
	 
# Add sub-directories
add_subdirectory(src)


# ============== Eclipse CDT Error Handling======
IF(CMAKE_COMPILER_IS_GNUCC)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmessage-length=0")
ENDIF(CMAKE_COMPILER_IS_GNUCC)
IF(CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)
set(CMAKE_VERBOSE_MAKEFILE ON)

# ===============================
 
# Add all targets to the build-tree export set
export(TARGETS gmx_hole FILE "${PROJECT_BINARY_DIR}/gmx_holeTargets.cmake")
    
# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE gmx_hole)


# Create the gmx_holeConfig.cmake and gmx_holeConfigVersion files
file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}" "${INSTALL_INCLUDE_DIR}")

# ... for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
configure_file(gmx_holeConfig.cmake.in "${PROJECT_BINARY_DIR}/gmx_holeConfig.cmake" @ONLY)
# ... for the install tree
set(CONF_INCLUDE_DIRS "\${gmx_hole_CMAKE_DIR}/${REL_INCLUDE_DIR}")
configure_file(gmx_holeConfig.cmake.in "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/gmx_holeConfig.cmake" @ONLY)
# ... for both
configure_file(gmx_holeConfigVersion.cmake.in "${PROJECT_BINARY_DIR}/gmx_holeConfigVersion.cmake" @ONLY)

# Install the gmx_holeConfig.cmake and gmx_holeConfigVersion.cmake
install(FILES "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/gmx_holeConfig.cmake" 
              "${PROJECT_BINARY_DIR}/gmx_holeConfigVersion.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)
	   
# Install the export set for use with the install-tree
install(EXPORT gmx_holeTargets DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)
