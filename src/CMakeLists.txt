
# Checking GROMACS PAth
if (NOT DEFINED GMX_PATH AND (NOT DEFINED GMX_LIB OR NOT DEFINED GMX_INCLUDE) )
message(FATAL_ERROR "\nGromacs path is not defined.....\nUSE:\n -DGMX_PATH=/path/to/gromacs \nOR\n -DGMX_LIB=/path/to/gromacs/lib -DGMX_INCLUDE=/path/to/gromacs/include")
endif()


if(DEFINED GMX_PATH)
	set(GMX_INCLUDE "${GMX_PATH}/include")
	set(GMX_LIB "${GMX_PATH}/lib")
endif()

include_directories("${GMX_INCLUDE}/gromacs")



# For library
# add_library(gmxhole SHARED ExtractData.h ExtractData.c gmx_hole.c)
# set_target_properties(gmxhole PROPERTIES PUBLIC_HEADER "ExtractData.h; ${CMAKE_CURRENT_BINARY_DIR}/config.h")

# Searching FFTW library
find_library(LIBFFTW fftw3f HINTS ${FFTW_LIB})
if(DEFINED LIBFFTW-NOTFOUND)
message(FATAL_ERROR "\nFFTW library file libfftw3f.so or libfftw3f.a not found at standard location\nUSE:\n -DFFTW_LIB=/path/to/fftw3/lib")
endif()

# Searching GROMACS Library
find_library(LIBGMX gmx HINTS ${GMX_LIB})
find_library(LIBGMXANA gmxana HINTS ${GMX_LIB})
find_library(LIBMD md HINTS ${GMX_LIB})

# For Library
#install(TARGETS gmxhole
#		EXPORT gmxholeTargets
#		RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
#		LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib
#		PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/gmx_hole"
#		COMPONENT dev)

add_executable(gmx_hole gmx_hole.c ExtractData.h ExtractData.c)


target_link_libraries(gmx_hole m ${LIBFFTW} ${LIBGMX} ${LIBGMXANA} ${LIBMD} pthread dl)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")

install(TARGETS gmx_hole
  # IMPORTANT: Add the bar executable to the "export-set"
        EXPORT gmx_holeTargets
	  	RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin)

# Adding external program hole
# add_executable(hole IMPORTED)
# install(PROGRAMS hole DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin)
