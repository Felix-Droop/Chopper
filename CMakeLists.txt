cmake_minimum_required (VERSION 3.8...3.18)

## CUSTOMISE

# Define the application name and version.
project (chopper VERSION 1.0.0)

# make C++20 the default standard
if (NOT CMAKE_CXX_FLAGS)
    set (CMAKE_CXX_FLAGS "-std=c++20" CACHE STRING "default C++ standard" FORCE)
endif ()

# add recommended C++ compiler flags for SIMDe
set (CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -fopenmp-simd -DSIMDE_ENABLE_OPENMP -march=native" CACHE STRING 
"recommended compiler flags for SIMDe" FORCE)

# add C++ compiler flags for OpenMP
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp" CACHE STRING "OpenMP to compiler flag" FORCE)

## BUILD

# Make Release default build type
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release CACHE STRING
         "Choose the type of build, options are: Debug Release RelWithDebInfo"
         FORCE)
endif ()

# Specify the directories where to store the built archives, libraries and executables
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Messages
string (ASCII 27 Esc)
set (FontBold "${Esc}[1m")
set (FontReset "${Esc}[m")

# Dependency: SeqAn3.
find_package (SeqAn3 QUIET REQUIRED HINTS lib/seqan3/build_system)

# Dependency: SeqAn
# Note: includes seqan-config.cmake manually to have
# SEQAN_DEFINITIONS, SEQAN_CXX_FLAGS_LIST, SEQAN_LIBRARIES, SEQAN_INCLUDE_DIRS
# defined. We don't want any seqan version but the submodule one!
set(CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/lib/seqan/util/cmake")
set(CMAKE_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/lib/seqan/include/")
include (${CMAKE_SOURCE_DIR}/lib/seqan/util/cmake/seqan-config.cmake)

# Dependency: xxHash
# We don't need the command line application
set(XXHASH_BUILD_XXHSUM OFF)
add_subdirectory("${CMAKE_SOURCE_DIR}/lib/xxHash/cmake_unofficial" "${CMAKE_BINARY_DIR}/lib/xxHash/" EXCLUDE_FROM_ALL)

add_subdirectory (src)
message (STATUS "${FontBold}You can run `make` to build the application.${FontReset}")

## DOCUMENTATION

add_subdirectory (doc EXCLUDE_FROM_ALL)

## TEST

enable_testing ()
add_subdirectory (test EXCLUDE_FROM_ALL)
