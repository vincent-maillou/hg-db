cmake_minimum_required(VERSION 3.16)
project(hypergraph_reorder VERSION 1.0.0 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Options
option(BUILD_SHARED_LIBS "Build shared library" ON)
option(BUILD_CLI "Build command-line tool" ON)
option(BUILD_TESTS "Build unit tests" OFF)
option(BUILD_BENCHMARKS "Build benchmarks" OFF)
option(USE_OPENMP "Enable OpenMP parallelization" ON)
option(ENABLE_LTO "Enable link-time optimization" ON)

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra -Wpedantic")

# Enable LTO for release builds
if(ENABLE_LTO AND CMAKE_BUILD_TYPE STREQUAL "Release")
    include(CheckIPOSupported)
    check_ipo_supported(RESULT lto_supported OUTPUT lto_error)
    if(lto_supported)
        message(STATUS "Link-time optimization enabled")
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    else()
        message(WARNING "LTO not supported: ${lto_error}")
    endif()
endif()

# Find OpenMP
if(USE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        message(STATUS "OpenMP found: ${OpenMP_CXX_VERSION}")
    else()
        message(WARNING "OpenMP not found, parallelization disabled")
    endif()
endif()

# MT-KaHyPar configuration
set(KAHYPAR_DOWNLOAD_TBB ON CACHE BOOL "Download TBB for MT-KaHyPar")
set(KAHYPAR_PYTHON OFF CACHE BOOL "Disable Python bindings")
set(KAHYPAR_ENABLE_HIGHEST_QUALITY_PRESET ON CACHE BOOL "Enable highest quality preset")
set(KAHYPAR_ENABLE_LARGE_K_PARTITIONING_PRESET ON CACHE BOOL "Enable large k preset")
set(KAHYPAR_ENABLE_GRAPH_PARTITIONING_PRESET ON CACHE BOOL "Enable graph partitioning")

# Add MT-KaHyPar as subdirectory
add_subdirectory(external/mt-kahypar EXCLUDE_FROM_ALL)
message(STATUS "MT-KaHyPar included from external/mt-kahypar")

# Find SuiteSparse components
# Allow specifying SUITESPARSE_ROOT or use environment's CONDA_PREFIX
set(SUITESPARSE_HINT_PATHS
    ${SUITESPARSE_ROOT}
    $ENV{CONDA_PREFIX}
    /usr/local
    /usr
)

find_library(AMD_LIBRARY amd HINTS ${SUITESPARSE_HINT_PATHS} PATH_SUFFIXES lib lib64 REQUIRED)
find_library(CHOLMOD_LIBRARY cholmod HINTS ${SUITESPARSE_HINT_PATHS} PATH_SUFFIXES lib lib64 REQUIRED)
find_library(CAMD_LIBRARY camd HINTS ${SUITESPARSE_HINT_PATHS} PATH_SUFFIXES lib lib64 REQUIRED)
find_library(COLAMD_LIBRARY colamd HINTS ${SUITESPARSE_HINT_PATHS} PATH_SUFFIXES lib lib64 REQUIRED)
find_library(SUITESPARSECONFIG_LIBRARY suitesparseconfig HINTS ${SUITESPARSE_HINT_PATHS} PATH_SUFFIXES lib lib64 REQUIRED)
find_library(CCOLAMD_LIBRARY ccolamd HINTS ${SUITESPARSE_HINT_PATHS} PATH_SUFFIXES lib lib64 REQUIRED)

find_path(SUITESPARSE_INCLUDE_DIR
    NAMES amd.h cholmod.h
    HINTS ${SUITESPARSE_HINT_PATHS}
    PATH_SUFFIXES include/suitesparse include
    REQUIRED)

message(STATUS "SuiteSparse include: ${SUITESPARSE_INCLUDE_DIR}")
message(STATUS "AMD library: ${AMD_LIBRARY}")
message(STATUS "CHOLMOD library: ${CHOLMOD_LIBRARY}")

# Source files
set(LIBRARY_SOURCES
    src/types.cpp
    src/sparse_matrix.cpp
    src/graph.cpp
    src/clique_cover.cpp
    src/hypergraph.cpp
    src/partitioner.cpp
    src/ordering.cpp
    src/reorderer.cpp
    src/io/mtx_reader.cpp
    src/io/metis_reader.cpp
    src/io/csr_io.cpp
    src/c_api.cpp
)

# Library target
add_library(hypergraph_reorder ${LIBRARY_SOURCES})

target_include_directories(hypergraph_reorder
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/external/mt-kahypar/include
        ${SUITESPARSE_INCLUDE_DIR}
)

target_link_libraries(hypergraph_reorder
    PRIVATE
        mtkahypar
    PUBLIC
        ${AMD_LIBRARY}
        ${CHOLMOD_LIBRARY}
        ${CAMD_LIBRARY}
        ${COLAMD_LIBRARY}
        ${CCOLAMD_LIBRARY}
        ${SUITESPARSECONFIG_LIBRARY}
)

# Link BLAS and LAPACK (required by CHOLMOD)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
target_link_libraries(hypergraph_reorder PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

# Add OpenMP if available
if(OpenMP_CXX_FOUND AND USE_OPENMP)
    target_link_libraries(hypergraph_reorder PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(hypergraph_reorder PUBLIC USE_OPENMP)
    message(STATUS "Building with OpenMP support")
endif()

# Set library properties
set_target_properties(hypergraph_reorder PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1
    PUBLIC_HEADER "include/hypergraph_reorder.h"
)

# CLI tool
if(BUILD_CLI)
    add_executable(hypergraph_reorder_cli cli/hypergraph_reorder_cli.cpp)
    target_link_libraries(hypergraph_reorder_cli PRIVATE hypergraph_reorder)
    set_target_properties(hypergraph_reorder_cli PROPERTIES OUTPUT_NAME hypergraph_reorder)

    # Install CLI
    install(TARGETS hypergraph_reorder_cli
            RUNTIME DESTINATION bin)
endif()

# Tests
if(BUILD_TESTS)
    enable_testing()
    # Add test subdirectory if needed
    # add_subdirectory(tests)
endif()

# Benchmarks
if(BUILD_BENCHMARKS)
    # add_subdirectory(benchmarks)
endif()

# Installation (basic - no CMake export due to MT-KaHyPar bundling)
install(TARGETS hypergraph_reorder
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        PUBLIC_HEADER DESTINATION include)

install(DIRECTORY include/hypergraph_reorder
        DESTINATION include
        FILES_MATCHING PATTERN "*.hpp")

# Print configuration summary
message(STATUS "")
message(STATUS "Configuration Summary:")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  Shared libs: ${BUILD_SHARED_LIBS}")
message(STATUS "  CLI tool: ${BUILD_CLI}")
message(STATUS "  OpenMP: ${USE_OPENMP}")
message(STATUS "  LTO: ${ENABLE_LTO}")
message(STATUS "")
