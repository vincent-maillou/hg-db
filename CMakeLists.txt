cmake_minimum_required(VERSION 3.16)
project(hypergraph_reorder VERSION 1.0.0 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(BUILD_SHARED_LIBS "Build shared library" ON)
option(BUILD_CLI "Build command-line tool" ON)
option(USE_OPENMP "Enable OpenMP parallelization" ON)
option(ENABLE_LTO "Enable link-time optimization" ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra -Wpedantic")

if(ENABLE_LTO AND CMAKE_BUILD_TYPE STREQUAL "Release")
    include(CheckIPOSupported)
    check_ipo_supported(RESULT lto_supported OUTPUT lto_error)
    if(lto_supported)
        message(STATUS "LTO enabled")
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    else()
        message(WARNING "LTO not supported: ${lto_error}")
    endif()
endif()

if(USE_OPENMP)
    find_package(OpenMP)
endif()

# MT-KaHyPar (master branch, no Boost dependency)
# On HPC pass -DKAHYPAR_DOWNLOAD_TBB=ON if TBB is not available system-wide
set(KAHYPAR_PYTHON OFF CACHE BOOL "" FORCE)
set(KAHYPAR_ENABLE_HIGHEST_QUALITY_PRESET ON CACHE BOOL "" FORCE)
set(KAHYPAR_ENABLE_LARGE_K_PARTITIONING_PRESET ON CACHE BOOL "" FORCE)
set(KAHYPAR_ENABLE_GRAPH_PARTITIONING_PRESET ON CACHE BOOL "" FORCE)
add_subdirectory(external/mt-kahypar EXCLUDE_FROM_ALL)

# SuiteSparse — search SUITESPARSE_ROOT, CONDA_PREFIX, then system paths
set(SS_HINTS ${SUITESPARSE_ROOT} $ENV{CONDA_PREFIX} /usr/local /usr)
find_library(AMD_LIBRARY            amd             HINTS ${SS_HINTS} PATH_SUFFIXES lib lib64 REQUIRED)
find_library(CHOLMOD_LIBRARY        cholmod         HINTS ${SS_HINTS} PATH_SUFFIXES lib lib64 REQUIRED)
find_library(CAMD_LIBRARY           camd            HINTS ${SS_HINTS} PATH_SUFFIXES lib lib64 REQUIRED)
find_library(COLAMD_LIBRARY         colamd          HINTS ${SS_HINTS} PATH_SUFFIXES lib lib64 REQUIRED)
find_library(CCOLAMD_LIBRARY        ccolamd         HINTS ${SS_HINTS} PATH_SUFFIXES lib lib64 REQUIRED)
find_library(SUITESPARSECONFIG_LIBRARY suitesparseconfig HINTS ${SS_HINTS} PATH_SUFFIXES lib lib64 REQUIRED)
find_path(SUITESPARSE_INCLUDE_DIR NAMES amd.h cholmod.h
    HINTS ${SS_HINTS} PATH_SUFFIXES include/suitesparse include REQUIRED)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Library
add_library(hypergraph_reorder
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

target_include_directories(hypergraph_reorder
    PUBLIC  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:include>
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/external/mt-kahypar/include
            ${SUITESPARSE_INCLUDE_DIR}
)

target_link_libraries(hypergraph_reorder
    PRIVATE mtkahypar
    PUBLIC  ${AMD_LIBRARY} ${CHOLMOD_LIBRARY} ${CAMD_LIBRARY}
            ${COLAMD_LIBRARY} ${CCOLAMD_LIBRARY} ${SUITESPARSECONFIG_LIBRARY}
            ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}
)

if(OpenMP_CXX_FOUND AND USE_OPENMP)
    target_link_libraries(hypergraph_reorder PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(hypergraph_reorder PUBLIC USE_OPENMP)
endif()

set_target_properties(hypergraph_reorder PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1
    PUBLIC_HEADER include/hypergraph_reorder.h
)

# CLI
if(BUILD_CLI)
    add_executable(hypergraph_reorder_cli cli/hypergraph_reorder_cli.cpp)
    target_link_libraries(hypergraph_reorder_cli PRIVATE hypergraph_reorder)
    set_target_properties(hypergraph_reorder_cli PROPERTIES OUTPUT_NAME hypergraph_reorder)
    install(TARGETS hypergraph_reorder_cli RUNTIME DESTINATION bin)
endif()

install(TARGETS hypergraph_reorder
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    PUBLIC_HEADER DESTINATION include)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}  |  Shared: ${BUILD_SHARED_LIBS}  |  OpenMP: ${OpenMP_CXX_FOUND}  |  LTO: ${ENABLE_LTO}")
