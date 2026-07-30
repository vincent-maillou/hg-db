cmake_minimum_required(VERSION 3.16)
project(hypergraph_reorder VERSION 1.0.0 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(BUILD_SHARED_LIBS "Build shared library" ON)
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

# ---- MT-KaHyPar discovery ---------------------------------------------------
# MTKAHYPAR_ROOT is provided by the parent Parallax build system.
# It points to the install prefix of a pre-built MT-KaHyPar.
find_path(MTKAHYPAR_INCLUDE_DIR mtkahypar.h
    HINTS ${MTKAHYPAR_ROOT}/include ${MTKAHYPAR_ROOT}/include/mt-kahypar
    REQUIRED)
find_library(MTKAHYPAR_LIBRARY mtkahypar
    HINTS ${MTKAHYPAR_ROOT}/lib ${MTKAHYPAR_ROOT}/lib64
    REQUIRED)

message(STATUS "MT-KaHyPar include: ${MTKAHYPAR_INCLUDE_DIR}")
message(STATUS "MT-KaHyPar library: ${MTKAHYPAR_LIBRARY}")

# ---- Library ----------------------------------------------------------------
add_library(hypergraph_reorder
    src/types.cpp
    src/sparse_matrix.cpp
    src/graph.cpp
    src/clique_cover.cpp
    src/hypergraph.cpp
    src/partitioner.cpp
    src/reorderer.cpp
    src/c_api.cpp
)

target_include_directories(hypergraph_reorder
    PUBLIC  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:include>
    PRIVATE ${MTKAHYPAR_INCLUDE_DIR}
)

target_link_libraries(hypergraph_reorder
    PRIVATE ${MTKAHYPAR_LIBRARY}
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

install(TARGETS hypergraph_reorder
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    PUBLIC_HEADER DESTINATION include)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}  |  Shared: ${BUILD_SHARED_LIBS}  |  OpenMP: ${OpenMP_CXX_FOUND}  |  LTO: ${ENABLE_LTO}")
