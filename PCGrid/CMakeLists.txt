cmake_minimum_required(VERSION 3.10)
project(PCGrid VERSION 0.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Find required packages
find_package(Eigen3 REQUIRED)
message(STATUS "Found Eigen3: ${EIGEN3_INCLUDE_DIR}")

# 指定 Boost 路径（自定义安装目录）
set(BOOST_ROOT "/usr/local")
set(Boost_INCLUDE_DIR "/usr/local/include")
set(Boost_LIBRARY_DIRS "/usr/local/lib")
set(Boost_NO_SYSTEM_PATHS ON)

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/src/nanoflann
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party/libLAS/include
    ${EIGEN3_INCLUDE_DIR}
    ${Boost_INCLUDE_DIR}
)

# Add libLAS as subdirectory
# Set options for libLAS build
set(WITH_UTILITIES OFF CACHE BOOL "Build libLAS utilities" FORCE)
set(WITH_TESTS OFF CACHE BOOL "Build libLAS tests" FORCE)
set(WITH_LASZIP OFF CACHE BOOL "Build with LASzip support" FORCE)
set(WITH_GDAL OFF CACHE BOOL "Build with GDAL support" FORCE)
set(WITH_GEOTIFF OFF CACHE BOOL "Build with GeoTIFF support" FORCE)

add_subdirectory(third_party/libLAS)

# Collect all source files
set(SOURCES
    src/Las2PowerLine.cpp
    src/LasFileIo.cpp
    src/PCGrid.cpp
    src/PointCloutFeature.cpp
)

# Collect all header files
set(HEADERS
    src/base.h
    src/DenseCluster.h
    src/LasFileIo.h
    src/PointCloutFeature.h
    src/nanoflann/nanoflann.hpp
)

# Create executable
add_executable(Las2PowerLine ${SOURCES} ${HEADERS})

# Link libraries
target_link_libraries(Las2PowerLine PRIVATE las)

# Add compiler flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(Las2PowerLine PRIVATE
        -Wall
        -Wextra
        -Wno-unused-parameter
        -Wno-unused-variable
    )
endif()

# Set output directory
set_target_properties(Las2PowerLine PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)


# Print summary
message(STATUS "")
message(STATUS "========================================")
message(STATUS "PCGrid Configuration Summary")
message(STATUS "========================================")
message(STATUS "  CMake version:          ${CMAKE_VERSION}")
message(STATUS "  System:                 ${CMAKE_SYSTEM_NAME}")
message(STATUS "  C++ Compiler:           ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "  C++ Standard:           ${CMAKE_CXX_STANDARD}")
message(STATUS "  Build type:             ${CMAKE_BUILD_TYPE}")
message(STATUS "  Eigen3 include dir:     ${EIGEN3_INCLUDE_DIR}")
message(STATUS "  Boost include dir:      ${Boost_INCLUDE_DIR}")
message(STATUS "========================================")
message(STATUS "")
