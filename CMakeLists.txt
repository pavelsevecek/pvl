cmake_minimum_required(VERSION 3.5)

project(pvl LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

find_package(TBB REQUIRED tbb)

find_package(Catch2 REQUIRED catch2)

#find_package( Boost 1.40 COMPONENTS bgl REQUIRED )
#include_directories( ${Boost_INCLUDE_DIR} )
include_directories("/usr/include/catch2/")

add_executable(pvl
    main.cpp
    Utils.hpp
    Assert.hpp
    Math.hpp
    Vector.hpp
    Matrix.hpp
    Svd.hpp
    Box.hpp
    KdTree.hpp KdTree.inl.hpp
    PlyWriter.hpp PlyReader.hpp
    Cloud.hpp
    CloudUtils.hpp
    Kernels.hpp
    TriangleMesh.hpp
    Refinement.hpp
    Simplification.hpp
    Range.hpp
    Graph.hpp Graph.inl.hpp
    UniformGrid.hpp
    OctreeGrid.hpp
    MarchingCubes.hpp MarchingCubes.cpp
    )

target_link_libraries(pvl PRIVATE tbb)
