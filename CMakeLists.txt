cmake_minimum_required(VERSION 3.5)

project(pvl LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(TBB REQUIRED tbb)

add_executable(pvl
    main.cpp
    Utils.hpp
    Assert.hpp
    Math.hpp
    Vector.hpp
    Matrix.hpp
    Svd.hpp
    Box.hpp
    KdTree.h KdTree.inl.h
    PlyWriter.hpp PlyReader.hpp
    Cloud.hpp
    Kernels.hpp
    TriangleMesh.hpp
    MeshUtils.hpp
    Range.hpp
    Graph.hpp Graph.inl.hpp
    UniformGrid.hpp
    OctreeGrid.hpp
    MarchingCubes.hpp MarchingCubes.cpp
    )

target_link_libraries(pvl PRIVATE tbb)
