#pragma once

#include "TriangleMesh.hpp"
#include "Utils.hpp"
#include <iostream>

namespace Pvl {

template <typename ConcurrencyTag = SequentialTag, typename Vec, typename Index>
void laplacianSmoothing(TriangleMesh<Vec, Index>& mesh, bool preserveBoundary = true) {
    std::vector<Vec> laplacian(mesh.numVertices(), Vec(0));
    ParallelForEach<ConcurrencyTag>()(
        mesh.vertexRange(), [&mesh, &laplacian, preserveBoundary](VertexHandle v1) {
            if (preserveBoundary && mesh.boundary(v1)) {
                return;
            }
            Vec3f delta(0, 0, 0);
            std::size_t neighCnt = 0;
            for (VertexHandle v2 : mesh.vertexRing(v1)) {
                delta += mesh.points[v2];
                ++neighCnt;
            }
            if (neighCnt > 0) {
                laplacian[v1] = delta / neighCnt - mesh.points[v1];
            }
        });
    std::vector<Vec> biharmonic(mesh.numVertices(), Vec(0));
    ParallelForEach<ConcurrencyTag>()(
        mesh.vertexRange(), [&mesh, &laplacian, &biharmonic, preserveBoundary](VertexHandle v1) {
            if (preserveBoundary && mesh.boundary(v1)) {
                return;
            }
            Vec3f delta(0, 0, 0);
            // float weight = 0.;
            std::size_t neighCnt = 0;
            for (VertexHandle v2 : mesh.vertexRing(v1)) {
                delta += laplacian[v2];
                ++neighCnt;
            }
            if (neighCnt > 0) {
                float norm = 1.; // (1. + )
                biharmonic[v1] = -delta / neighCnt + laplacian[v1];
            }
        });

    ParallelFor<ConcurrencyTag>()(Index(0), Index(mesh.numVertices()), [&mesh, &biharmonic](std::size_t i) {
        mesh.points[i] += 0.5f * biharmonic[i];
    });
}

} // namespace Pvl
