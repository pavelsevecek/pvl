#pragma once

#include "TriangleMesh.hpp"
#include "Utils.hpp"
#include <iostream>

namespace Pvl {

template <typename ConcurrencyTag = SequentialTag, typename Vec>
void laplacianSmoothing(TriangleMesh<Vec>& mesh, bool preserveBoundary = true, float rer = 1.f) {
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
    if (rer > 0.f) {
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
                    biharmonic[v1] = -delta / neighCnt + laplacian[v1];
                }
            });
    }

    ParallelFor<ConcurrencyTag>()(
        std::size_t(0), mesh.numVertices(), [&mesh, &laplacian, &biharmonic, &rer](std::size_t i) {
            mesh.points[i] += 0.5 * (rer * biharmonic[i] + (1.f - rer) * laplacian[i]);
        });
}

} // namespace Pvl
