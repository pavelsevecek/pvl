#pragma once

#include "Matrix.hpp"
#include "TriangleMesh.hpp"
#include <iostream>
#include <map>
#include <queue>

namespace Pvl {

/// \todo make callable with both Context and EdgeHandle? (distinguish using traits)
class QuadricDecimator {
    using Mesh = TriangleMesh<Vec3f, int>;

    std::map<VertexHandle, Mat44f> quadrics_;

public:
    QuadricDecimator(Mesh& mesh) {
        for (VertexHandle vh : mesh.vertexRange()) {
            quadrics_[vh] = Mat44f::null();
        }

        for (FaceHandle fh : mesh.faceRange()) {
            Vec3f n = mesh.normal(fh);
            Vec3f p0 = mesh.triangle(fh)[0];
            Vec4f plane;
            plane[0] = n[0];
            plane[1] = n[1];
            plane[2] = n[2];
            plane[3] = -dot(p0, n);

            /// \todo optimize - compute together with normal
            const float area = mesh.area(fh);
            PVL_ASSERT(area > 0.f);

            Mat44f Q = outerProd(plane, plane) * area;

            std::array<int, 3> idxs = mesh.faceIndices(fh);
            for (int i : idxs) {
                quadrics_[VertexHandle(i)] += Q;
            }
        }
    }

    float cost(const Mesh& mesh, const Graph::CollapseContext& context) const {
        /*Vec3f v1 = mesh.point(context.remaining);
        Vec3f v2 = mesh.point(context.removed);*/
        Mat44f Q = quadrics_.at(context.remaining) + quadrics_.at(context.removed);
        Vec3f target = placement(mesh, context);
        Vec4f v = homogeneous(target);
        float result = dot(Q.transform(v), v);
        // PVL_ASSERT(result >= 0.f);
        return result;
    }

    Vec3f placement(const Mesh& mesh, const Graph::CollapseContext& context) const {
        Mat44f Q = quadrics_.at(context.remaining) + quadrics_.at(context.removed);
        float eps = 1.e-6 * pow<4>(trace(Q));

        Q(0, 3) = Q(1, 3) = Q(2, 3) = 0;
        Q(3, 3) = 1;
        // bit of magic to get dimensionless eps
        float det = std::abs(determinant(Q));
        if (det > eps) {
            Mat44f Qinv = invert(Q);
            return Vec3f(Qinv(3, 0), Qinv(3, 1), Qinv(3, 2));
        } else {
            Vec3f v1 = mesh.point(context.remaining);
            Vec3f v2 = mesh.point(context.removed);
            return 0.5 * (v1 + v2);
        }
    }

    void postprocess(const Mesh&, const Graph::CollapseContext& context) {
        quadrics_[context.remaining] += quadrics_[context.removed];
    }
};

} // namespace Pvl
