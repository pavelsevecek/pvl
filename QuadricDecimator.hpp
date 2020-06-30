#pragma once

#include "Matrix.hpp"
#include "Optional.hpp"
#include "TriangleMesh.hpp"
#include <iostream>
#include <map>
#include <queue>

namespace Pvl {

/// \todo make callable with both Context and EdgeHandle? (distinguish using traits)
template <typename Mesh>
class QuadricDecimator {
public:
    using Float = typename Mesh::Float;
    using Point = Vector<Float, 3>;
    using Plane = Vector<Float, 4>;
    using Quadric = Matrix<Float, 4, 4>;

private:
    std::map<VertexHandle, Quadric> quadrics_;

public:
    QuadricDecimator(Mesh& mesh) {
        for (VertexHandle vh : mesh.vertexRange()) {
            quadrics_[vh] = Quadric::null();
        }

        for (FaceHandle fh : mesh.faceRange()) {
            Point n = mesh.areaNormal(fh);
            const Float area = norm(n);
            if (area == 0) {
                continue;
            }
            n /= area;
            Point p0 = mesh.triangle(fh)[0]; /// \todo fix
            Plane plane;
            plane[0] = n[0];
            plane[1] = n[1];
            plane[2] = n[2];
            plane[3] = -dotProd(p0, n);


            Quadric Q = outerProd(plane, plane) * area;

            for (VertexHandle vh : mesh.vertexRing(fh)) {
                quadrics_[vh] += Q;
            }
        }
    }

    float cost(const Mesh& mesh, const Graph::CollapseContext& context) const {
        Quadric Q = quadrics_.at(context.remaining) + quadrics_.at(context.removed);
        Point target = placement(mesh, context);
        Plane v = homogeneous(target);
        Float result = dotProd(prod(Q, v), v);
        PVL_ASSERT(std::isfinite(result));
        return result;
    }

    Point placement(const Mesh& mesh, const Graph::CollapseContext& context) const {
        Quadric Q = quadrics_.at(context.remaining) + quadrics_.at(context.removed);
        Float eps = 1.e-6 * pow<4>(trace(Q)); /// \todo eps dependent on type of Float?

        Q(0, 3) = Q(1, 3) = Q(2, 3) = 0;
        Q(3, 3) = 1;
        // bit of magic to get dimensionless eps
        Float det = std::abs(determinant(Q));
        if (det > eps) {
            Quadric Qinv = invert(Q);
            return Point(Qinv(3, 0), Qinv(3, 1), Qinv(3, 2));
        } else {
            // select from the endpoints
            Plane v1 = homogeneous(mesh.point(context.remaining));
            Plane v2 = homogeneous(mesh.point(context.removed));
            Float cost1 = dotProd(prod(Q, v1), v1);
            Float cost2 = dotProd(prod(Q, v2), v2);
            Plane v = (cost1 < cost2) ? v1 : v2;
            return euclidean(v);
        }
    }

    void postprocess(const Mesh&, const Graph::CollapseContext& context) {
        quadrics_[context.remaining] += quadrics_[context.removed];
    }
};

} // namespace Pvl
