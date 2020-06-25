#pragma once

#include "Graph.hpp"
#include "Vector.hpp"
#include <vector>

namespace Pvl {

template <typename Vec, typename Index>
class TriangleMesh : public Graph {
public:
    using Point = Vec;
    using Float = typename Vec::Float;
    static constexpr int Dim = Vec::size();

public:
    TriangleMesh() = default;

    std::vector<Vec> points;

    Point point(VertexHandle vh) const {
        return points[vh];
    }

    std::array<Point, 3> triangle(FaceHandle fh) const {
        std::array<VertexHandle, 3> idxs = faceVertices(fh);
        Point p0 = point(idxs[0]);
        Point p1 = point(idxs[1]);
        Point p2 = point(idxs[2]);
        return { p0, p1, p2 };
    }

    Point normal(FaceHandle fh) const {
        std::array<Point, 3> tr = triangle(fh);
        return normalize(crossProd(tr[1] - tr[0], tr[2] - tr[0]));
    }

    float area(FaceHandle fh) const {
        std::array<Vec3f, 3> tr = triangle(fh);
        return 0.5f * norm(crossProd(tr[1] - tr[0], tr[2] - tr[0]));
    }

    // to -> from
    void collapse(EdgeHandle eh, const Vec& placement) {
        HalfEdgeHandle heh = halfEdge(eh);
        points[to(heh)] = points[from(heh)] = placement;
        Graph::collapse(eh);
    }
};

} // namespace Pvl
