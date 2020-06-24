#pragma once

#include "Graph.hpp"
#include "Vector.hpp"
#include <vector>

namespace Pvl {

template <typename Vec, typename Index>
class TriangleMesh : public Graph {
public:
    using Float = typename Vec::Float;
    static constexpr int Dim = Vec::size();

public:
    TriangleMesh() = default;

    std::vector<Vec> points;

    Vec3f point(VertexHandle vh) const {
        return points[vh];
    }

    std::array<Vec3f, 3> triangle(FaceHandle fh) const {
        std::array<int, 3> idxs = faceIndices(fh);
        Vec3f p0 = points[idxs[0]];
        Vec3f p1 = points[idxs[1]];
        Vec3f p2 = points[idxs[2]];
        return { p0, p1, p2 };
    }

    Vec3f normal(FaceHandle fh) const {
        std::array<Vec3f, 3> tr = triangle(fh);
        return normalize(cross(tr[1] - tr[0], tr[2] - tr[0]));
    }

    float area(FaceHandle fh) const {
        std::array<Vec3f, 3> tr = triangle(fh);
        return 0.5f * norm(cross(tr[1] - tr[0], tr[2] - tr[0]));
    }

    // to -> from
    void collapse(HalfEdgeHandle eh, const Vec& placement) {
        points[to(eh)] = points[from(eh)] = placement;
        Graph::collapse(eh);
    }
};

} // namespace Pvl
