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

    // to -> from
    void collapse(HalfEdgeHandle eh, const Vec& placement) {
        points[to(eh)] = points[from(eh)] = placement;
        Graph::collapse(eh);
    }
};

} // namespace Pvl
