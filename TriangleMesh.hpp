#pragma once

#include "Graph.hpp"
#include "Vector.hpp"
#include <vector>

namespace Pvl {

template <typename Vec, typename Index>
class TriangleMesh : public Graph {
public:
    using Float = typename Vec::Type;
    static constexpr int Dim = Vec::size();

public:
    TriangleMesh() = default;

    std::vector<Vec> points;
};

} // namespace Pvl
