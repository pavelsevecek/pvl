#pragma once

#include "Vector.hpp"
#include <memory>

namespace Pvl {

template <typename Point>
class OctreeGrid {
    class Node {
        std::unique_ptr<std::array<Node, 8>> children_;
    };

    Node root_;
};

} // namespace Pvl
