#pragma once

#include "PlyWriter.hpp"
#include "TriangleMesh.hpp"
#include <queue>

namespace Pvl {

template <typename Vec>
class MidpointPlacement {
public:
    Vec operator()(const Vec& v1, const Vec& v2) const {
        return 0.5 * (v1 + v2);
    }
};

template <typename Vec>
class EdgeLengthCost {
public:
    typename Vec::Float operator()(const Vec& v1, const Vec& v2) const {
        return norm(v1 - v2);
    }
};

class EdgeCountStop {
    std::size_t count_;

public:
    EdgeCountStop(std::size_t count)
        : count_(count) {}

    bool operator()(std::size_t collapsed) const {
        return collapsed > count_;
    }
};

class CollapseQueue {
    using EdgeCost = std::pair<float, HalfEdgeHandle>;
    std::map<HalfEdgeHandle, EdgeCost> edges_;
    std::set<EdgeCost> queue_;

public:
    void insert(HalfEdgeHandle eh, const float cost) {
        EdgeCost ec{ cost, eh };
        PVL_ASSERT(edges_.find(eh) == edges_.end());
        queue_.insert(ec);
        edges_[eh] = ec;
    }

    void update(HalfEdgeHandle eh, const float cost) {
        /// \todo fix PVL_ASSERT(edges_.find(eh) != edges_.end());
        if (edges_.find(eh) == edges_.end()) {
            insert(eh, cost);
        } else {
            EdgeCost oldEc = edges_[eh];
            queue_.erase(oldEc);

            EdgeCost newEc{ cost, eh };
            queue_.insert(newEc);
        }
    }

    std::pair<HalfEdgeHandle, float> pop() {
        EdgeCost ec = *queue_.begin();
        queue_.erase(ec);
        edges_.erase(ec.second);
        return std::make_pair(ec.second, ec.first);
    }

    void remove(HalfEdgeHandle eh) {
        PVL_ASSERT(edges_.find(eh) != edges_.end());
        EdgeCost ec = edges_[eh];
        queue_.erase(ec);
        edges_.erase(ec.second);
    }

    bool empty() const {
        return queue_.empty();
    }
};

template <typename Vec, typename Index>
void savePatch(const TriangleMesh<Vec, Index>& mesh, HalfEdgeHandle eh) {
    TriangleMesh<Vec, Index> patch;
    std::set<FaceHandle> faces;
    for (FaceHandle fh : mesh.faceRing(mesh.from(eh))) {
        faces.insert(fh);
    }
    for (FaceHandle fh : mesh.faceRing(mesh.to(eh))) {
        faces.insert(fh);
    }

    std::set<VertexHandle> vertices;
    for (FaceHandle fh : faces) {
        VertexHandle v1 = patch.addVertex();
        VertexHandle v2 = patch.addVertex();
        VertexHandle v3 = patch.addVertex();
        patch.addFace(v1, v2, v3);

        for (VertexHandle vh : mesh.vertexRing(fh)) {
            vertices.insert(vh);
        }

        std::array<int, 3> is = mesh.faceIndices(fh);
        patch.points.push_back(mesh.points[is[0]]);
        patch.points.push_back(mesh.points[is[1]]);
        patch.points.push_back(mesh.points[is[2]]);
    }
    std::ofstream ofs("patch-faces.ply");
    PlyWriter writer(ofs);
    writer << patch;

    std::ofstream out("patch-points.ply");
    out << "ply\n";
    out << "format ascii 1.0\n";
    out << "comment Created by PVL library\n";
    out << "element vertex " << vertices.size() << "\n";
    out << "property float x\n";
    out << "property float y\n";
    out << "property float z\n";
    out << "property uchar red\n";
    out << "property uchar green\n";
    out << "property uchar blue\n";
    out << "element face 0\n";
    out << "property list uchar int vertex_index\n";
    out << "end_header\n";
    for (VertexHandle vh : vertices) {
        Vec p = mesh.points[vh];
        Vec c;
        if (vh == mesh.from(eh)) {
            c = Vec(0, 0, 255);
        } else if (vh == mesh.to(eh)) {
            c = Vec(255, 0, 0);
        } else {
            c = Vec(128, 128, 128);
        }
        out << p[0] << " " << p[1] << " " << p[2] << " " << c[0] << " " << c[1] << " " << c[2]
            << "\n";
    }
}

template <typename Vec, typename Index, typename Cost, typename Placement, typename Stop>
void simplify(TriangleMesh<Vec, Index>& mesh,
    const Cost& cost,
    const Placement& placement,
    const Stop& stop) {

    CollapseQueue queue;
    for (EdgeHandle eh : mesh.edgeRange()) {
        HalfEdgeHandle heh = mesh.halfEdge(eh);
        float c = cost(mesh.points[mesh.from(heh)], mesh.points[mesh.to(heh)]);
        queue.insert(heh, c);
    }

    /*for (HalfEdgeHandle heh : mesh.halfEdgeRange()) {
        EdgeHandle eh = mesh.edge(heh);
        HalfEdgeHandle h = mesh.halfEdge(eh);
        float c = cost(mesh.points[mesh.from(h)], mesh.points[mesh.to(h)]);
        queue.update(h, c);
    }*/

    std::size_t cnt = 0;
    HalfEdgeHandle collapsedEdge;
    float c;
    while (!queue.empty()) {
        std::tie(collapsedEdge, c) = queue.pop();
        if (mesh.removed(collapsedEdge)) {
            continue;
        }
        if (mesh.collapseAllowed(collapsedEdge)) {
            Vec p = placement(
                mesh.points[mesh.from(collapsedEdge)], mesh.points[mesh.to(collapsedEdge)]);
            std::vector<HalfEdgeHandle> ring;
            for (HalfEdgeHandle eh : mesh.halfEdgeRing(mesh.from(collapsedEdge))) {
                ring.push_back(eh);
            }
            for (HalfEdgeHandle eh : mesh.halfEdgeRing(mesh.to(collapsedEdge))) {
                ring.push_back(eh);
            }
            std::cout << "# " << cnt << " Collapsing " << mesh.to(collapsedEdge) << " into "
                      << mesh.from(collapsedEdge) << ", cost = " << c << std::endl;
            /*std::cout << "ring:\n";
            for (HalfEdgeHandle eh : ring) {
                std::cout << mesh.from(eh) << "-" << mesh.to(eh) << "\n";
            }
            std::cout << "faces:\n";
            for (HalfEdgeHandle eh : ring) {
                std::array<int, 3> is = mesh.faceIndices(mesh.left(eh));
                std::cout << is[0] << " " << is[1] << " " << is[2] << std::endl;
            }
            std::cout << std::endl;
             savePatch(mesh, collapsedEdge);*/
            mesh.collapse(collapsedEdge, p);
            for (HalfEdgeHandle heh : ring) {
                /// \todo fix this
                if (!mesh.valid(heh)) {
                    continue;
                }
                HalfEdgeHandle heh1 = mesh.halfEdge(mesh.edge(heh));
                if (!mesh.removed(heh1)) {
                    float c = cost(mesh.points[mesh.from(heh1)], mesh.points[mesh.to(heh1)]);
                    queue.update(heh1, c);
                }
            }
            if (stop(cnt++)) {
                return;
            }
        }
    }
}

} // namespace Pvl
