#pragma once

#include "Assert.hpp"
#include "Range.hpp"
#include <algorithm>
#include <array>
#include <deque>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <vector>


namespace Pvl {

template <typename TObject, typename Index = int>
class Handle {
    Index idx_;

public:
    using Object = TObject;

    Handle() = default;

    explicit Handle(const Index idx)
        : idx_(idx) {}

    Handle& operator++() {
        ++idx_;
        return *this;
    }

    Index index() const {
        return idx_;
    }

    operator Index() const {
        return idx_;
    }
};

class HalfEdge;
class Vertex;
class Face;
class Edge;
class Graph;

using HalfEdgeHandle = Handle<HalfEdge>;
using VertexHandle = Handle<Vertex>;
using FaceHandle = Handle<Face>;
using EdgeHandle = Handle<Edge>;

struct HalfEdge {
    using Handle = HalfEdgeHandle;

    VertexHandle to;
    FaceHandle left;
    HalfEdgeHandle opposite;
    HalfEdgeHandle next;
    HalfEdgeHandle prev;

    HalfEdge(VertexHandle to)
        : to(to) {
        left = FaceHandle(-1);
    }

    bool boundary() const {
        return opposite < 0;
    }

    /* bool removed() const {
          return left < 0;
      }*/
};

struct Vertex {
    using Handle = VertexHandle;

    // boundary emanating edge for boundary vertices, any emanating vertex for internal
    // vertices
    HalfEdgeHandle edge;

    Vertex() {
        edge = HalfEdgeHandle(-1);
    }

    Vertex(HalfEdgeHandle eh)
        : edge(eh) {}

    /*bool valid() const {
        return edge >= 0;
    }*/

    struct EndTag {};

    class IteratorBase {
    protected:
        const Graph& graph_;
        HalfEdgeHandle eh_;

        enum class Status {
            BEGIN,
            ITERS,
            END,
        } status_;

    public:
        IteratorBase(const Graph& graph, VertexHandle vh);

        IteratorBase(const Graph& graph, VertexHandle vh, EndTag);

        IteratorBase& operator++();

        bool operator!=(const IteratorBase& other) const;
    };

    class HalfEdgeIterator : public IteratorBase {
    public:
        using IteratorBase::IteratorBase;

        using iterator_category = std::forward_iterator_tag;
        using value_type = HalfEdgeHandle;
        using difference_type = int;
        using reference = HalfEdgeHandle&;
        using pointer = HalfEdgeHandle*;

        HalfEdgeHandle operator*() const;
    };

    class VertexIterator : public IteratorBase {
    public:
        using IteratorBase::IteratorBase;

        using iterator_category = std::forward_iterator_tag;
        using value_type = VertexHandle;
        using difference_type = int;
        using reference = VertexHandle&;
        using pointer = VertexHandle*;

        VertexHandle operator*() const;
    };

    class FaceIterator : public IteratorBase {
    public:
        using IteratorBase::IteratorBase;

        using iterator_category = std::forward_iterator_tag;
        using value_type = FaceHandle;
        using difference_type = int;
        using reference = FaceHandle&;
        using pointer = FaceHandle*;

        FaceHandle operator*() const;
    };

    class EdgeIterator : public IteratorBase {
    public:
        using IteratorBase::IteratorBase;

        using iterator_category = std::forward_iterator_tag;
        using value_type = EdgeHandle;
        using difference_type = int;
        using reference = EdgeHandle&;
        using pointer = EdgeHandle*;

        EdgeIterator& operator++();

        EdgeHandle operator*() const;
    };


    using HalfEdgeRange = Range<HalfEdgeIterator>;
    using EdgeRange = Range<EdgeIterator>;
    using VertexRange = Range<VertexIterator>;
    using FaceRange = Range<FaceIterator>;
};

static_assert(sizeof(Vertex) == sizeof(Handle<Vertex>), "error");

struct Face {
    HalfEdgeHandle edge;

    Face()
        : edge(-1) {}

    Face(HalfEdgeHandle eh)
        : edge(eh) {}

    /* bool valid() const {
         return edge >= 0;
     }*/

    struct EndTag {};

    class IteratorBase {
    protected:
        const Graph& graph_;
        HalfEdgeHandle eh_;
        int steps_;

    public:
        IteratorBase(const Graph& graph, FaceHandle fh);

        IteratorBase(const Graph& graph, FaceHandle fh, EndTag);

        IteratorBase& operator++();

        bool operator!=(const IteratorBase& other) const;
    };

    class HalfEdgeIterator : public IteratorBase {
    public:
        using IteratorBase::IteratorBase;

        using iterator_category = std::forward_iterator_tag;
        using value_type = HalfEdgeHandle;
        using difference_type = int;
        using reference = HalfEdgeHandle&;
        using pointer = HalfEdgeHandle*;

        HalfEdgeHandle operator*() const;
    };

    class VertexIterator : public IteratorBase {
    public:
        using IteratorBase::IteratorBase;

        using iterator_category = std::forward_iterator_tag;
        using value_type = VertexHandle;
        using difference_type = int;
        using reference = VertexHandle&;
        using pointer = VertexHandle*;

        VertexHandle operator*() const;
    };

    class FaceIterator : public IteratorBase {
    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = FaceHandle;
        using difference_type = int;
        using reference = FaceHandle&;
        using pointer = FaceHandle*;

        FaceIterator(const Graph& graph, FaceHandle fh);

        FaceIterator(const Graph& graph, FaceHandle fh, EndTag);

        FaceHandle operator*() const;

        FaceIterator& operator++();
    };


    using HalfEdgeRange = Range<HalfEdgeIterator>;
    using VertexRange = Range<VertexIterator>;
    using FaceRange = Range<FaceIterator>;
};

struct Edge {};

class Graph {
    friend class Vertex;
    friend class Face;

    std::vector<Vertex> vertices_;
    std::vector<Face> faces_;
    std::vector<HalfEdge> halfEdges_;

    // maps vertex to the set of emanating boundary edges
    std::map<VertexHandle, std::set<HalfEdgeHandle>> boundaryEdges_;

public:
    HalfEdgeHandle next(HalfEdgeHandle eh) const {
        PVL_ASSERT(valid(eh));
        return halfEdges_[eh].next;
    }
    HalfEdgeHandle prev(HalfEdgeHandle eh) const {
        PVL_ASSERT(valid(eh));
        return halfEdges_[eh].prev;
    }
    HalfEdgeHandle opposite(HalfEdgeHandle eh) const {
        PVL_ASSERT(valid(eh));
        PVL_ASSERT(!boundary(eh));
        return halfEdges_[eh].opposite;
    }
    bool boundary(HalfEdgeHandle eh) const {
        PVL_ASSERT(valid(eh));
        return halfEdges_[eh].boundary();
    }
    bool boundary(EdgeHandle eh) const {
        PVL_ASSERT(valid(eh));
        return boundary(halfEdge(eh));
    }
    bool boundary(VertexHandle vh) const {
        PVL_ASSERT(valid(vh));
        return boundary(vertices_[vh].edge);
    }
    VertexHandle to(HalfEdgeHandle eh) const {
        if (!valid(eh)) {
            std::cout << "Invalid he handle " << halfEdges_[halfEdges_[eh].prev].to << "-"
                      << halfEdges_[eh].to << std::endl;
        }
        PVL_ASSERT(valid(eh));
        return halfEdges_[eh].to;
    }
    VertexHandle from(HalfEdgeHandle eh) const {
        PVL_ASSERT(valid(eh));
        return halfEdges_[prev(eh)].to;
    }
    FaceHandle left(HalfEdgeHandle eh) const {
        PVL_ASSERT(valid(eh));
        return halfEdges_[eh].left;
    }
    FaceHandle right(HalfEdgeHandle eh) const {
        PVL_ASSERT(!boundary(eh));
        return halfEdges_[opposite(eh)].left;
    }
    HalfEdgeHandle from(VertexHandle vh) const {
        PVL_ASSERT(valid(vh));
        return vertices_[vh].edge;
    }
    HalfEdgeHandle to(VertexHandle vh) const {
        PVL_ASSERT(valid(vh));
        return prev(from(vh));
    }
    // returns halfedge from vh1 to vh2, or -1 if not exists
    HalfEdgeHandle halfEdge(VertexHandle vh1, VertexHandle vh2) const {
        PVL_ASSERT(valid(vh1) && valid(vh2));
        for (HalfEdgeHandle eh : halfEdgeRing(vh1)) {
            if (to(eh) == vh2) {
                return eh;
            }
        }
        return HalfEdgeHandle(-1);
    }
    // returns any half edge in a face
    HalfEdgeHandle halfEdge(FaceHandle fh) const {
        PVL_ASSERT(valid(fh));
        return faces_[fh].edge;
    }
    // returns any half edge for given edge
    HalfEdgeHandle halfEdge(EdgeHandle eh) const {
        PVL_ASSERT(valid(eh));
        return HalfEdgeHandle(eh.index());
    }
    EdgeHandle edge(HalfEdgeHandle heh) const {
        PVL_ASSERT(valid(heh));
        EdgeHandle eh;
        if (boundary(heh) || from(heh) < to(heh)) {
            eh = EdgeHandle(heh.index());
        } else {
            eh = EdgeHandle(opposite(heh).index());
        }
        PVL_ASSERT(valid(eh));
        return eh;
    }
    EdgeHandle edge(VertexHandle vh1, VertexHandle vh2) const {
        PVL_ASSERT(valid(vh1) && valid(vh2));
        if (vh1 > vh2) {
            std::swap(vh1, vh2);
        }
        HalfEdgeHandle eh = halfEdge(vh1, vh2);
        if (!valid(eh)) {
            // might still be a boundary edge
            eh = halfEdge(vh2, vh1);
            PVL_ASSERT(!valid(eh) || boundary(eh));
        }
        return edge(eh);
    }
    bool valid(FaceHandle fh) const {
        PVL_ASSERT(fh < int(faces_.size()));
        return fh != FaceHandle(-1) && faces_[fh].edge != HalfEdgeHandle(-1);
    }
    bool valid(HalfEdgeHandle heh) const {
        PVL_ASSERT(heh < int(halfEdges_.size()));
        return heh != HalfEdgeHandle(-1) && halfEdges_[heh].left != FaceHandle(-1);
    }
    bool valid(VertexHandle vh) const {
        PVL_ASSERT(vh < int(vertices_.size()));
        return vh != VertexHandle(-1) && vertices_[vh].edge != HalfEdgeHandle(-1);
    }
    bool valid(EdgeHandle eh) const {
        PVL_ASSERT(eh < int(halfEdges_.size()));
        HalfEdgeHandle heh(eh.index());
        return valid(heh) && (boundary(heh) || from(heh) < to(heh));
    }

    template <typename THandle>
    class HandleIterator {
        THandle handle_;

    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = THandle;
        using difference_type = int;
        using reference = THandle&;
        using pointer = THandle*;

        HandleIterator(THandle handle)
            : handle_(handle) {}

        THandle operator*() const {
            return handle_;
        }

        HandleIterator& operator++() {
            ++handle_;
            return *this;
        }

        bool operator==(const HandleIterator& other) const {
            return handle_ == other.handle_;
        }

        bool operator!=(const HandleIterator& other) const {
            return handle_ != other.handle_;
        }
    };


    class EdgeIterator {
        const Graph& graph_;
        HalfEdgeHandle eh_;

    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = HalfEdgeHandle;
        using difference_type = int;
        using reference = HalfEdgeHandle&;
        using pointer = HalfEdgeHandle*;

        EdgeIterator(const Graph& graph, HalfEdgeHandle eh)
            : graph_(graph)
            , eh_(eh) {
            while (!end() && !dereferencable()) {
                ++eh_;
            }
        }

        EdgeHandle operator*() const {
            return EdgeHandle(eh_);
        }

        EdgeIterator& operator++() {
            do {
                ++eh_;
            } while (!end() && !dereferencable());
            return *this;
        }

        bool operator==(const EdgeIterator& other) const {
            return eh_ == other.eh_;
        }

        bool operator!=(const EdgeIterator& other) const {
            return eh_ != other.eh_;
        }

    private:
        bool end() const {
            return eh_ >= int(graph_.halfEdges_.size());
        }

        bool dereferencable() const {
            if (!graph_.valid(eh_)) {
                // other iterators visit removed simplices, so let's be consistent
                return true;
            }
            return graph_.boundary(eh_) || graph_.from(eh_) < graph_.to(eh_);
        }
    };

    using HalfEdgeIterator = HandleIterator<HalfEdgeHandle>;
    using VertexIterator = HandleIterator<VertexHandle>;
    using FaceIterator = HandleIterator<FaceHandle>;

    using HalfEdgeRange = Range<HalfEdgeIterator>;
    using VertexRange = Range<VertexIterator>;
    using FaceRange = Range<FaceIterator>;
    using EdgeRange = Range<EdgeIterator>;

    /// \todo rename to Iterators and Circulators
    HalfEdgeRange halfEdgeRange() const {
        return { HalfEdgeHandle(0), HalfEdgeHandle(halfEdges_.size()) };
    }

    EdgeRange edgeRange() const {
        return { EdgeIterator(*this, HalfEdgeHandle(0)),
            EdgeIterator(*this, HalfEdgeHandle(halfEdges_.size())) };
    }

    VertexRange vertexRange() const {
        return { VertexHandle(0), VertexHandle(vertices_.size()) };
    }

    FaceRange faceRange() const {
        return { FaceHandle(0), FaceHandle(faces_.size()) };
    }

    Vertex::HalfEdgeRange halfEdgeRing(VertexHandle vh) const {
        return { Vertex::HalfEdgeIterator(*this, vh),
            Vertex::HalfEdgeIterator(*this, vh, Vertex::EndTag{}) };
    }

    Vertex::EdgeRange edgeRing(VertexHandle vh) const {
        return { Vertex::EdgeIterator(*this, vh),
            Vertex::EdgeIterator(*this, vh, Vertex::EndTag{}) };
    }

    Vertex::VertexRange vertexRing(VertexHandle vh) const {
        return { Vertex::VertexIterator(*this, vh),
            Vertex::VertexIterator(*this, vh, Vertex::EndTag{}) };
    }

    Vertex::FaceRange faceRing(VertexHandle vh) const {
        return { Vertex::FaceIterator(*this, vh),
            Vertex::FaceIterator(*this, vh, Vertex::EndTag{}) };
    }

    Face::HalfEdgeRange halfEdgeRing(FaceHandle fh) const {
        return { Face::HalfEdgeIterator(*this, fh),
            Face::HalfEdgeIterator(*this, fh, Face::EndTag{}) };
    }

    Face::VertexRange vertexRing(FaceHandle fh) const {
        return { Face::VertexIterator(*this, fh),
            Face::VertexIterator(*this, fh, Face::EndTag{}) };
    }

    Face::FaceRange faceRing(FaceHandle fh) const {
        return { Face::FaceIterator(*this, fh),
            Face::FaceIterator(*this, fh, Face::EndTag{}) };
    }

    std::array<int, 3> faceIndices(FaceHandle fh) const {
        if (!valid(fh)) {
            return { 0, 0, 0 };
        }
        std::array<int, 3> indices;
        int i = 0;
        for (VertexHandle vh : vertexRing(fh)) {
            indices[i] = vh.index();
            ++i;
        }
        return indices;
    }

    VertexHandle addVertex() {
        vertices_.emplace_back(HalfEdgeHandle(-1));
        return VertexHandle(vertices_.size() - 1);
    }

    FaceHandle addFace(VertexHandle vh1, VertexHandle vh2, VertexHandle vh3) {
        // std::cout << "Adding face " << vh1 << ", " << vh2 << ", " << vh3 << std::endl;
        faces_.emplace_back();
        FaceHandle fh = FaceHandle(faces_.size() - 1);
        Face& f = faces_[fh];

        halfEdges_.emplace_back(vh2);
        halfEdges_.emplace_back(vh3);
        halfEdges_.emplace_back(vh1);

        std::array<HalfEdgeHandle, 3> eh = {
            HalfEdgeHandle(halfEdges_.size() - 3),
            HalfEdgeHandle(halfEdges_.size() - 2),
            HalfEdgeHandle(halfEdges_.size() - 1),
        };
        std::array<VertexHandle, 3> vh = { vh1, vh2, vh3 };

        // set face edge
        f.edge = eh[0];

        for (int i = 0; i < 3; ++i) {
            // set emanating edge
            if (!valid(vh[i])) {
                vertices_[vh[i]].edge = eh[i];
            }

            // set next and previous edge
            HalfEdge& e = halfEdges_[eh[i]];
            e.left = fh;
            e.next = eh[(i + 1) % 3];
            e.prev = eh[(i + 2) % 3];
            e.opposite = HalfEdgeHandle(-1);
        }

        for (int i = 0; i < 3; ++i) {
            for (HalfEdgeHandle neh : halfEdgeRing(to(eh[i]))) {
                // std::cout << "Visiting " << neh << "\n";
                if (to(neh) == from(eh[i])) {
                    if (!halfEdges_[neh].boundary()) {
                        std::cout << "Complex edge!\n";
                        break;
                    }
                    halfEdges_[neh].opposite = eh[i];
                    halfEdges_[eh[i]].opposite = neh;
                    break;
                }
            }

            if (halfEdges_[eh[i]].boundary()) {
                // opposite not found by circulation
                std::set<HalfEdgeHandle> removed;
                for (HalfEdgeHandle neh : boundaryEdges_[to(eh[i])]) {
                    if (to(neh) == from(eh[i])) {
                        if (!halfEdges_[neh].boundary()) {
                            std::cout << "Complex edge!\n";
                        }
                        halfEdges_[neh].opposite = eh[i];
                        halfEdges_[eh[i]].opposite = neh;
                        removed.insert(neh);
                        break;
                    }
                }
                for (HalfEdgeHandle neh : removed) {
                    boundaryEdges_[to(eh[i])].erase(neh);
                }
                if (boundaryEdges_[to(eh[i])].empty()) {
                    boundaryEdges_.erase(to(eh[i]));
                }
            }
        }
        for (int i = 0; i < 3; ++i) {
            if (halfEdges_[eh[i]].boundary()) {
                boundaryEdges_[from(eh[i])].insert(eh[i]);
            }

            {
                VertexHandle vh = from(eh[i]);
                if (!boundary(eh[i])) {
                    // find boundary
                    HalfEdgeHandle e = eh[i];
                    do {
                        e = next(opposite(e));
                    } while (e != eh[i] && !boundary(e));
                    vertices_[vh].edge = e;
                }
            }
            {
                HalfEdgeHandle e = from(vh[i]);
                HalfEdgeHandle e0 = from(vh[i]);
                if (!boundary(e)) {
                    // find boundary
                    do {
                        e = next(opposite(e));
                    } while (e != e0 && !boundary(e));
                    vertices_[vh[i]].edge = e;
                }
            }
        }

        return fh;
    }

    struct CollapseContext {
        HalfEdgeHandle edge; // halfedge from "remaining" to "removed"
        VertexHandle removed;
        VertexHandle remaining;
        FaceHandle left;  // left of the halfedge
        FaceHandle right; // right of the halfedge

        CollapseContext(const Graph& graph, const HalfEdgeHandle heh)
            : edge(heh) {
            removed = graph.to(edge);
            remaining = graph.from(edge);
            left = graph.left(edge);
            right = graph.right(edge);
        }
    };


    // collapse to to from
    void collapse(HalfEdgeHandle heh) {
        PVL_ASSERT(collapseAllowed(heh));
        // check();
        CollapseContext context(*this, heh);
        HalfEdgeHandle oheh = opposite(heh);

        /*if (vertices_[context.remaining].edge == heh) {
            HalfEdgeHandle l = opposite(prev(heh));
            std::cout << "Moving emanating edge for vertex " << context.remaining << " to "
                      << from(l) << "-" << to(l) << std::endl;
            vertices_[context.remaining].edge = opposite(prev(heh));
        }*/

        VertexHandle v0 = context.remaining;
        VertexHandle vl = to(next(heh));
        VertexHandle vr = to(next(oheh));

        /// \todo move only if necessary?
        /// \todo fix boundary case
        vertices_[v0].edge = opposite(prev(heh));
        vertices_[vl].edge = next(opposite(prev(heh)));
        vertices_[vr].edge = opposite(next(oheh));

        /*for (HalfEdgeHandle neh : halfEdgeRing(to(heh))) {
            if (to(neh) == from(heh)) {
                continue;
            }
            halfEdges_[prev(neh)].to = from(heh);
        }*/
        HalfEdgeHandle neh = heh;
        do {
            neh = opposite(next(neh));

            /*std::cout << "Reassigning 'to' for edge " << from(neh) << "-" << to(neh)
                      << " from " << halfEdges_[neh].to << " to " << from(heh) << std::endl;*/
            halfEdges_[neh].to = from(heh);
        } while (neh != heh);
        /*edges_[prev(eh)].next = next(opposite(next(eh)));
        edges_[prev(opposite(next(eh)))].next = from(eh);
        edges_[next(opposite(prev(oeh)))].next = next(oeh);*/
        connect(opposite(next(heh)), opposite(prev(heh)));
        connect(opposite(next(oheh)), opposite(prev(oheh)));

        vertices_[context.removed].edge = HalfEdgeHandle(-1);
        remove(context.left);
        remove(context.right);

        // check();
        /*remove(opposite(next(heh)));
        remove(prev(opposite(heh)));*/
        //   remove(heh);
        // remove(to(eh));
    }

    bool collapseAllowed(HalfEdgeHandle eh) {
        PVL_ASSERT(valid(eh));
        if (boundary(eh)) { /// \todo enable boundary collapse
            return false;
        }

        if (boundary(to(eh)) != boundary(from(eh))) {
            // disallow contraction of boundary vertex and inner vertex
            return false;
        }
        std::set<VertexHandle> ring1;
        for (VertexHandle vh : vertexRing(to(eh))) {
            // std::cout << "circ " << vh << std::endl;
            if (vh != from(eh)) {
                ring1.insert(vh);
            }
            for (HalfEdgeHandle neh : halfEdgeRing(vh)) {
                if (boundary(neh)) {
                    return false;
                }
            }
        }
        std::set<VertexHandle> ring2;
        for (VertexHandle vh : vertexRing(from(eh))) {
            if (vh != to(eh)) {
                ring2.insert(vh);
            }
            for (HalfEdgeHandle neh : halfEdgeRing(vh)) {
                if (boundary(neh)) {
                    return false;
                }
            }
        }
        std::vector<VertexHandle> is;
        std::set_intersection(
            ring1.begin(), ring1.end(), ring2.begin(), ring2.end(), std::back_inserter(is));
        if (is.size() != 2) {
            return false;
        }
        VertexHandle vl = to(next(eh));
        VertexHandle vr = to(next(opposite(eh)));
        return ((is[0] == vl && is[1] == vr) || (is[1] == vl && is[0] == vr));
    }

    void remove(FaceHandle fh) {
        PVL_ASSERT(valid(fh));
        HalfEdgeHandle eh = halfEdge(fh);
        HalfEdgeHandle eh0 = eh;
        do {
            // careful not to call next with invalid handle
            HalfEdgeHandle eh1 = eh;
            eh = next(eh);
            halfEdges_[eh1].left = FaceHandle(-1);
        } while (eh != eh0);
        faces_[fh].edge = HalfEdgeHandle(-1);

        /// \todo can create non-manifold graph, re-add boundary edges
    }

    /*void remove(HalfEdgeHandle heh) {
        PVL_ASSERT(valid(heh));
        if (!boundary(heh)) {
            halfEdges_[opposite(heh)].left = FaceHandle(-1);
        }
        halfEdges_[heh].left = FaceHandle(-1);
    }*/

    bool removed(HalfEdgeHandle heh) {
        return halfEdges_[heh].left == -1;
    }

    bool removed(EdgeHandle heh) {
        return halfEdges_[heh].left == -1;
    }

    void collectGarbage() {}

    std::size_t numVertices() const {
        return vertices_.size();
    }

    std::size_t numFaces() const {
        return faces_.size();
    }

    // warning: O(n) complexity
    std::size_t numEdges() const {
        EdgeRange edges = edgeRange();
        return std::distance(edges.begin(), edges.end());
    }

private:
    void connect(HalfEdgeHandle eh1, HalfEdgeHandle eh2) {
        halfEdges_[eh1].opposite = eh2;
        halfEdges_[eh2].opposite = eh1;
    }

    void check() const {
        // consistency checks
        for (VertexHandle vh : vertexRange()) {
            if (!valid(vh)) {
                continue; // removed or isolated
            }

            for (VertexHandle nvh : vertexRing(vh)) {
                if (!valid(nvh)) {
                    std::cout << "vertex " << vh << " connected to invalid vertex " << nvh
                              << std::endl;
                }
            }
            for (HalfEdgeHandle heh : halfEdgeRing(vh)) {
                if (!valid(heh)) {
                    std::cout << "vertex " << vh << " connected to invalid halfedge " << heh
                              << std::endl;
                }
            }
            for (FaceHandle fh : faceRing(vh)) {
                if (!valid(fh)) {
                    std::cout << "vertex " << vh << " connected to invalid face " << fh
                              << std::endl;
                }
            }
            for (EdgeHandle eh : edgeRing(vh)) {
                if (!valid(eh)) {
                    std::cout << "vertex " << vh << " connected to invalid edge " << eh
                              << std::endl;
                }
            }
        }
    }
};

} // namespace Pvl

#include "Graph.inl.hpp"
