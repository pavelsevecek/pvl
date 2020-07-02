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

template <typename TObject>
class Handle {
    using Index = std::uint32_t;

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
        opposite = HalfEdgeHandle(-1);
    }

    bool boundary() const {
        return opposite == HalfEdgeHandle(-1);
    }
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

// template<typename Vertex, typename HalfEdge, typename Face>
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

    /*HalfEdgeHandle from(VertexHandle vh) const {
        PVL_ASSERT(valid(vh));
        return vertices_[vh].edge;
    }
    HalfEdgeHandle to(VertexHandle vh) const {
        PVL_ASSERT(valid(vh));
        return prev(from(vh));
    }*/
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
        /*  PVL_ASSERT(valid(vh1) && valid(vh2));
          if (vh1 > vh2) {
              std::swap(vh1, vh2);
          }
          HalfEdgeHandle eh = halfEdge(vh1, vh2);
          if (!valid(eh)) {
              // might still be a boundary edge
              eh = halfEdge(vh2, vh1);
              PVL_ASSERT(!valid(eh) || boundary(eh));
              if (boundary(eh)) {
                  return edge(eh);
              }
          }
          return edge(eh);*/
        PVL_ASSERT(valid(vh1) && valid(vh2));
        for (EdgeHandle eh : edgeRing(vh1)) {
            HalfEdgeHandle heh = halfEdge(eh);
            if ((from(heh) == vh1 && to(heh) == vh2) || (to(heh) == vh1 && from(heh) == vh2)) {
                return eh;
            }
        }
        return EdgeHandle(-1);
    }
    HalfEdgeHandle emanating(VertexHandle vh) const {
        PVL_ASSERT(valid(vh));
        return vertices_[vh].edge;
    }
    HalfEdgeHandle incoming(VertexHandle vh) const {
        return prev(emanating(vh));
    }
    bool valid(FaceHandle fh) const {
        PVL_ASSERT(fh < faces_.size());
        return fh != FaceHandle(-1) && faces_[fh].edge != HalfEdgeHandle(-1);
    }
    bool valid(HalfEdgeHandle heh) const {
        PVL_ASSERT(heh < halfEdges_.size());
        return heh != HalfEdgeHandle(-1) && halfEdges_[heh].left != FaceHandle(-1);
    }
    bool valid(VertexHandle vh) const {
        PVL_ASSERT(vh < vertices_.size());
        return vh != VertexHandle(-1) && vertices_[vh].edge != HalfEdgeHandle(-1);
    }
    bool valid(EdgeHandle eh) const {
        PVL_ASSERT(eh < halfEdges_.size());
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
            return eh_ >= graph_.halfEdges_.size();
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
        return { Vertex::HalfEdgeIterator(*this, vh), Vertex::HalfEdgeIterator(*this, vh, Vertex::EndTag{}) };
    }

    Vertex::EdgeRange edgeRing(VertexHandle vh) const {
        return { Vertex::EdgeIterator(*this, vh), Vertex::EdgeIterator(*this, vh, Vertex::EndTag{}) };
    }

    Vertex::VertexRange vertexRing(VertexHandle vh) const {
        return { Vertex::VertexIterator(*this, vh), Vertex::VertexIterator(*this, vh, Vertex::EndTag{}) };
    }

    Vertex::FaceRange faceRing(VertexHandle vh) const {
        return { Vertex::FaceIterator(*this, vh), Vertex::FaceIterator(*this, vh, Vertex::EndTag{}) };
    }

    Face::HalfEdgeRange halfEdgeRing(FaceHandle fh) const {
        return { Face::HalfEdgeIterator(*this, fh), Face::HalfEdgeIterator(*this, fh, Face::EndTag{}) };
    }

    Face::VertexRange vertexRing(FaceHandle fh) const {
        return { Face::VertexIterator(*this, fh), Face::VertexIterator(*this, fh, Face::EndTag{}) };
    }

    Face::FaceRange faceRing(FaceHandle fh) const {
        return { Face::FaceIterator(*this, fh), Face::FaceIterator(*this, fh, Face::EndTag{}) };
    }

    std::array<VertexHandle, 3> faceVertices(FaceHandle fh) const {
        PVL_ASSERT(valid(fh));
        std::array<VertexHandle, 3> vertices;
        Face::VertexRange ring = vertexRing(fh);
        std::copy(ring.begin(), ring.end(), vertices.begin());
        return vertices;
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
            int safetyCntr = 0;
            for (HalfEdgeHandle neh : halfEdgeRing(to(eh[i]))) {
                // std::cout << "Visiting " << neh << "\n";
                if (to(neh) == from(eh[i])) {
                    if (!halfEdges_[neh].boundary()) {
                        std::cout << "Complex edge! when setting opposites" << std::endl;
                        break;
                    }
                    halfEdges_[neh].opposite = eh[i];
                    halfEdges_[eh[i]].opposite = neh;
                    break;
                }

                if (++safetyCntr > 50) {
                    std::cout << "Iterated 50times without finding the end, terminating" << std::endl;
                    break;
                }
            }

            if (halfEdges_[eh[i]].boundary()) {
                // opposite not found by circulation
                std::set<HalfEdgeHandle> removed;
                int safetyCntr = 0;

                for (HalfEdgeHandle neh : boundaryEdges_[to(eh[i])]) {
                    if (to(neh) == from(eh[i])) {
                        if (!halfEdges_[neh].boundary()) {
                            std::cout << "Complex edge! when finding boundary" << std::endl;
                            break;
                        }
                        halfEdges_[neh].opposite = eh[i];
                        halfEdges_[eh[i]].opposite = neh;
                        removed.insert(neh);
                        break;
                    }
                    if (++safetyCntr > 50) {
                        std::cout << "Iterated 50times without finding the end, terminating" << std::endl;
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
                HalfEdgeHandle e = vertices_[vh[i]].edge;
                HalfEdgeHandle e0 = e;
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
        HalfEdgeHandle edge;
        VertexHandle removed;
        VertexHandle remaining;
        FaceHandle left;  // left of the halfedge
        FaceHandle right; // right of the halfedge

        CollapseContext(const Graph& graph, const EdgeHandle eh)
            : edge(graph.halfEdge(eh)) {
            removed = graph.to(edge);
            remaining = graph.from(edge);
            left = graph.left(edge);
            if (!graph.boundary(edge)) {
                right = graph.right(edge);
            } else {
                right = FaceHandle(-1);
            }
        }
    };


    // collapse to to from
    void collapse(EdgeHandle edge) {
        PVL_ASSERT(collapseAllowed(edge));
        // check();
        CollapseContext context(*this, edge);

        HalfEdgeHandle heh = context.edge;
        bool onesided = boundary(heh);
        HalfEdgeHandle oheh = !onesided ? opposite(heh) : HalfEdgeHandle(-1);

        VertexHandle v0 = context.remaining;
        VertexHandle v1 = context.removed;
        VertexHandle vL = to(next(heh));
        VertexHandle vR = !onesided ? to(next(oheh)) : VertexHandle(-1);

        /*if (vertices_[context.remaining].edge == heh) {
            HalfEdgeHandle l = opposite(prev(heh));
            std::cout << "Moving emanating edge for vertex " << context.remaining << " to "
                      << from(l) << "-" << to(l) << std::endl;
            vertices_[context.remaining].edge = opposite(prev(heh));
        }*/

        HalfEdgeHandle ev0, evL, evR;
        /// \todo fix boundary case
        PVL_ASSERT(!boundary(next(heh)));
        if (vertices_[v1].edge != next(heh) && vertices_[v1].edge != oheh) {
            // unless the halfedge gets removed, simply use the emanating edge of the vertex to
            // be removed
            ev0 = vertices_[v1].edge;
        } else {
            PVL_ASSERT(!onesided);
            ev0 = opposite(prev(heh));
        }
        PVL_ASSERT(valid(ev0));

        if (vertices_[vL].edge == prev(heh)) {
            // move to any other halfedge, it cannot be a boundary
            PVL_ASSERT(!boundary(prev(heh)));
            evL = next(opposite(prev(heh)));
        } else {
            // keep the edge as it might be boundary
            evL = vertices_[vL].edge;
        }
        PVL_ASSERT(valid(evL));

        if (!onesided) {
            if (vertices_[vR].edge == prev(oheh)) {
                PVL_ASSERT(!boundary(prev(oheh)));
                evR = opposite(next(oheh));
            } else {
                evR = vertices_[vR].edge;
            }
            PVL_ASSERT(valid(evR));
        }


        /*for (HalfEdgeHandle neh : halfEdgeRing(to(heh))) {
            if (to(neh) == from(heh)) {
                continue;
            }
            halfEdges_[prev(neh)].to = from(heh);
        }*/
        HalfEdgeHandle neh = heh;
        do {
            neh = next(neh);
            if (boundary(neh)) {
                PVL_ASSERT(boundary(heh));
                break;
            }
            neh = opposite(neh);

            /*std::cout << "Reassigning 'to' for edge " << from(neh) << "-" << to(neh)
                      << " from " << halfEdges_[neh].to << " to " << from(heh) << std::endl;*/
            halfEdges_[neh].to = v0;
        } while (neh != heh);
        /*edges_[prev(eh)].next = next(opposite(next(eh)));
        edges_[prev(opposite(next(eh)))].next = from(eh);
        edges_[next(opposite(prev(oeh)))].next = next(oeh);*/
        connect(opposite(next(heh)), opposite(prev(heh)));
        if (!onesided) {
            connect(opposite(next(oheh)), opposite(prev(oheh)));
        }

        vertices_[v1].edge = HalfEdgeHandle(-1);
        vertices_[v0].edge = ev0;
        vertices_[vL].edge = evL;

        remove(context.left);
        if (!onesided) {
            vertices_[vR].edge = evR;
            remove(context.right);
        }

        PVL_ASSERT(valid(v0));
        PVL_ASSERT(valid(vL));
        PVL_ASSERT(onesided || valid(vR));
        PVL_ASSERT(from(vertices_[v0].edge) == v0);
        PVL_ASSERT(to(vertices_[v0].edge) != v1);
        PVL_ASSERT(from(vertices_[vL].edge) == vL);
        PVL_ASSERT(to(vertices_[vL].edge) != v1);
        if (!onesided) {
            PVL_ASSERT(from(vertices_[vR].edge) == vR);
            PVL_ASSERT(to(vertices_[vR].edge) != v1);
        }

        // check();
        /*remove(opposite(next(heh)));
        remove(prev(opposite(heh)));*/
        //   remove(heh);
        // remove(to(eh));
    }

    bool collapseAllowed(EdgeHandle edge) {
        PVL_ASSERT(valid(edge));
        HalfEdgeHandle eh = halfEdge(edge);

        /*if (boundary(eh)) { /// \todo enable boundary collapse
            return false;
        }*/

        if (boundary(to(eh)) != boundary(from(eh))) {
            // disallow contraction of boundary vertex and inner vertex
            return false;
        }
        if (boundary(next(eh)) || boundary(prev(eh))) {
            // disallow contraction of triangle with >1 boundary edge
            return false;
        }
        std::set<VertexHandle> ring1;
        for (VertexHandle vh : vertexRing(to(eh))) {
            if (vh != from(eh)) {
                // std::cout << "ring1 " << vh << std::endl;
                ring1.insert(vh);
            }
            /*if (!boundary(eh)) {
                for (HalfEdgeHandle neh : halfEdgeRing(vh)) {
                    if (boundary(neh)) {
                        return false;
                    }
                }
            }*/
        }
        std::set<VertexHandle> ring2;
        for (VertexHandle vh : vertexRing(from(eh))) {
            if (vh != to(eh)) {
                // std::cout << "ring2 " << vh << std::endl;
                ring2.insert(vh);
            }
            /*if (!boundary(eh)) {
                for (HalfEdgeHandle neh : halfEdgeRing(vh)) {
                    if (boundary(neh)) {
                        return false;
                    }
                }
            }*/
        }
        std::vector<VertexHandle> is;
        std::set_intersection(ring1.begin(), ring1.end(), ring2.begin(), ring2.end(), std::back_inserter(is));
        VertexHandle vl = to(next(eh));
        VertexHandle vr = !boundary(eh) ? to(next(opposite(eh))) : VertexHandle(-1);
        return std::all_of(is.begin(), is.end(), [vl, vr](VertexHandle vh) { return vh == vl || vh == vr; });
        // return is.size() == 2;
        // if (is.size() != 2) {
        //  return false;
        //}
        // VertexHandle vl = to(next(eh));
        // VertexHandle vr = to(next(opposite(eh)));
        // return ((is[0] == vl && is[1] == vr) || (is[1] == vl && is[0] == vr));
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
        return halfEdges_[heh].left == HalfEdgeHandle(-1);
    }

    bool removed(EdgeHandle eh) {
        return removed(halfEdge(eh));
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
                    std::cout << "vertex " << vh << " connected to invalid vertex " << nvh << std::endl;
                }
            }
            for (HalfEdgeHandle heh : halfEdgeRing(vh)) {
                if (!valid(heh)) {
                    std::cout << "vertex " << vh << " connected to invalid halfedge " << heh << std::endl;
                }
            }
            for (FaceHandle fh : faceRing(vh)) {
                if (!valid(fh)) {
                    std::cout << "vertex " << vh << " connected to invalid face " << fh << std::endl;
                }
            }
            for (EdgeHandle eh : edgeRing(vh)) {
                if (!valid(eh)) {
                    std::cout << "vertex " << vh << " connected to invalid edge " << eh << std::endl;
                }
            }
        }
    }
};

} // namespace Pvl

#include "Graph.inl.hpp"
