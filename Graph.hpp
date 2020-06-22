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
class Graph;

using HalfEdgeHandle = Handle<HalfEdge>;
using VertexHandle = Handle<Vertex>;
using FaceHandle = Handle<Face>;

struct HalfEdge {
    using Handle = HalfEdgeHandle;

    VertexHandle to;
    FaceHandle left;
    HalfEdgeHandle opposite;
    HalfEdgeHandle next;
    HalfEdgeHandle prev;


    HalfEdge(VertexHandle to)
        : to(to) {}

    bool boundary() const {
        return opposite < 0;
    }
};

struct Vertex {
    using Handle = VertexHandle;

    HalfEdgeHandle edge; // emanating vertex

    Vertex() {
        edge = HalfEdgeHandle(-1);
    }

    Vertex(HalfEdgeHandle eh)
        : edge(eh) {}

    bool valid() const {
        return edge >= 0;
    }

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

        HalfEdgeHandle operator*() const;
    };

    class VertexIterator : public IteratorBase {
    public:
        using IteratorBase::IteratorBase;

        VertexHandle operator*() const;
    };

    class FaceIterator : public IteratorBase {
    public:
        using IteratorBase::IteratorBase;

        FaceHandle operator*() const;
    };


    using HalfEdgeRange = Range<HalfEdgeIterator>;
    using VertexRange = Range<VertexIterator>;
    using FaceRange = Range<FaceIterator>;
};

static_assert(sizeof(Vertex) == sizeof(Handle<Vertex>), "error");

struct Face {
    HalfEdgeHandle edge;

    Face() = default;

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

        HalfEdgeHandle operator*() const;
    };

    class VertexIterator : public IteratorBase {
    public:
        using IteratorBase::IteratorBase;

        VertexHandle operator*() const;
    };

    class FaceIterator : public IteratorBase {
    public:
        FaceIterator(const Graph& graph, FaceHandle fh);

        FaceIterator(const Graph& graph, FaceHandle fh, EndTag);

        FaceHandle operator*() const;

        FaceIterator& operator++();
    };


    using HalfEdgeRange = Range<HalfEdgeIterator>;
    using VertexRange = Range<VertexIterator>;
    using FaceRange = Range<FaceIterator>;
};

class Graph {
    friend class Vertex;
    friend class Face;

    std::vector<Vertex> vertices_;
    std::vector<Face> faces_;
    std::vector<HalfEdge> edges_;

    // maps vertex to the set of emanating boundary edges
    std::map<VertexHandle, std::set<HalfEdgeHandle>> boundaryEdges_;


public:
    HalfEdgeHandle next(HalfEdgeHandle eh) const {
        return edges_[eh].next;
    }
    HalfEdgeHandle prev(HalfEdgeHandle eh) const {
        return edges_[eh].prev;
    }
    HalfEdgeHandle opposite(HalfEdgeHandle eh) const {
        ASSERT(!boundary(eh));
        return edges_[eh].opposite;
    }
    bool boundary(HalfEdgeHandle eh) const {
        return edges_[eh].boundary();
    }
    bool boundary(VertexHandle vh) const {
        return boundary(vertices_[vh].edge);
    }
    VertexHandle to(HalfEdgeHandle eh) const {
        return edges_[eh].to;
    }
    VertexHandle from(HalfEdgeHandle eh) const {
        return edges_[prev(eh)].to;
    }
    FaceHandle left(HalfEdgeHandle eh) const {
        return edges_[eh].left;
    }
    FaceHandle right(HalfEdgeHandle eh) const {
        ASSERT(!edges_[eh].boundary());
        return edges_[opposite(eh)].left;
    }
    HalfEdgeHandle from(VertexHandle vh) const {
        return vertices_[vh].edge;
    }
    HalfEdgeHandle to(VertexHandle vh) const {
        return prev(from(vh));
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
            handle_ = THandle(handle_ + 1);
            return *this;
        }

        bool operator==(const HandleIterator& other) const {
            return handle_ == other.handle_;
        }

        bool operator!=(const HandleIterator& other) const {
            return handle_ != other.handle_;
        }
    };

    using HalfEdgeIterator = HandleIterator<HalfEdgeHandle>;
    using VertexIterator = HandleIterator<VertexHandle>;
    using FaceIterator = HandleIterator<FaceHandle>;

    using HalfEdgeRange = Range<HalfEdgeIterator>;
    using VertexRange = Range<VertexIterator>;
    using FaceRange = Range<FaceIterator>;

    /// \todo rename to Iterators and Circulators
    HalfEdgeRange halfEdgeRange() const {
        return { HalfEdgeHandle(0), HalfEdgeHandle(edges_.size()) };
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

    std::array<int, 3> faceIndices(FaceHandle fh) const {
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

        edges_.emplace_back(vh2);
        edges_.emplace_back(vh3);
        edges_.emplace_back(vh1);

        std::array<HalfEdgeHandle, 3> eh = {
            HalfEdgeHandle(edges_.size() - 3),
            HalfEdgeHandle(edges_.size() - 2),
            HalfEdgeHandle(edges_.size() - 1),
        };
        std::array<VertexHandle, 3> vh = { vh1, vh2, vh3 };

        // set face edge
        f.edge = eh[0];

        for (int i = 0; i < 3; ++i) {
            // set emanating edge
            if (!vertices_[vh[i]].valid()) {
                vertices_[vh[i]].edge = eh[i];
            }

            // set next and previous edge
            HalfEdge& e = edges_[eh[i]];
            e.left = fh;
            e.next = eh[(i + 1) % 3];
            e.prev = eh[(i + 2) % 3];
            e.opposite = HalfEdgeHandle(-1);
        }

        for (int i = 0; i < 3; ++i) {
            for (HalfEdgeHandle neh : halfEdgeRing(to(eh[i]))) {
                // std::cout << "Visiting " << neh << "\n";
                if (to(neh) == from(eh[i])) {
                    if (!edges_[neh].boundary()) {
                        std::cout << "Complex edge!\n";
                        break;
                    }
                    edges_[neh].opposite = eh[i];
                    edges_[eh[i]].opposite = neh;
                    break;
                }
            }

            if (edges_[eh[i]].boundary()) {
                // opposite not found by circulation
                std::set<HalfEdgeHandle> removed;
                for (HalfEdgeHandle neh : boundaryEdges_[to(eh[i])]) {
                    if (to(neh) == from(eh[i])) {
                        if (!edges_[neh].boundary()) {
                            std::cout << "Complex edge!\n";
                        }
                        edges_[neh].opposite = eh[i];
                        edges_[eh[i]].opposite = neh;
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
            if (edges_[eh[i]].boundary()) {
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


    std::size_t numVertices() const {
        return vertices_.size();
    }

    std::size_t numFaces() const {
        return faces_.size();
    }
}; // namespace Nevim

} // namespace Pvl

#include "Graph.inl.hpp"
