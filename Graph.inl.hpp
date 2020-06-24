#pragma once

#include "Graph.hpp"

namespace Pvl {

inline Vertex::IteratorBase::IteratorBase(const Graph& graph, VertexHandle vh)
    : graph_(graph) {
    eh_ = graph_.vertices_[vh].edge;
    if (eh_ != HalfEdgeHandle(-1)) {
        status_ = Status::BEGIN;
    } else {
        status_ = Status::END; // unreferenced vertex
    }
}

inline Vertex::IteratorBase::IteratorBase(const Graph& graph, VertexHandle vh, EndTag)
    : graph_(graph) {
    eh_ = graph_.vertices_[vh].edge;
    status_ = Status::END;
}

inline Vertex::IteratorBase& Vertex::IteratorBase::operator++() {
    const HalfEdgeHandle prev = graph_.prev(eh_);
    if (!graph_.boundary(prev)) {
        eh_ = graph_.opposite(prev);
        status_ = Status::ITERS;
    } else {
        status_ = Status::END;
    }
    PVL_ASSERT(eh_ != HalfEdgeHandle(-1));
    return *this;
}

inline bool Vertex::IteratorBase::operator!=(const IteratorBase& other) const {
    /// \todo proper
    // std::cout << "Comparing " << eh_ << " and " << other.eh_ << std::endl;
    return (status_ != Status::END || other.status_ != Status::END) &&
           (eh_ != other.eh_ || status_ == Status::BEGIN);
}

inline HalfEdgeHandle Vertex::HalfEdgeIterator::operator*() const {
    return eh_;
}

inline VertexHandle Vertex::VertexIterator::operator*() const {
    return graph_.to(eh_);
}

inline FaceHandle Vertex::FaceIterator::operator*() const {
    return graph_.left(eh_);
}

inline Vertex::EdgeIterator& Vertex::EdgeIterator::operator++() {
    if (status_ == Status::ITERS && graph_.boundary(eh_)) {
        // already reached the end
        status_ = Status::END;
    } else {
        const HalfEdgeHandle prev = graph_.prev(eh_);
        if (!graph_.boundary(prev)) {
            eh_ = graph_.opposite(prev);
        } else {
            eh_ = prev;
        }
        status_ = Status::ITERS;
    }
    PVL_ASSERT(eh_ != HalfEdgeHandle(-1));
    return *this;
}

inline EdgeHandle Vertex::EdgeIterator::operator*() const {
    return graph_.edge(eh_);
}


inline Face::IteratorBase::IteratorBase(const Graph& graph, FaceHandle fh)
    : graph_(graph) {
    eh_ = graph_.faces_[fh].edge;
    steps_ = 0;
}

inline Face::IteratorBase::IteratorBase(const Graph& graph, FaceHandle fh, EndTag)
    : graph_(graph) {
    eh_ = graph_.faces_[fh].edge;
    steps_ = 3;
}

inline Face::IteratorBase& Face::IteratorBase::operator++() {
    eh_ = graph_.next(eh_);
    ++steps_;
    return *this;
}

inline bool Face::IteratorBase::operator!=(const IteratorBase& other) const {
    return eh_ != other.eh_ || steps_ != other.steps_;
}

inline HalfEdgeHandle Face::HalfEdgeIterator::operator*() const {
    return eh_;
}

inline VertexHandle Face::VertexIterator::operator*() const {
    return graph_.to(eh_);
}

inline Face::FaceIterator::FaceIterator(const Graph& graph, FaceHandle fh)
    : IteratorBase(graph, fh) {
    while (graph_.boundary(eh_) && steps_ < 3) {
        eh_ = graph_.next(eh_);
        ++steps_;
    }
}

inline Face::FaceIterator::FaceIterator(const Graph& graph, FaceHandle fh, EndTag t)
    : IteratorBase(graph, fh, t) {}

inline FaceHandle Face::FaceIterator::operator*() const {
    return graph_.right(eh_);
}

inline Face::FaceIterator& Face::FaceIterator::operator++() {
    do {
        eh_ = graph_.next(eh_);
        ++steps_;
    } while (graph_.boundary(eh_) && steps_ < 3);
    return *this;
}

} // namespace Pvl
