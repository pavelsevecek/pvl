#pragma once
#include "Vector.hpp"
#include <vector>

namespace Pvl {

template <typename Object, int Dim>
class GridIterator {
public:
    using Idxs = Vector<int, Dim>;

private:
    Object* ptr_;
    std::size_t index_;
    Idxs dims_;


public:
    GridIterator(Object* ptr, const std::size_t index, const Idxs& dims)
        : ptr_(ptr)
        , index_(index)
        , dims_(dims) {}

    GridIterator& operator++() {
        ptr_++;
        index_++;
        return *this;
    }
    Object& operator*() const {
        return *ptr_;
    }
    bool operator==(const GridIterator& other) const {
        return ptr_ == other.ptr_;
    }
    bool operator!=(const GridIterator& other) const {
        return ptr_ != other.ptr_;
    }

    Idxs idxs() const {
        Idxs c;
        std::size_t k = index_;
        for (int i = 0; i < Dim; ++i) {
            c[i] = k % dims_[i];
            k /= dims_[i];
        }
    }
};

template <typename Object, int Dim>
class UniformGrid {
public:
    // using Point = Point;
    using Idxs = Vector<int, Dim>;

public:
    std::vector<Object> data_;
    Idxs dims_;

public:
    UniformGrid() = default;

    UniformGrid(const Idxs& dims)
        : dims_(dims) {
        data_.resize(voxelCount());
    }

    Object& operator()(const Idxs& idxs) {
        return data_[map(idxs)];
    }

    const Object& operator()(const Idxs& idxs) const {
        return data_[map(idxs)];
    }

    std::size_t voxelCount() const {
        std::size_t cnt = 1;
        for (int i = 0; i < Dim; ++i) {
            cnt *= dims_[i];
        }
        return cnt;
    }

    GridIterator<Object, Dim> begin() {
        return { data_.data(), 0, dims_ };
    }

    GridIterator<Object, Dim> end() {
        return { data_.data() + data_.size(), 0, dims_ };
    }

    GridIterator<const Object, Dim> begin() const {
        return { data_.data(), 0, dims_ };
    }

    GridIterator<const Object, Dim> end() const {
        return { data_.data() + data_.size(), 0, dims_ };
    }

private:
    std::size_t map(const Idxs& idxs) const {
        std::size_t linear = 0;
        std::size_t stride = 1;
        for (int i = 0; i < Dim; ++i) {
            linear += idxs[i] * stride;
            stride *= dims_[i];
        }
        return linear;
    }
};


} // namespace Pvl
