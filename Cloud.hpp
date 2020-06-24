#pragma once

#include "Box.hpp"
#include <vector>

namespace Pvl {

template <typename Type>
class CloudProperty {
private:
    std::vector<Type> values_;

public:
};


template <typename Vec>
class PointInsideBoxPredicate {
    BoundingBox<Vec> box_;

public:
    bool operator()(const Vec& p) const {
        return box_.contains(p);
    }
};

template <typename Iterator, typename Predicate>
class SubsetIterator {
    Iterator iter_;
    Predicate predicate_;

public:
    SubsetIterator(const Predicate& predicate)
        : predicate_(predicate) {}

    SubsetIterator& operator++() {
        do {
            ++iter_;
        } while (!predicate_());
        return *this;
    }

    bool operator!=(const SubsetIterator& other) const {
        return iter_ != other.iter_;
    }
};

template <typename Predicate>
class CloudSubset {
public:
    /*SubsetIterator<Predicate> begin() {}

    SubsetIterator<Predicate> end() {}*/
};

template <typename... Properties>
class Cloud : public Properties... {};

// using MyCloud = Cloud<Velocity, Normals, Colors, Custom>;

} // namespace Pvl
