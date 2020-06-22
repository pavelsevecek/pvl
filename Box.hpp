#pragma once

#include "Vector.hpp"

namespace Pvl {

template <typename T, int Dim>
class BoundingBox {
    Vector<T, Dim> lower_;
    Vector<T, Dim> upper_;

public:
    BoundingBox()
        : lower_(std::numeric_limits<T>::max())
        , upper_(std::numeric_limits<T>::lowest()) {}

    BoundingBox(const Vector<T, Dim>& lower, const Vector<T, Dim>& upper)
        : lower_(lower)
        , upper_(upper) {}

    Vector<T, Dim>& lower() {
        return lower_;
    }

    const Vector<T, Dim>& lower() const {
        return lower_;
    }

    Vector<T, Dim>& upper() {
        return upper_;
    }

    const Vector<T, Dim>& upper() const {
        return upper_;
    }

    Vector<T, Dim> size() const {
        return upper_ - lower_;
    }

    Vector<T, Dim> center() const {
        return T(0.5) * (upper_ + lower_);
    }


    bool contains(const Vector<T, Dim>& p) const {
        for (int i = 0; i < Dim; ++i) {
            if (p[i] < lower_[i] || p[i] > upper_[i]) {
                return false;
            }
        }
        return true;
    }

    void extend(const Vector<T, Dim>& p) {
        lower_ = min(lower_, p);
        upper_ = max(upper_, p);
    }
};

/// \brief Splits the box along given coordinate.
///
/// The splitting plane must pass through the box.
template <typename T, int Dim>
std::pair<BoundingBox<T, Dim>, BoundingBox<T, Dim>> splitBox(const BoundingBox<T, Dim>& box,
    const int dim,
    const T x) {
    /*ASSERT(isValid());*/
    ASSERT(dim < Dim);
    ASSERT(x >= box.lower()[dim] && x <= box.upper()[dim]);
    BoundingBox<T, Dim> b1 = box;
    BoundingBox<T, Dim> b2 = box;
    b1.lower()[dim] = x;
    b2.upper()[dim] = x;
    return std::make_pair(b1, b2);
}


} // namespace Pvl
