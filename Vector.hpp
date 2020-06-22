#pragma once

#include "Assert.hpp"
#include "Math.hpp"
#include <array>

namespace Pvl {

template <typename T, int Dim>
class Vector {
    std::array<T, Dim> values_;

public:
    using Type = T;

    Vector() = default;

    Vector(const T value) {
        for (int i = 0; i < Dim; ++i) {
            values_[i] = value;
        }
    }

    Vector(const T x, const T y) /// \todo
        : values_{ x, y } {}

    Vector(const T x, const T y, const T z)
        : values_{ x, y, z } {}

    Vector(const T x, const T y, const T z, const T w)
        : values_{ x, y, z, w } {}

    T& operator[](const int idx) {
        ASSERT(unsigned(idx) < unsigned(Dim), idx, Dim);
        return values_[idx];
    }

    const T& operator[](const int idx) const {
        ASSERT(unsigned(idx) < unsigned(Dim), idx, Dim);
        return values_[idx];
    }

    Vector operator-() const {
        Vector res;
        for (int i = 0; i < Dim; ++i) {
            res[i] = -values_[i];
        }
        return res;
    }

    bool operator==(const Vector& other) const {
        for (int i = 0; i < Dim; ++i) {
            if (values_[i] != other[i]) {
                return false;
            }
        }
        return true;
    }

    bool operator!=(const Vector& other) const {
        return !(*this == other);
    }

    Vector& operator+=(const Vector& other) {
        *this = *this + other;
        return *this;
    }

    Vector& operator-=(const Vector& other) {
        *this = *this - other;
        return *this;
    }

    Vector operator+(const Vector& other) const {
        Vector res;
        for (int i = 0; i < Dim; ++i) {
            res[i] = values_[i] + other[i];
        }
        return res;
    }

    Vector operator-(const Vector& other) const {
        Vector res;
        for (int i = 0; i < Dim; ++i) {
            res[i] = values_[i] - other[i];
        }
        return res;
    }

    Vector operator*(const T f) const {
        Vector res;
        for (int i = 0; i < Dim; ++i) {
            res[i] = values_[i] * f;
        }
        return res;
    }

    Vector operator/(const T f) const {
        Vector res;
        for (int i = 0; i < Dim; ++i) {
            res[i] = values_[i] / f;
        }
        return res;
    }

    friend Vector operator*(const T f, const Vector& v) {
        return v * f;
    }

    static constexpr int size() {
        return Dim;
    }

    T* begin() {
        return values_.begin();
    }

    T* end() {
        return values_.end();
    }

    const T* begin() const {
        return values_.begin();
    }

    const T* end() const {
        return values_.end();
    }

    const T* data() const {
        return values_.data();
    }
}; // namespace Universe

using Vec2f = Vector<float, 2>;
using Vec3f = Vector<float, 3>;
using Vec4f = Vector<float, 4>;

using Vec2i = Vector<int, 2>;
using Vec3i = Vector<int, 3>;
using Vec4i = Vector<int, 4>;

template <typename T, int Dim>
T dot(const Vector<T, Dim>& v1, const Vector<T, Dim>& v2) {
    T res = 0.f;
    for (int i = 0; i < Dim; ++i) {
        res += v1[i] * v2[i];
    }
    return res;
}

template <typename T, int Dim>
T normSqr(const Vector<T, Dim>& v) {
    return dot(v, v);
}

template <typename T, int Dim>
T norm(const Vector<T, Dim>& v) {
    return sqrt(normSqr(v));
}

template <typename T, int Dim>
T normL1(const Vector<T, Dim>& v) {
    T sum = T(0);
    for (int i = 0; i < Dim; ++i) {
        sum += std::abs(v[i]);
    }
    return sum;
}

template <typename T, int Dim>
Vector<T, Dim> normalize(const Vector<T, Dim>& v) {
    return v / norm(v);
}

inline Vec3f cross(const Vec3f& v1, const Vec3f& v2) {
    return Vec3f(v1[1] * v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2], v1[0] * v2[1] - v1[1] * v2[0]);
}

template <typename T1, typename T2, int Dim>
Vector<T1, Dim> vectorCast(const Vector<T2, Dim>& v) {
    Vector<T1, Dim> res;
    for (int i = 0; i < Dim; ++i) {
        res[i] = static_cast<T1>(v[i]);
    }
    return res;
}

template <typename T, int Dim>
Vector<T, Dim> max(const Vector<T, Dim>& v1, const Vector<T, Dim>& v2) {
    Vector<T, Dim> res;
    for (int i = 0; i < Dim; ++i) {
        res[i] = std::max(v1[i], v2[i]);
    }
    return res;
}

template <typename T, int Dim>
Vector<T, Dim> min(const Vector<T, Dim>& v1, const Vector<T, Dim>& v2) {
    Vector<T, Dim> res;
    for (int i = 0; i < Dim; ++i) {
        res[i] = std::min(v1[i], v2[i]);
    }
    return res;
}

template <typename T, int Dim>
int argMax(const Vector<T, Dim>& v) {
    return int(std::max_element(v.begin(), v.end()) - v.begin());
}

template <typename T, int Dim>
int argMin(const Vector<T, Dim>& v) {
    return int(std::min_element(v.begin(), v.end()) - v.begin());
}

template <typename T, int Dim>
Vector<T, Dim> floor(const Vector<T, Dim>& v) {
    Vector<T, Dim> res;
    for (int i = 0; i < Dim; ++i) {
        res[i] = std::floor(v[i]);
    }
    return res;
}

template <typename T, int Dim>
Vector<T, Dim> round(const Vector<T, Dim>& v) {
    Vector<T, Dim> res;
    for (int i = 0; i < Dim; ++i) {
        res[i] = std::round(v[i]);
    }
    return res;
}

template <typename T, int Dim>
Vector<T, Dim> sqr(const Vector<T, Dim>& v) {
    Vector<T, Dim> res;
    for (int i = 0; i < Dim; ++i) {
        res[i] = Pvl::sqr(v[i]);
    }
    return res;
}

inline Vec4f homogeneous(const Vec3f& v) {
    return Vec4f(v[0], v[1], v[2], 1.f);
}

} // namespace Pvl
