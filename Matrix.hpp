#pragma once

#include "Vector.hpp"

namespace Pvl {

template <typename T, int Cols, int Rows>
class Matrix {
    std::array<Vector<T, Cols>, Rows> rows_;

public:
    Matrix() = default;

    Matrix(const Vector<T, Cols>& r1, const Vector<T, Cols>& r2, const Vector<T, Cols>& r3)
        : rows_{ r1, r2, r3 } {}

    Matrix(const Vector<T, Cols>& r1,
        const Vector<T, Cols>& r2,
        const Vector<T, Cols>& r3,
        const Vector<T, Cols>& r4)
        : rows_{ r1, r2, r3, r4 } {}

    T& operator()(const int c, const int r) {
        PVL_ASSERT(r < Rows);
        PVL_ASSERT(c < Cols);
        return rows_[r][c];
    }

    const T& operator()(const int c, const int r) const {
        PVL_ASSERT(r < Rows);
        PVL_ASSERT(c < Cols);
        return rows_[r][c];
    }

    const Vector<T, Cols>& row(const int idx) const {
        PVL_ASSERT(idx < Rows);
        return rows_[idx];
    }

    Vector<T, Cols>& row(const int idx) {
        PVL_ASSERT(idx < Rows);
        return rows_[idx];
    }

    Vector<T, Rows> column(const int idx) const {
        PVL_ASSERT(idx < Cols);
        Vector<T, Rows> col;
        for (int i = 0; i < Rows; ++i) {
            col[i] = row(i)[idx];
        }
        return col;
    }

    Vector<T, Rows> transform(const Vector<T, Cols>& v) const {
        Vector<T, Rows> res;
        for (int i = 0; i < Rows; ++i) {
            res[i] = dot(row(i), v);
        }
        return res;
    }

    static Matrix null() {
        return Matrix(Vector<T, Cols>(0), Vector<T, Cols>(0), Vector<T, Cols>(0));
    }

    static Matrix identity() {
        Matrix m;
        for (int j = 0; j < Rows; ++j) {
            for (int i = 0; i < Cols; ++i) {
                m(i, j) = int(i == j);
            }
        }
        return m;
    }

    Matrix operator+(const Matrix& other) const {
        Matrix m;
        for (int j = 0; j < Rows; ++j) {
            for (int i = 0; i < Cols; ++i) {
                m(i, j) = (*this)(i, j) + other(i, j);
            }
        }
        return m;
    }

    Matrix operator-(const Matrix& other) const {
        Matrix m;
        for (int j = 0; j < Rows; ++j) {
            for (int i = 0; i < Cols; ++i) {
                m(i, j) = (*this)(i, j) - other(i, j);
            }
        }
        return m;
    }

    Matrix operator*(const T f) const {
        Matrix m;
        for (int j = 0; j < Rows; ++j) {
            for (int i = 0; i < Cols; ++i) {
                m(i, j) = (*this)(i, j) * f;
            }
        }
        return m;
    }

    Matrix operator/(const T f) const {
        PVL_ASSERT(f != 0);
        Matrix m;
        for (int j = 0; j < Rows; ++j) {
            for (int i = 0; i < Cols; ++i) {
                m(i, j) = (*this)(i, j) / f;
            }
        }
        return m;
    }

    Matrix& operator+=(const Matrix& other) {
        *this = *this + other;
        return *this;
    }

    Matrix& operator-=(const Matrix& other) {
        *this = *this - other;
        return *this;
    }

    template <int N>
    Matrix<T, Rows, N> operator*(const Matrix<T, N, Cols>& other) const {
        Matrix<T, Rows, N> result;
        for (int r = 0; r < N; ++r) {
            for (int c = 0; c < Rows; ++c) {
                result(c, r) = dot(row(r), other.column(c));
            }
        }
        return result;
    }
};

using Mat33f = Matrix<float, 3, 3>;
using Mat44f = Matrix<float, 4, 4>;

template <typename T, int N, int M>
Matrix<T, N, M> outerProd(const Vector<T, N>& v1, const Vector<T, M>& v2) {
    Matrix<T, N, M> res;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            res(i, j) = v1[i] * v2[j];
        }
    }
    return res;
}

template <typename T, int N>
inline T trace(const Matrix<T, N, N>& m) {
    T tr = 0.;
    for (int i = 0; i < N; ++i) {
        tr += m(i, i);
    }
    return tr;
}

inline Mat33f invert(const Mat33f& m) {
    double det = m(0, 0) * (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)) -
                 m(0, 1) * (m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0)) +
                 m(0, 2) * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0));

    double invdet = 1 / det;

    Mat33f minv; // inverse of matrix m
    minv(0, 0) = (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)) * invdet;
    minv(0, 1) = (m(0, 2) * m(2, 1) - m(0, 1) * m(2, 2)) * invdet;
    minv(0, 2) = (m(0, 1) * m(1, 2) - m(0, 2) * m(1, 1)) * invdet;
    minv(1, 0) = (m(1, 2) * m(2, 0) - m(1, 0) * m(2, 2)) * invdet;
    minv(1, 1) = (m(0, 0) * m(2, 2) - m(0, 2) * m(2, 0)) * invdet;
    minv(1, 2) = (m(1, 0) * m(0, 2) - m(0, 0) * m(1, 2)) * invdet;
    minv(2, 0) = (m(1, 0) * m(2, 1) - m(2, 0) * m(1, 1)) * invdet;
    minv(2, 1) = (m(2, 0) * m(0, 1) - m(0, 0) * m(2, 1)) * invdet;
    minv(2, 2) = (m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1)) * invdet;
    return minv;
}

inline float determinant(const Mat44f& m) {
    Mat44f minv; /// \todo
    minv(0, 0) = m(1, 1) * m(2, 2) * m(3, 3) - m(1, 1) * m(3, 2) * m(2, 3) -
                 m(1, 2) * m(2, 1) * m(3, 3) + m(1, 2) * m(3, 1) * m(2, 3) +
                 m(1, 3) * m(2, 1) * m(3, 2) - m(1, 3) * m(3, 1) * m(2, 2);

    minv(0, 1) = -m(0, 1) * m(2, 2) * m(3, 3) + m(0, 1) * m(3, 2) * m(2, 3) +
                 m(0, 2) * m(2, 1) * m(3, 3) - m(0, 2) * m(3, 1) * m(2, 3) -
                 m(0, 3) * m(2, 1) * m(3, 2) + m(0, 3) * m(3, 1) * m(2, 2);

    minv(0, 2) = m(0, 1) * m(1, 2) * m(3, 3) - m(0, 1) * m(3, 2) * m(1, 3) -
                 m(0, 2) * m(1, 1) * m(3, 3) + m(0, 2) * m(3, 1) * m(1, 3) +
                 m(0, 3) * m(1, 1) * m(3, 2) - m(0, 3) * m(3, 1) * m(1, 2);

    minv(0, 3) = -m(0, 1) * m(1, 2) * m(2, 3) + m(0, 1) * m(2, 2) * m(1, 3) +
                 m(0, 2) * m(1, 1) * m(2, 3) - m(0, 2) * m(2, 1) * m(1, 3) -
                 m(0, 3) * m(1, 1) * m(2, 2) + m(0, 3) * m(2, 1) * m(1, 2);

    const float det = m(0, 0) * minv(0, 0) + m(1, 0) * minv(0, 1) + m(2, 0) * minv(0, 2) +
                      m(3, 0) * minv(0, 3);
    return det;
}


inline Mat44f invert(const Mat44f& m) {
    Mat44f minv;
    minv(0, 0) = m(1, 1) * m(2, 2) * m(3, 3) - m(1, 1) * m(3, 2) * m(2, 3) -
                 m(1, 2) * m(2, 1) * m(3, 3) + m(1, 2) * m(3, 1) * m(2, 3) +
                 m(1, 3) * m(2, 1) * m(3, 2) - m(1, 3) * m(3, 1) * m(2, 2);

    minv(0, 1) = -m(0, 1) * m(2, 2) * m(3, 3) + m(0, 1) * m(3, 2) * m(2, 3) +
                 m(0, 2) * m(2, 1) * m(3, 3) - m(0, 2) * m(3, 1) * m(2, 3) -
                 m(0, 3) * m(2, 1) * m(3, 2) + m(0, 3) * m(3, 1) * m(2, 2);

    minv(0, 2) = m(0, 1) * m(1, 2) * m(3, 3) - m(0, 1) * m(3, 2) * m(1, 3) -
                 m(0, 2) * m(1, 1) * m(3, 3) + m(0, 2) * m(3, 1) * m(1, 3) +
                 m(0, 3) * m(1, 1) * m(3, 2) - m(0, 3) * m(3, 1) * m(1, 2);

    minv(0, 3) = -m(0, 1) * m(1, 2) * m(2, 3) + m(0, 1) * m(2, 2) * m(1, 3) +
                 m(0, 2) * m(1, 1) * m(2, 3) - m(0, 2) * m(2, 1) * m(1, 3) -
                 m(0, 3) * m(1, 1) * m(2, 2) + m(0, 3) * m(2, 1) * m(1, 2);

    minv(1, 0) = -m(1, 0) * m(2, 2) * m(3, 3) + m(1, 0) * m(3, 2) * m(2, 3) +
                 m(1, 2) * m(2, 0) * m(3, 3) - m(1, 2) * m(3, 0) * m(2, 3) -
                 m(1, 3) * m(2, 0) * m(3, 2) + m(1, 3) * m(3, 0) * m(2, 2);

    minv(1, 1) = m(0, 0) * m(2, 2) * m(3, 3) - m(0, 0) * m(3, 2) * m(2, 3) -
                 m(0, 2) * m(2, 0) * m(3, 3) + m(0, 2) * m(3, 0) * m(2, 3) +
                 m(0, 3) * m(2, 0) * m(3, 2) - m(0, 3) * m(3, 0) * m(2, 2);

    minv(1, 2) = -m(0, 0) * m(1, 2) * m(3, 3) + m(0, 0) * m(3, 2) * m(1, 3) +
                 m(0, 2) * m(1, 0) * m(3, 3) - m(0, 2) * m(3, 0) * m(1, 3) -
                 m(0, 3) * m(1, 0) * m(3, 2) + m(0, 3) * m(3, 0) * m(1, 2);

    minv(1, 3) = m(0, 0) * m(1, 2) * m(2, 3) - m(0, 0) * m(2, 2) * m(1, 3) -
                 m(0, 2) * m(1, 0) * m(2, 3) + m(0, 2) * m(2, 0) * m(1, 3) +
                 m(0, 3) * m(1, 0) * m(2, 2) - m(0, 3) * m(2, 0) * m(1, 2);

    minv(2, 0) = m(1, 0) * m(2, 1) * m(3, 3) - m(1, 0) * m(3, 1) * m(2, 3) -
                 m(1, 1) * m(2, 0) * m(3, 3) + m(1, 1) * m(3, 0) * m(2, 3) +
                 m(1, 3) * m(2, 0) * m(3, 1) - m(1, 3) * m(3, 0) * m(2, 1);

    minv(2, 1) = -m(0, 0) * m(2, 1) * m(3, 3) + m(0, 0) * m(3, 1) * m(2, 3) +
                 m(0, 1) * m(2, 0) * m(3, 3) - m(0, 1) * m(3, 0) * m(2, 3) -
                 m(0, 3) * m(2, 0) * m(3, 1) + m(0, 3) * m(3, 0) * m(2, 1);

    minv(2, 2) = m(0, 0) * m(1, 1) * m(3, 3) - m(0, 0) * m(3, 1) * m(1, 3) -
                 m(0, 1) * m(1, 0) * m(3, 3) + m(0, 1) * m(3, 0) * m(1, 3) +
                 m(0, 3) * m(1, 0) * m(3, 1) - m(0, 3) * m(3, 0) * m(1, 1);

    minv(2, 3) = -m(0, 0) * m(1, 1) * m(2, 3) + m(0, 0) * m(2, 1) * m(1, 3) +
                 m(0, 1) * m(1, 0) * m(2, 3) - m(0, 1) * m(2, 0) * m(1, 3) -
                 m(0, 3) * m(1, 0) * m(2, 1) + m(0, 3) * m(2, 0) * m(1, 1);

    minv(3, 0) = -m(1, 0) * m(2, 1) * m(3, 2) + m(1, 0) * m(3, 1) * m(2, 2) +
                 m(1, 1) * m(2, 0) * m(3, 2) - m(1, 1) * m(3, 0) * m(2, 2) -
                 m(1, 2) * m(2, 0) * m(3, 1) + m(1, 2) * m(3, 0) * m(2, 1);

    minv(3, 1) = m(0, 0) * m(2, 1) * m(3, 2) - m(0, 0) * m(3, 1) * m(2, 2) -
                 m(0, 1) * m(2, 0) * m(3, 2) + m(0, 1) * m(3, 0) * m(2, 2) +
                 m(0, 2) * m(2, 0) * m(3, 1) - m(0, 2) * m(3, 0) * m(2, 1);

    minv(3, 2) = -m(0, 0) * m(1, 1) * m(3, 2) + m(0, 0) * m(3, 1) * m(1, 2) +
                 m(0, 1) * m(1, 0) * m(3, 2) - m(0, 1) * m(3, 0) * m(1, 2) -
                 m(0, 2) * m(1, 0) * m(3, 1) + m(0, 2) * m(3, 0) * m(1, 1);

    minv(3, 3) = m(0, 0) * m(1, 1) * m(2, 2) - m(0, 0) * m(2, 1) * m(1, 2) -
                 m(0, 1) * m(1, 0) * m(2, 2) + m(0, 1) * m(2, 0) * m(1, 2) +
                 m(0, 2) * m(1, 0) * m(2, 1) - m(0, 2) * m(2, 0) * m(1, 1);
    const float det = m(0, 0) * minv(0, 0) + m(1, 0) * minv(0, 1) + m(2, 0) * minv(0, 2) +
                      m(3, 0) * minv(0, 3);
    const float norm = 1.f / det;

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            minv(i, j) *= norm;
        }
    }

    return minv;
}

inline Mat44f homogeneous(const Mat33f& m) {
    Mat44f res = Mat44f::identity();
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            res(i, j) = m(i, j);
        }
    }
    return res;
}

inline Mat44f getTranslationMatrix(const Vec3f& pos) {
    Mat44f m = Mat44f::identity();
    for (int i = 0; i < 3; ++i) {
        m(3, i) = pos[i];
    }
    return m;
}

inline Mat33f getRotationMatrix(const Vec3f& axis, const float angle) {
    const float u = axis[0];
    const float v = axis[1];
    const float w = axis[2];
    const float s = sin(angle);
    const float c = cos(angle);
    return {
        Vec3f(u * u + (v * v + w * w) * c, u * v * (1 - c) - w * s, u * w * (1 - c) + v * s),
        Vec3f(u * v * (1 - c) + w * s, v * v + (u * u + w * w) * c, v * w * (1 - c) - u * s),
        Vec3f(u * w * (1 - c) - v * s, v * w * (1 - c) + u * s, w * w + (u * u + v * v) * c),
    };
}

} // namespace Pvl
