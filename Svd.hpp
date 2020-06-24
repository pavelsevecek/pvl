#pragma once

#include "Matrix.hpp"

namespace Pvl {

template <typename T>
struct Svd {
    Matrix<T, 3, 3> U;
    Vector<T, 3> S;
    Matrix<T, 3, 3> V;
};


inline double SIGN(double a, double b) {
    return b >= 0.0 ? fabs(a) : -fabs(a);
}

inline double PYTHAG(double a, double b) {
    double at = fabs(a), bt = fabs(b), ct, result;

    if (at > bt) {
        ct = bt / at;
        // ASSERT(isReal(ct));
        result = at * sqrt(1.0 + ct * ct);
    } else if (bt > 0.0) {
        ct = at / bt;
        // ASSERT(isReal(ct));
        result = bt * sqrt(1.0 + ct * ct);
    } else {
        result = 0.0;
    }
    return result;
}


int dsvd(float (&a)[3][3], float* w, float (&v)[3][3]) {
    const int m = 3, n = 3;
    int flag, i, its, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;
    double* rv1;

    PVL_ASSERT(m >= n, "#rows must be > #cols");

    rv1 = (double*)malloc((unsigned int)n * sizeof(double));

    /* Householder reduction to bidiagonal form */
    for (i = 0; i < n; i++) {
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m) {
            for (k = i; k < m; k++)
                scale += fabs((double)a[k][i]);
            if (scale) {
                for (k = i; k < m; k++) {
                    a[k][i] = (float)((double)a[k][i] / scale);
                    s += ((double)a[k][i] * (double)a[k][i]);
                }
                PVL_ASSERT(s >= 0., s);
                f = (double)a[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][i] = (float)(f - g);
                if (i != n - 1) {
                    for (j = l; j < n; j++) {
                        for (s = 0.0, k = i; k < m; k++)
                            s += ((double)a[k][i] * (double)a[k][j]);
                        f = s / h;
                        for (k = i; k < m; k++)
                            a[k][j] += (float)(f * (double)a[k][i]);
                    }
                }
                for (k = i; k < m; k++)
                    a[k][i] = (float)((double)a[k][i] * scale);
            }
        }
        w[i] = (float)(scale * g);

        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != n - 1) {
            for (k = l; k < n; k++)
                scale += fabs((double)a[i][k]);
            if (scale) {
                for (k = l; k < n; k++) {
                    a[i][k] = (float)((double)a[i][k] / scale);
                    s += ((double)a[i][k] * (double)a[i][k]);
                }
                PVL_ASSERT(s >= 0., s);
                f = (double)a[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][l] = (float)(f - g);
                for (k = l; k < n; k++)
                    rv1[k] = (double)a[i][k] / h;
                if (i != m - 1) {
                    for (j = l; j < m; j++) {
                        for (s = 0.0, k = l; k < n; k++)
                            s += ((double)a[j][k] * (double)a[i][k]);
                        for (k = l; k < n; k++)
                            a[j][k] += (float)(s * rv1[k]);
                    }
                }
                for (k = l; k < n; k++)
                    a[i][k] = (float)((double)a[i][k] * scale);
            }
        }
        anorm = std::max(anorm, (fabs((double)w[i]) + fabs(rv1[i])));
    }

    /* accumulate the right-hand transformation */
    for (i = n - 1; i >= 0; i--) {
        if (i < n - 1) {
            if (g) {
                for (j = l; j < n; j++)
                    v[j][i] = (float)(((double)a[i][j] / (double)a[i][l]) / g);
                /* double division to avoid underflow */
                for (j = l; j < n; j++) {
                    for (s = 0.0, k = l; k < n; k++)
                        s += ((double)a[i][k] * (double)v[k][j]);
                    for (k = l; k < n; k++)
                        v[k][j] += (float)(s * (double)v[k][i]);
                }
            }
            for (j = l; j < n; j++)
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }

    /* accumulate the left-hand transformation */
    for (i = n - 1; i >= 0; i--) {
        l = i + 1;
        g = (double)w[i];
        if (i < n - 1)
            for (j = l; j < n; j++)
                a[i][j] = 0.0;
        if (g) {
            g = 1.0 / g;
            if (i != n - 1) {
                for (j = l; j < n; j++) {
                    for (s = 0.0, k = l; k < m; k++)
                        s += ((double)a[k][i] * (double)a[k][j]);
                    f = (s / (double)a[i][i]) * g;
                    for (k = i; k < m; k++)
                        a[k][j] += (float)(f * (double)a[k][i]);
                }
            }
            for (j = i; j < m; j++)
                a[j][i] = (float)((double)a[j][i] * g);
        } else {
            for (j = i; j < m; j++)
                a[j][i] = 0.0;
        }
        ++a[i][i];
    }

    /* diagonalize the bidiagonal form */
    for (k = n - 1; k >= 0; k--) {       /* loop over singular values */
        for (its = 0; its < 30; its++) { /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--) { /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm) {
                    flag = 0;
                    break;
                }
                if (fabs((double)w[nm]) + anorm == anorm)
                    break;
            }
            if (flag) {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++) {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm) {
                        g = (double)w[i];
                        h = PYTHAG(f, g);
                        w[i] = (float)h;
                        h = 1.0 / h;
                        c = g * h;
                        s = (-f * h);
                        for (j = 0; j < m; j++) {
                            y = (double)a[j][nm];
                            z = (double)a[j][i];
                            a[j][nm] = (float)(y * c + z * s);
                            a[j][i] = (float)(z * c - y * s);
                        }
                    }
                }
            }
            z = (double)w[k];
            if (l == k) {      /* convergence */
                if (z < 0.0) { /* make singular value nonnegative */
                    w[k] = (float)(-z);
                    for (j = 0; j < n; j++)
                        v[j][k] = (-v[j][k]);
                }
                break;
            }
            PVL_ASSERT(its < 30, "No convergence after 30,000! iterations ");


            /* shift from bottom 2 x 2 minor */
            x = (double)w[l];
            nm = k - 1;
            y = (double)w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = PYTHAG(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++) {
                i = j + 1;
                g = rv1[i];
                y = (double)w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < n; jj++) {
                    x = (double)v[jj][j];
                    z = (double)v[jj][i];
                    v[jj][j] = (float)(x * c + z * s);
                    v[jj][i] = (float)(z * c - x * s);
                }
                z = PYTHAG(f, h);
                w[j] = (float)z;
                if (z) {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++) {
                    y = (double)a[jj][j];
                    z = (double)a[jj][i];
                    a[jj][j] = (float)(y * c + z * s);
                    a[jj][i] = (float)(z * c - y * s);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = (float)x;
        }
    }
    free((void*)rv1);
    return (1);
}


template <typename T>
Svd<T> singularValueDecomposition(const Matrix<T, 3, 3>& t) {
    float V[3][3];
    float S[3];

    // clang-format off
    float U[3][3] = {
        { float(t(0, 0)), float(t(1, 0)), float(t(2, 0)) },
        { float(t(0, 1)), float(t(1, 1)), float(t(2, 1)) },
        { float(t(0, 2)), float(t(1, 2)), float(t(2, 2)) },
    };

    dsvd(U, S, V);

    Svd<T> result;
    result.S = Vector<T, 3>(S[0], S[1], S[2]);
    result.U = Matrix<T, 3, 3>(Vector<T, 3>(U[0][0], U[0][1], U[0][2]),
                               Vector<T, 3>(U[1][0], U[1][1], U[1][2]),
                               Vector<T, 3>(U[2][0], U[2][1], U[2][2]));
    result.V = Matrix<T, 3, 3>(Vector<T, 3>(V[0][0], V[0][1], V[0][2]),
                               Vector<T, 3>(V[1][0], V[1][1], V[1][2]),
                               Vector<T, 3>(V[2][0], V[2][1], V[2][2]));
    // clang-format on
    return result;
}

} // namespace Pvl
