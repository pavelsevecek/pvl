#pragma once

#include "Math.hpp"
#include "Utils.hpp"
#include "Vector.hpp"

namespace Pvl {

/// \brief Base class for all SPH kernels.
///
/// Provides an interface for computing kernel values and gradients. All derived class must implement method
/// <code>valueImpl</code> and <code>gradImpl</code>. Both function take SQUARED value of dimensionless
/// distance q as a parameter. Function value returns the kernel value, grad returns gradient DIVIDED BY q.
template <typename TDerived, typename Float, int Dim>
class Kernel : public Noncopyable {
public:
    using Vec = Vector<Float, Dim>;

    Kernel() = default;

    /// Value of kernel at given point
    /// this should be called only once for a pair of particles as there is expensive division
    /// \todo Potentially precompute the 3rd power ...
    Float value(const Vec& r, const Float h) const noexcept {
        ASSERT(h > 0.f);
        const Float hInv = 1.f / h;
        return pow<Dim>(hInv) * impl().valueImpl(normSqr(r) * sqr(hInv));
    }

    Vec grad(const Vec& r, const Float h) const noexcept {
        ASSERT(h > 0.f);
        const Float hInv = 1.f / h;
        return r * pow<Dim + 2>(hInv) * impl().gradImpl(normSqr(r) * sqr(hInv));
    }

private:
    const TDerived& impl() const noexcept {
        return static_cast<const TDerived&>(*this);
    }
};


/// \brief A look-up table approximation of the kernel.
///
/// Can be constructed from any SPH kernel. Use this class exclusively for any high-performance computations,
/// it is always faster than using kernel functions directly (except for trivial kerneles, such as \ref
/// TriangleKernel). The precision difference is about 1.e-6.
template <typename Float, int Dim>
class LutKernel : public Kernel<LutKernel<Float, Dim>, Float, Dim> {
private:
    static constexpr int NEntries = 40000;

    /// \todo replace with StaticArray?
    Float values[NEntries];
    Float grads[NEntries];

    Float rad = 0.f;
    Float qSqrToIdx = 0.f;

public:
    LutKernel() {
        /// \todo initialize, otherwise compiler complains about using uninitialized values
        qSqrToIdx = NAN;
        for (Float& v : values) {
            v = NAN;
        }
        for (Float& g : grads) {
            g = NAN;
        }
    }

    LutKernel(const LutKernel& other) {
        *this = other;
    }

    LutKernel& operator=(const LutKernel& other) {
        rad = other.rad;
        qSqrToIdx = other.qSqrToIdx;
        for (int i = 0; i < NEntries; ++i) {
            values[i] = other.values[i];
            grads[i] = other.grads[i];
        }
        return *this;
    }

    /// \brief Constructs LUT kernel given an exact SPH kernel.
    template <typename TKernel,
        typename = std::enable_if_t<!std::is_same<std::decay_t<TKernel>, LutKernel>::value>>
    LutKernel(TKernel&& source) {
        rad = source.radius();

        ASSERT(rad > 0.f);
        const Float radInvSqr = 1.f / (rad * rad);
        qSqrToIdx = Float(NEntries) * radInvSqr;
        for (int i = 0; i < NEntries; ++i) {
            const Float qSqr = Float(i) / qSqrToIdx;
            values[i] = source.valueImpl(qSqr);
            grads[i] = source.gradImpl(qSqr);
        }
        /// \todo re-normalize?
    }

    bool initialized() const {
        return rad > 0.f;
    }

    Float radius() const noexcept {
        return rad;
    }

    Float valueImpl(const Float qSqr) const noexcept {
        ASSERT(qSqr >= 0.f);
        ASSERT(initialized());
        if (qSqr >= sqr(rad)) {
            // outside of kernel support
            return 0.f;
        }
        // linear interpolation of stored values
        const Float floatIdx = qSqrToIdx * qSqr;
        const int idx1 = Size(floatIdx);
        ASSERT(idx1 < NEntries);
        const int idx2 = idx1 + 1;
        const Float ratio = floatIdx - Float(idx1);
        ASSERT(ratio >= 0.f && ratio < 1.f);

        return values[idx1] * (1.f - ratio) + (int(idx2 < NEntries) * values[idx2]) * ratio;
    }

    Float gradImpl(const Float qSqr) const noexcept {
        ASSERT(qSqr >= 0.f);
        ASSERT(initialized());
        if (qSqr >= sqr(rad)) {
            // outside of kernel support
            return 0.f;
        }
        const Float floatIdx = qSqrToIdx * qSqr;
        const int idx1 = Size(floatIdx);
        ASSERT(unsigned(idx1) < unsigned(NEntries));
        const int idx2 = idx1 + 1;
        const Float ratio = floatIdx - Float(idx1);
        ASSERT(ratio >= 0.f && ratio < 1.f);

        return grads[idx1] * (1.f - ratio) + (int(idx2 < NEntries) * grads[idx2]) * ratio;
    }
};


/// \brief Cubic spline (M4) kernel
template <typename Float, int Dim>
class CubicSpline : public Kernel<CubicSpline<Float, Dim>, Float, Dim> {
private:
    const Float normalization[3] = { 2.f / 3.f, 10.f / (7.f * PI), 1.f / PI };

public:
    Float radius() const {
        return 2.f;
    }

    // template <bool TApprox = false>
    Float valueImpl(const Float qSqr) const {
        const Float q = sqrt(qSqr);
        ASSERT(q >= 0);
        if (q < 1.f) {
            return normalization[Dim - 1] * (0.25f * pow<3>(2.f - q) - pow<3>(1.f - q));
        }
        if (q < 2.f) {
            return normalization[Dim - 1] * (0.25f * pow<3>(2.f - q));
        }
        return 0.f; // compact within 2r radius
    }

    // template <bool TApprox = false>
    Float gradImpl(const Float qSqr) const {
        const Float q = sqrt(qSqr);
        if (q == 0.f) {
            // gradient of kernel is 0 at q = 0, but here we divide by q,
            // the expression grad/q has a finite limit for q->0
            return -3.f * normalization[Dim - 1];
        }
        if (q < 1.f) {
            return (1.f / q) * normalization[Dim - 1] * (-0.75f * pow<2>(2.f - q) + 3.f * pow<2>(1.f - q));
        }
        if (q < 2.f) {
            return (1.f / q) * normalization[Dim - 1] * (-0.75f * pow<2>(2.f - q));
        }
        return 0.f;
    }
};
/*
/// \brief Gaussian kernel
///
/// Clamped to zero at radius 5, the error is therefore about exp(-5^2) = 10^-11.
template <Size D>
class Gaussian : public Kernel<Gaussian<D>, D> {
private:
    const Float normalization[3] = { 1.f / sqrt(PI), 1.f / PI, 1.f / (PI * sqrt(PI)) };

public:
    Float radius() const {
        return 5.f;
    }

    Float valueImpl(const Float qSqr) const {
        if (qSqr >= sqr(radius())) {
            return 0.f;
        }
        return normalization[D - 1] * exp(-qSqr);
    }

    Float gradImpl(const Float qSqr) const {
        if (qSqr >= sqr(radius())) {
            return 0.f;
        }
        if (qSqr == 0.f) {
            return -2.f * normalization[D - 1];
        }
        const Float q = sqrt(qSqr);
        return normalization[D - 1] / q * exp(-qSqr) * (-2.f * q);
    }
};

/// \brief Triangular (piecewise linear) kernel.
///
/// Does not have continuous derivatives, mainly for testing purposes and non-SPH applications.
template <Size D>
class TriangleKernel : public Kernel<TriangleKernel<D>, D> {
private:
    const Float normalization[3] = { 1.f, 3.f / PI, 3.f / PI };

public:
    Float radius() const {
        return 1.f;
    }

    Float valueImpl(const Float qSqr) const {
        if (qSqr >= sqr(radius())) {
            return 0.f;
        }
        const Float q = sqrt(qSqr);
        return normalization[D - 1] * (1.f - q);
    }

    Float gradImpl(const Float qSqr) const {
        if (qSqr >= sqr(radius())) {
            return 0.f;
        }
        // unfortunately this gradient is nonzero at q->0, so grad/q diverges;
        // let's return a reasonable value to avoid numerical problems
        if (qSqr == 0.f) {
            return -100.f;
        }
        const Float q = sqrt(qSqr);
        return -normalization[D - 1] / q;
    }
};
*/

} // namespace Pvl
