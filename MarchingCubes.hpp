#pragma once

#include "TriangleMesh.hpp"
#include <tbb/enumerable_thread_specific.h>

namespace Pvl {

template <typename Float>
class Cell {
private:
    std::array<Vector<Float, 3>, 8> points_;
    std::array<Float, 8> values_;

public:
    Float& value(const int idx) {
        return values_[idx];
    }

    Vector<Float, 3>& node(const int idx) {
        return points_[idx];
    }
};


/*template <typename Volume, typename Float>
void extractIsosurface(const Volume& volume, TriangleMesh& mesh, const float isovalue) {
    typename Volume::ConstIterator iter;
    typename Volume::ConstIterator end(volume.end());
    for (iter = volume.begin(); iter != end; ++iter) {
        Cell cell;
        bool allBelow = true, allAbove = true;
        int i = 0;
        for (int z = 0; z <= 1; ++z) {
            for (int y = 0; y <= 1; ++y) {
                for (int x = 0; x <= 1; ++x) {
                    // the mapping convention is (0, 0) - (1, 0) - (1, 1) - (0, 1)
                    int actX = (x + y) % 2;
                    cell.node(i) = v + dr * Vector(actX, y, z);
                    cell.value(i) = cached.phi[mapping(idxs + Indices(actX, y, z))];
                    if (cell.value(i) > surfaceLevel) {
                        allBelow = false;
                    } else if (cell.value(i) < surfaceLevel) {
                        allAbove = false;
                    }
                    i++;
                }
            }
        }

        if (!allBelow && !allAbove) {
            this->intersectCell(cell, tri.local());
        }
    }

    for (Array<Triangle>& triTl : tri) {
        triangles.pushAll(triTl);
    }
}


template <typename TFunctor>
bool MarchingCubes::iterateWithIndices(const Box& box, const Vector& step, TFunctor&& functor) {
    MEASURE_SCOPE("MC - evaluating field");
    ASSERT(box != Box::EMPTY());

    Size reportCnt = max(Size(box.size()[Z] / step[Z]), 1u);
    Size reportStep = max(reportCnt / 100, 1u);
    std::atomic_int counter{ 0 };
    std::atomic_bool shouldContinue{ true };
    auto task = [this, &step, &box, &functor, reportStep, reportCnt, &counter, &shouldContinue](
                    const Size k) {
        const Float z = box.lower()[Z] + k * step[Z];
        Size i = 0;
        Size j = 0;
        for (Float y = box.lower()[Y]; y <= box.upper()[Y]; y += step[Y], j++) {
            i = 0;
            for (Float x = box.lower()[X]; x <= box.upper()[X]; x += step[X], i++) {
                functor(Indices(i, j, k), Vector(x, y, z));
            }
        }

        if (progressCallback && (++counter % reportStep == 0)) {
            shouldContinue = shouldContinue && progressCallback(Float(counter) / Float(reportCnt));
        }
    };
    parallelFor(scheduler, 0, Size(box.size()[Z] / step[Z]) + 1, task);
    return shouldContinue;
}

MarchingCubes::MarchingCubes(IScheduler& scheduler,
    const Float surfaceLevel,
    const SharedPtr<IScalarField>& field,
    Function<bool(Float progress)> progressCallback)
    : scheduler(scheduler)
    , surfaceLevel(surfaceLevel)
    , field(field)
    , progressCallback(progressCallback) {}

void MarchingCubes::addComponent(const Box& box, const Float gridResolution) {
    MEASURE_SCOPE("MC addComponent");

    const Vector dr = min(Vector(gridResolution), box.size() * (1._f - EPS));
    cached.phi.clear();
    // multiply by (1 + EPS) to handle case where box size is divisible by dr
    Indices cnts((1._f + EPS) * box.size() / dr);
    ASSERT(cnts[X] >= 1 && cnts[Y] >= 1 && cnts[Z] >= 1);

    // find values of grid nodes
    auto mapping = [&cnts](const Indices& idxs) {
        ASSERT(idxs[X] >= 0 && idxs[X] <= cnts[X], idxs[X], cnts[X]);
        ASSERT(idxs[Y] >= 0 && idxs[Y] <= cnts[Y], idxs[Y], cnts[Y]);
        ASSERT(idxs[Z] >= 0 && idxs[Z] <= cnts[Z], idxs[Z], idxs[Z]);
        return idxs[X] + (cnts[X] + 1) * idxs[Y] + (cnts[X] + 1) * (cnts[Y] + 1) * idxs[Z];
    };
    cached.phi.resize((cnts[X] + 1) * (cnts[Y] + 1) * (cnts[Z] + 1));
    bool shouldContinue =
        this->iterateWithIndices(box, dr, [this, &mapping](const Indices& idxs, const Vector& v) { //
            cached.phi[mapping(idxs)] = (*field)(v);
        });
    if (!shouldContinue) {
        return;
    }

    // for each non-empty grid, find all intersecting triangles
    Box boxWithoutLast(box.lower(), box.upper() - dr);
    ThreadLocal<Array<Triangle>> tri(scheduler);

    auto intersect = [this, &dr, &mapping, &tri](const Indices& idxs, const Vector& v) {
        Cell cell;
        bool allBelow = true, allAbove = true;
        Size i = 0;
        for (Size z = 0; z <= 1; ++z) {
            for (Size y = 0; y <= 1; ++y) {
                for (Size x = 0; x <= 1; ++x) {
                    Size actX = (x + y) % 2; // the mapping convention is (0, 0) - (1, 0) - (1, 1) - (0, 1)
                    cell.node(i) = v + dr * Vector(actX, y, z);
                    cell.value(i) = cached.phi[mapping(idxs + Indices(actX, y, z))];
                    if (cell.value(i) > surfaceLevel) {
                        allBelow = false;
                    } else if (cell.value(i) < surfaceLevel) {
                        allAbove = false;
                    }
                    i++;
                }
            }
        }

        if (!allBelow && !allAbove) {
            this->intersectCell(cell, tri.local());
        }
    };
    shouldContinue = this->iterateWithIndices(boxWithoutLast, dr, intersect);
    if (!shouldContinue) {
        return;
    }

    for (Array<Triangle>& triTl : tri) {
        triangles.pushAll(triTl);
    }
}

void MarchingCubes::intersectCell(Cell& cell, Array<Triangle>& tri) {
    Size cubeIdx = 0;
    for (Size i = 0; i < 8; ++i) {
        if (cell.value(i) <= surfaceLevel) {
            cubeIdx |= 1 << i;
        }
    }

    if (MC_EDGES[cubeIdx] == 0) {
        // cube is entirely in/out of the surface
        return;
    }

    // find the vertices where the surface intersects the cube
    StaticArray<Vector, 12> vertices;
    for (Size i = 0; i < 12; ++i) {
        if (MC_EDGES[cubeIdx] & (1 << i)) {
            const Size k = IDXS1[i];
            const Size l = IDXS2[i];
            vertices[i] = this->interpolate(cell.node(k), cell.value(k), cell.node(l), cell.value(l));
        } else {
            vertices[i] = Vector(NAN);
        }
    }

    for (Size i = 0; MC_TRIANGLES[cubeIdx][i] != -1; i += 3) {
        Triangle t;
        t[0] = vertices[MC_TRIANGLES[cubeIdx][i + 0]];
        t[1] = vertices[MC_TRIANGLES[cubeIdx][i + 1]];
        t[2] = vertices[MC_TRIANGLES[cubeIdx][i + 2]];
        if (!t.isValid()) {
            // skip degenerated triangles
            continue;
        }
        tri.push(t);
    }
}

INLINE Vector MarchingCubes::interpolate(const Vector& v1,
    const Float p1,
    const Vector& v2,
    const Float p2) const {
    if (almostEqual(p1, surfaceLevel)) {
        return v1;
    }
    if (almostEqual(p2, surfaceLevel)) {
        return v2;
    }
    if (almostEqual(p1, p2)) {
        // small difference between values, just return the center to avoid instabilities
        return 0.5_f * (v1 + v2);
    }

    Float mu = (p1 - surfaceLevel) / (p1 - p2);
    ASSERT(mu >= 0._f && mu <= 1._f);
    return v1 + mu * (v2 - v1);
}

namespace {

    class ColorField : public IScalarField {
    private:
        LutKernel<3> kernel;
        AutoPtr<IBasicFinder> finder;

        ArrayView<const Vector> r;
        ArrayView<const Float> m, rho;
        ArrayView<const Size> flag;
        ArrayView<const SymmetricTensor> G;
        Float maxH = 0._f;

        ThreadLocal<Array<NeighbourRecord>> neighs;

    public:
        /// \brief Creates the number density field.
        ///
        /// This implementation uses anisotropic kernel to reduce the perturbations of the boundary, see
        /// https://www.cc.gatech.edu/~turk/my_papers/sph_surfaces.pdf.
        /// \param storage Storage containing particle masses, densities and flags used to distinguish
        /// different
        ///                bodies (we don't want to blend together their surfaces)
        /// \param scheduler Scheduler used for parallelization.
        /// \param r Particle positions, generally different than the ones stored in the storage.
        /// \param aniso Particle anisotropy matrix, for isotropic distribution equals to I/h
        /// \param kernel SPH kernel used for particle smoothing
        /// \param finder Neighbour finder
        ColorField(const Storage& storage,
            IScheduler& scheduler,
            const ArrayView<const Vector> r,
            const ArrayView<const SymmetricTensor> aniso,
            const Float maxH,
            LutKernel<3>&& kernel,
            AutoPtr<IBasicFinder>&& finder)
            : kernel(std::move(kernel))
            , finder(std::move(finder))
            , r(r)
            , G(aniso)
            , maxH(maxH)
            , neighs(scheduler) {
            tie(m, rho) = storage.getValues<Float>(QuantityId::MASS, QuantityId::DENSITY);
            flag = storage.getValue<Size>(QuantityId::FLAG);

            // we have to re-build the tree since we are using different positions (in general)
            this->finder->build(scheduler, r);
        }

        virtual Float operator()(const Vector& pos) override {
            ASSERT(maxH > 0._f);
            Array<NeighbourRecord>& neighsTl = neighs.local();
            /// \todo for now let's just search some random multiple of smoothing length, we should use the
            /// largest singular value here
            finder->findAll(pos, maxH * kernel.radius(), neighsTl);
            Float phi = 0._f;

            // find average h of neighbours and the flag of the closest particle
            Size closestFlag = 0;
            Float flagDistSqr = INFTY;
            for (NeighbourRecord& n : neighsTl) {
                const Size j = n.index;
                if (n.distanceSqr < flagDistSqr) {
                    closestFlag = flag[j];
                    flagDistSqr = n.distanceSqr;
                }
            }

            // interpolate values of neighbours
            for (NeighbourRecord& n : neighsTl) {
                const Size j = n.index;
                if (flag[j] != closestFlag) {
                    continue;
                }
                phi +=
                    m[j] / rho[j] * G[j].determinant() * kernel.valueImpl(getSqrLength(G[j] * (pos - r[j])));
            }
            return phi;
        }
    };


    class FallbackField : public IScalarField {
    private:
        LutKernel<3> kernel;
        AutoPtr<IBasicFinder> finder;

        ArrayView<const Vector> r;
        ArrayView<const SymmetricTensor> G;
        Float maxH = 0._f;

        ThreadLocal<Array<NeighbourRecord>> neighs;

    public:
        FallbackField(IScheduler& scheduler,
            const ArrayView<const Vector> r,
            const ArrayView<const SymmetricTensor> aniso,
            const Float maxH,
            LutKernel<3>&& kernel,
            AutoPtr<IBasicFinder>&& finder)
            : kernel(std::move(kernel))
            , finder(std::move(finder))
            , r(r)
            , G(aniso)
            , maxH(maxH)
            , neighs(scheduler) {

            // we have to re-build the tree since we are using different positions (in general)
            this->finder->build(scheduler, r);
        }

        virtual Float operator()(const Vector& pos) override {
            ASSERT(maxH > 0._f);
            Array<NeighbourRecord>& neighsTl = neighs.local();
            /// \todo for now let's just search some random multiple of smoothing length, we should use the
            /// largest singular value here
            finder->findAll(pos, maxH * kernel.radius(), neighsTl);
            Float phi = 0._f;

            // interpolate values of neighbours
            for (NeighbourRecord& n : neighsTl) {
                const Size j = n.index;
                phi += sphereVolume(0.5_f * r[j][H]) * G[j].determinant() *
                       kernel.valueImpl(getSqrLength(G[j] * (pos - r[j])));
            }
            return phi;
        }
    };

} // namespace

INLINE Float weight(const Vector& r1, const Vector& r2) {
    const Float lengthSqr = getSqrLength(r1 - r2);
    // Eq. (11)
    if (lengthSqr < sqr(2._f * r1[H])) {
        return 1._f - pow<3>(sqrt(lengthSqr) / (2._f * r1[H]));
    } else {
        return 0._f;
    }
}

Array<Triangle> getSurfaceMesh(IScheduler& scheduler, const Storage& storage, const McConfig& config) {
    MEASURE_SCOPE("getSurfaceMesh");

    // (according to http://www.cc.gatech.edu/~turk/my_papers/sph_surfaces.pdf)

    ArrayView<const Vector> r = storage.getValue<Vector>(QuantityId::POSITION);
    RunSettings settings;
    LutKernel<3> kernel = Factory::getKernel<3>(settings);
    AutoPtr<IBasicFinder> finder = Factory::getFinder(settings);

    finder->build(scheduler, r);

    Array<Vector> r_bar(r.size());
    Array<SymmetricTensor> G(r.size()); // anisotropy matrix

    ThreadLocal<Array<NeighbourRecord>> neighsData(scheduler);
    parallelFor(scheduler, neighsData, 0, r.size(), [&](const Size i, Array<NeighbourRecord>& neighs) {
        /// \todo point cloud denoising?
        r_bar[i] = r[i];
        r_bar[i][H] = r[i][H] * config.smoothingMult;

        if (config.useAnisotropicKernels) {
            Vector r_center = Vector(0._f);
            finder->findAll(r_bar[i], 2 * r_bar[i][H], neighs);
            for (const NeighbourRecord& n : neighs) {
                r_center += r_bar[n.index];
            }
            r_center /= neighs.size();
            SymmetricTensor C = SymmetricTensor::null();
            for (const NeighbourRecord& n : neighs) {
                C += symmetricOuter(r[n.index] - r_center, r[n.index] - r_center);
            }

            Svd svd = singularValueDecomposition(C);
            const Float maxSigma = maxElement(svd.S);
            for (Size i = 0; i < 3; ++i) {
                svd.S[i] = 1._f / std::max(svd.S[i], 0.125_f * maxSigma);
            }

            AffineMatrix sigma = convert<AffineMatrix>(SymmetricTensor(svd.S, Vector(0._f)));
            G[i] = convert<SymmetricTensor>(svd.V * sigma * svd.U.transpose());
        } else {
            G[i] = SymmetricTensor(Vector(1._f / r[i][H]), Vector(0._f));
        }
    });
    // 5. find bounding box and maximum h (we need to search neighbours of arbitrary point in space)

    Float maxH = 0._f;
    for (Size i = 0; i < r_bar.size(); ++i) {
        maxH = max(maxH, r_bar[i][H]);
    }
    SharedPtr<IScalarField> field;

    if (storage.has(QuantityId::MASS) && storage.has(QuantityId::DENSITY) && storage.has(QuantityId::FLAG)) {
        field =
            makeShared<ColorField>(storage, scheduler, r_bar, G, maxH, std::move(kernel), std::move(finder));
    } else {
        field = makeShared<FallbackField>(scheduler, r_bar, G, maxH, std::move(kernel), std::move(finder));
    }

    MarchingCubes mc(scheduler, config.surfaceLevel, field, config.progressCallback);

    Array<Size> components;
    const Size numComponents = Post::findComponents(storage, 2._f, Post::ComponentFlag::OVERLAP, components);

    // 6. find the surface using marching cubes for each component
    Array<Box> boxes(numComponents);
    Array<Size> counts(numComponents);
    counts.fill(0);
    for (Size j = 0; j < components.size(); ++j) {
        const Vector padding(max(2._f * r_bar[j][H], 2._f * config.gridResolution));
        boxes[components[j]].extend(r_bar[j] + padding);
        boxes[components[j]].extend(r_bar[j] - padding);
        counts[components[j]]++;
    }
    for (Size i = 0; i < numComponents; ++i) {
        if (counts[i] > 10) {
            mc.addComponent(boxes[i], config.gridResolution);
        }
    }

    return std::move(mc.getTriangles());
}

*/
} // namespace Pvl
