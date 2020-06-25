#pragma once

#include "Matrix.hpp"
#include "Optional.hpp"
#include "TriangleMesh.hpp"

namespace Pvl {

static constexpr double ALPHA = 1. * M_PI / 180.; // 1 degree
static double SQR_COS_ALPHA = sqr(std::cos(ALPHA));
static double SQR_SIN_ALPHA = sqr(std::sin(ALPHA));

static constexpr double MAX_PHI = 90. * M_PI / 180.;
static double MIN_COS_PHI = std::cos(MAX_PHI);

template <typename Float>
class LindstromTurkConstraints {
    using Vec = Vector<Float, 3>;
    using Mat = Matrix<Float, 3, 3>;

    Mat A_;
    Vec b_;
    int n_ = 0;

public:
    bool addConstraint(const Vec& p, const Float b) {
        if (full()) {
            return false;
        }

        if (isCompatible(p)) {
            const Float distSqr = dotProd(p, p);
            const Float dist = sqrt(distSqr);
            const Vec pn = p / dist;
            const Float bn = b / dist;

            A_.row(n_) = pn;
            b_[n_] = bn;
            n_++;
            return true;
        } else {
            // std::cout << "Constraint not compatible" << std::endl;
            return false;
        }
    }

    /// Used for problems formulated as quadratic minimization,
    /// i.e. finds point where gradient Hv + c = 0.
    int addQuadraticConstraint(const Mat& H, const Vec& c) {
        int added = 0;
        switch (n_) {
        case 0:
            added += int(addConstraint(H.row(0), -c[0]));
            added += int(addConstraint(H.row(1), -c[1]));
            added += int(addConstraint(H.row(2), -c[2]));
            break;
        case 1: {
            const Vec r0 = A_.row(0);
            std::array<Float, 3> abs_r0 = {
                std::abs(r0[0]), std::abs(r0[1]), std::abs(r0[2])
            };
            Vec q0;
            const int maxIdx =
                int(std::max_element(abs_r0.begin(), abs_r0.end()) - abs_r0.begin());
            switch (maxIdx) {
            case 0:
                q0 = { -r0[2] / r0[0], 0., 1. };
                break;
            case 1:
                q0 = { 0., -r0[2] / r0[1], 1. };
                break;
            case 2:
                q0 = { 1., 0., -r0[0] / r0[2] };
                break;
            }
            const Vec q1 = crossProd(r0, q0);
            const Vec p0 = prod(H, q0);
            const Vec p1 = prod(H, q1);
            const Float b0 = -dotProd(q0, c);
            const Float b1 = -dotProd(q1, c);

            added += int(addConstraint(p0, b0));
            added += int(addConstraint(p1, b1));
            break;
        }
        case 2: {
            const Vec r0 = A_.row(0);
            const Vec r1 = A_.row(1);
            const Vec n = crossProd(r0, r1);
            const Vec p = prod(H, n);
            const Float b = -dotProd(n, c);
            added += int(addConstraint(p, b));
            break;
        }
        default:
            std::cout << "Overconstrained linear system, skipping" << std::endl;
            break;
        }
        return added;
    }

    Vec getPlacement() const {
        const Mat Ainv = invert(A_);
        return prod(Ainv, b_);
    }

    bool full() const {
        return n_ == 3;
    }

    int count() const {
        return n_;
    }

private:
    /// Contraint compatibility rules (see Sec. 4.2 in LT paper)
    bool isCompatible(const Vec& p) const {
        const Float distSqr = normSqr(p);
        if (distSqr < 1.e-20) {
            return false;
        }
        switch (n_) {
        case 0:
            return true;
        case 1: {
            const Vec r0 = A_.row(0);
            const Float pr0 = dotProd(r0, p);
            const Float r0sqr = dotProd(r0, r0);
            const Float lhs = sqr(pr0);
            const Float rhs = r0sqr * distSqr * SQR_COS_ALPHA;
            return lhs < rhs;
        }
        case 2: {
            const Vec r0 = A_.row(0);
            const Vec r1 = A_.row(1);
            const Vec n = crossProd(r0, r1);
            const Float np = dotProd(n, p);
            const Float nsqr = dotProd(n, n);
            const Float lhs = sqr(np);
            const Float rhs = nsqr * distSqr * SQR_SIN_ALPHA;
            return lhs > rhs;
        }
        default:
            throw;
        }
    }
};

template <typename Mesh>
class MemorylessDecimator {
    using Float = typename Mesh::Float;
    using Vec = Vector<Float, 3>;
    using Mat = Matrix<Float, 3, 3>;
    using Constraints = LindstromTurkConstraints<Float>;

public:
    Optional<float> cost(const Mesh& mesh, const Graph::CollapseContext& context) const {
        const Vec p0 = mesh.point(context.remaining);
        const Vec p1 = mesh.point(context.removed);
        const Float lengthSqr = normSqr(p0 - p1);
        std::set<FaceHandle> triangles = getTriangles(mesh, context);
        std::set<VertexHandle> ring = getVertexRing(mesh, context);
        std::set<HalfEdgeHandle> boundary; //= getBoundary(mesh, context);
        Optional<Vec> dp = computePlacement(mesh, p0, triangles, ring, boundary, lengthSqr);
        if (dp) {
            Vec p = p0 + dp.value();

            const Float volumeCost = computeVolumeCost(mesh, triangles, p);
            const Float boundaryCost = computeBoundaryCost(mesh, boundary, p);

            const Float totalCost = 0.5 * boundaryCost * lengthSqr + 0.5 * volumeCost;
            PVL_ASSERT(totalCost >= 0.);

            return totalCost;
        } else {
            return NONE;
        }
    }

    Optional<Vec> placement(const Mesh& mesh, const Graph::CollapseContext& context) const {
        const Vec p0 = mesh.point(context.remaining);
        const Vec p1 = mesh.point(context.removed);
        const Float lengthSqr = normSqr(p0 - p1);
        std::set<FaceHandle> faces = getTriangles(mesh, context);
        std::set<VertexHandle> ring = getVertexRing(mesh, context);
        std::set<HalfEdgeHandle> boundary; //= getBoundary(ci);
        Optional<Vec> dp = computePlacement(mesh, p0, faces, ring, boundary, lengthSqr);
        if (dp) {
            return p0 + dp.value();
        } else {
            return NONE;
        }
    }

    void postprocess(const Mesh&, const Graph::CollapseContext&) {}

private:
    inline Mat crossProductMatrixSqr(const Vec& p) const {
        return {
            Vec(p[1] * p[1] + p[2] * p[2], -p[0] * p[1], -p[0] * p[2]),
            Vec(-p[0] * p[1], p[0] * p[0] + p[2] * p[2], -p[1] * p[2]),
            Vec(-p[0] * p[2], -p[1] * p[2], p[0] * p[0] + p[1] * p[1]),
        };
    }

    Vec computePlacement(const Mesh& mesh,
        const Vec& pivot,
        const std::set<FaceHandle>& triangles,
        const std::set<VertexHandle>& ring,
        const std::set<HalfEdgeHandle>& boundary,
        const Float lengthSqr) const {
        Constraints constraints;

        addBoundaryConstraint(mesh, constraints, boundary);
        addVolumeConstraint(mesh, pivot, constraints, triangles);
        addVolumeAndBoundaryOptimizationConstraint(
            mesh, pivot, constraints, triangles, boundary, lengthSqr);
        addShapeConstraint(mesh, pivot, constraints, ring);

        return constraints.getPlacement();
    }

    void addVolumeConstraint(const Mesh& mesh,
        const Vec& pivot,
        Constraints& constraints,
        const std::set<FaceHandle>& faces) const {
        if (constraints.full()) {
            return;
        }
        Vec sumN(0);
        Float sumL = 0.;

        for (FaceHandle fh : faces) {
            /// \todo add pivot to triangle(fh) ?
            std::array<Vec, 3> tri = mesh.triangle(fh);
            const Vec n = crossProd(tri[1] - tri[0], tri[2] - tri[0]);
            const Float l = dotProd(crossProd(tri[0] - pivot, tri[1] - pivot), tri[2] - pivot);
            sumN += n;
            sumL += l;
        }

        constraints.addConstraint(sumN, sumL);
    }

    void addBoundaryConstraint(const Mesh& mesh,
        Constraints& constraints,
        const std::set<HalfEdgeHandle>& edges) const {
        if (edges.empty() || constraints.full()) {
            return;
        }

        Vec e1(0);
        Vec e2(0);

        for (HalfEdgeHandle eh : edges) {
            const Vec from = mesh.point(mesh.from(eh));
            const Vec to = mesh.point(mesh.to(eh));
            e1 += from - to;
            e2 += crossProd(from, to);
        }

        const Mat H = crossProductMatrixSqr(e1);
        const Vec c = crossProd(e1, e2);

        constraints.addQuadraticConstraint(H, c);
    }

    void addShapeConstraint(const Mesh& mesh,
        const Vec& pivot,
        Constraints& constraints,
        const std::set<VertexHandle>& ring) const {
        if (constraints.full()) {
            return;
        }
        Mat H = Mat::identity() * Float(ring.size());
        Vec c(0);
        for (VertexHandle v : ring) {
            c -= (mesh.point(v) - pivot);
        }
        constraints.addQuadraticConstraint(H, c);
    }

    void addVolumeAndBoundaryOptimizationConstraint(const Mesh& mesh,
        const Vec& pivot,
        Constraints& constraints,
        const std::set<FaceHandle>& faces,
        const std::set<HalfEdgeHandle>& boundary,
        const Float lengthSqr) const {
        if (constraints.full()) {
            std::cout << "Skipping optimization contraint" << std::endl;
            return;
        }

        Mat H = Mat::null();
        Vec c(0, 0, 0);

        for (FaceHandle face : faces) {
            std::array<Vec, 3> tri = mesh.triangle(face);
            const Vec n = crossProd(tri[1] - tri[0], tri[2] - tri[0]);
            const Float l = dotProd(crossProd(tri[0] - pivot, tri[1] - pivot), tri[2] - pivot);

            H += outerProd(n, n);
            c -= n * l;
        }

        if (!boundary.empty()) {
            Mat Hb = Mat::null();
            Vec cb(0, 0, 0);
            for (HalfEdgeHandle edge : boundary) {
                const Vec from = mesh.point(mesh.from(edge));
                const Vec to = mesh.point(mesh.to(edge));
                const Vec n = crossProd(from, to);
                const Vec dir = from - to;

                Hb += crossProductMatrixSqr(dir);
                cb += crossProd(dir, n);
            }
            // 9 * boundary weight * homogenizing factor
            const Float volumeWeight = 0.5;
            const Float boundaryWeight = 0.5;
            const Float w = 9 * boundaryWeight * lengthSqr;

            H = H * volumeWeight + Hb * w;
            c = c * volumeWeight + cb * w;
        }
        constraints.addQuadraticConstraint(H, c);
    }

    std::set<FaceHandle> getTriangles(const Mesh& mesh,
        const Graph::CollapseContext& context) const {
        std::set<FaceHandle> faces;
        for (FaceHandle fh : mesh.faceRing(context.remaining)) {
            faces.insert(fh);
        }
        for (FaceHandle fh : mesh.faceRing(context.removed)) {
            faces.insert(fh);
        }
        faces.erase(context.left);
        faces.erase(context.right);
        return faces;
    }

    /*std::set<HalfedgeHandle> getBoundary(const CollapseInfo& ci) const {
        std::set<HalfedgeHandle> handles;
        for (VertexHandle v : { ci.v0, ci.v1 }) {
            // incomming halfedges
            for (auto iter = mesh_.vih_iter(v); iter.is_valid(); ++iter) {
                HalfedgeHandle eh1 = *iter;
                HalfedgeHandle eh2 = mesh_.opposite_halfedge_handle(eh1);
                if (mesh_.is_boundary(eh1)) {
                    handles.insert(eh1);
                }
                if (mesh_.is_boundary(eh2)) {
                    handles.insert(eh2);
                }
            }
        }
        return handles;
    }*/

    std::set<VertexHandle> getVertexRing(const Mesh& mesh,
        const Graph::CollapseContext& context) const {
        std::set<VertexHandle> ring;
        for (VertexHandle v1 : { context.removed, context.remaining }) {
            // incomming halfedges
            for (VertexHandle v2 : mesh.vertexRing(v1)) {
                if (v2 != context.removed && v2 != context.remaining) {
                    ring.insert(v2);
                }
            }
        }
        return ring;
    }

    Float computeVolumeCost(const Mesh& mesh,
        const std::set<FaceHandle>& faces,
        const Vec& p) const {
        Float cost = 0.;
        for (FaceHandle face : faces) {
            std::array<Vec, 3> tri = mesh.triangle(face);

            const Vec v01 = tri[1] - tri[0];
            const Vec v02 = tri[2] - tri[0];
            const Vec n = crossProd(v01, v02);
            const Float l = dotProd(crossProd(tri[0], tri[1]), tri[2]);
            cost += sqr(dotProd(n, p) - l);
        }
        return cost / 36.;
    }

    Float computeBoundaryCost(const Mesh& mesh,
        const std::set<HalfEdgeHandle>& edges,
        const Vec& p) const {
        Float cost = 0.;
        for (HalfEdgeHandle edge : edges) {
            Vec from = mesh.point(mesh.from(edge));
            Vec to = mesh.point(mesh.to(edge));

            const Vec dir = from - to;
            const Vec c = crossProd(dir, from - p);
            cost += dotProd(c, c);
        }
        return cost / 4.;
    }

    /*  bool isCollapseAllowed(const CollapseInfo& ci, const Point& p) const {
          // simulate collapse
          mesh_.set_point(ci.v0, p);
          mesh_.set_point(ci.v1, p);
          OnScopeExit guard([&] {
              // rollback the collapse
              mesh_.set_point(ci.v0, ci.p0);
              mesh_.set_point(ci.v1, ci.p1);
          });

          for (VertexHandle vh : { ci.v0, ci.v1 }) {
              for (auto iter = mesh_.vf_iter(vh); iter.is_valid(); ++iter) {
                  FaceHandle face = *iter;
                  // these faces will be removed, so no need to test them
                  if (face != ci.fl && face != ci.fr) {
                      typename MeshT::Normal n1 = mesh_.normal(face);
                      typename MeshT::Normal n2 = mesh_.calc_face_normal(face);
                      const Float cosPhi = dot(n1, n2);
                      if (cosPhi < MIN_COS_PHI) {
                          return false;
                      }
                  }
              }
          }
          return true;
      }*/
};

} // namespace Pvl
