#pragma once

#include "Matrix.hpp"
#include "TriangleMesh.hpp"

namespace Pvl {

static constexpr double ALPHA = 1. * M_PI / 180.; // 1 degree
static double SQR_COS_ALPHA = sqr(std::cos(ALPHA));
static double SQR_SIN_ALPHA = sqr(std::sin(ALPHA));

static constexpr double MAX_PHI = 90. * M_PI / 180.;
static double MIN_COS_PHI = std::cos(MAX_PHI);

/*class OnScopeExit : boost::noncopyable {
    std::function<void()> closure_;

public:
    OnScopeExit(std::function<void()>&& closure)
        : closure_(std::move(closure)) {}

    ~OnScopeExit() {
        closure_();
    }
};*/

class LindstromTurkConstraints {
    Mat33f A_;
    Vec3f b_;
    int n_ = 0;

public:
    bool addConstraint(const Vec3f& p, const double b) {
        if (full()) {
            return false;
        }

        if (isCompatible(p)) {
            const double distSqr = dot(p, p);
            const double dist = sqrt(distSqr);
            const Vec3f pn = p / dist;
            const double bn = b / dist;

            A_.row(n_) = pn;
            b_[n_] = bn;
            n_++;
            return true;
        } else {
            std::cout << "Constraint not compatible" << std::endl;
            return false;
        }
    }

    /// Used for problems formulated as quadratic minimization,
    /// i.e. finds point where gradient Hv + c = 0.
    int addQuadraticConstraint(const Mat33f& H, const Vec3f& c) {
        int added = 0;
        switch (n_) {
        case 0:
            added += int(addConstraint(H.row(0), -c[0]));
            added += int(addConstraint(H.row(1), -c[1]));
            added += int(addConstraint(H.row(2), -c[2]));
            break;
        case 1: {
            const Vec3f r0 = A_.row(0);
            std::array<double, 3> abs_r0 = {
                std::abs(r0[0]), std::abs(r0[1]), std::abs(r0[2])
            };
            Vec3f q0;
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
            const Vec3f q1 = cross(r0, q0);
            const Vec3f p0 = H.transform(q0);
            const Vec3f p1 = H.transform(q1);
            const double b0 = -dot(q0, c);
            const double b1 = -dot(q1, c);

            added += int(addConstraint(p0, b0));
            added += int(addConstraint(p1, b1));
            break;
        }
        case 2: {
            const Vec3f r0 = A_.row(0);
            const Vec3f r1 = A_.row(1);
            const Vec3f n = cross(r0, r1);
            const Vec3f p = H.transform(n);
            const double b = -dot(n, c);
            added += int(addConstraint(p, b));
            break;
        }
        default:
            std::cout << "Overconstrained linear system, skipping" << std::endl;
            break;
        }
        return added;
    }

    Vec3f getPlacement() const {
        const Mat33f Ainv = invert(A_);
        return Ainv.transform(b_);
    }

    bool full() const {
        return n_ == 3;
    }

    int count() const {
        return n_;
    }

private:
    /// Contraint compatibility rules (see Sec. 4.2 in LT paper)
    bool isCompatible(const Vec3f& p) const {
        const double distSqr = normSqr(p);
        if (distSqr < 1.e-20) {
            return false;
        }
        switch (n_) {
        case 0:
            return true;
        case 1: {
            const Vec3f r0 = A_.row(0);
            const double pr0 = dot(r0, p);
            const double r0sqr = dot(r0, r0);
            const double lhs = sqr(pr0);
            const double rhs = r0sqr * distSqr * SQR_COS_ALPHA;
            return lhs < rhs;
        }
        case 2: {
            const Vec3f r0 = A_.row(0);
            const Vec3f r1 = A_.row(1);
            const Vec3f n = cross(r0, r1);
            const double np = dot(n, p);
            const double nsqr = dot(n, n);
            const double lhs = sqr(np);
            const double rhs = nsqr * distSqr * SQR_SIN_ALPHA;
            return lhs > rhs;
        }
        default:
            throw;
        }
    }
};

class MemorylessDecimator {
    using Mesh = TriangleMesh<Vec3f, int>;

public:
    float cost(const Mesh& mesh, const Graph::CollapseContext& context) const {
        const Vec3f p0 = mesh.point(context.remaining);
        const Vec3f p1 = mesh.point(context.removed);
        const double lengthSqr = normSqr(p0 - p1);
        std::set<FaceHandle> triangles = getTriangles(mesh, context);
        std::set<VertexHandle> ring = getVertexRing(mesh, context);
        std::set<HalfEdgeHandle> boundary; //= getBoundary(mesh, context);
        Vec3f p = computePlacement(mesh, triangles, ring, boundary, lengthSqr);
        /*if (!p || !isCollapseAllowed(ci, *p)) {
            return Base::ILLEGAL_COLLAPSE;
        }*/

        const double volumeCost = computeVolumeCost(mesh, triangles, p);
        const double boundaryCost = computeBoundaryCost(mesh, boundary, p);

        const double totalCost = 0.5 * boundaryCost * lengthSqr + 0.5 * volumeCost;
        PVL_ASSERT(totalCost >= 0.);

        return totalCost;
    }

    Vec3f placement(const Mesh& mesh, const Graph::CollapseContext& context) const {
        const Vec3f p0 = mesh.point(context.remaining);
        const Vec3f p1 = mesh.point(context.removed);
        const double lengthSqr = normSqr(p0 - p1);
        std::set<FaceHandle> faces = getTriangles(mesh, context);
        std::set<VertexHandle> ring = getVertexRing(mesh, context);
        std::set<HalfEdgeHandle> boundary; //= getBoundary(ci);
        return computePlacement(mesh, faces, ring, boundary, lengthSqr);
    }

    void postprocess(const Mesh&, const Graph::CollapseContext&) {}

private:
    inline Mat33f crossProductMatrixSqr(const Vec3f& p) const {
        return {
            Vec3f(p[1] * p[1] + p[2] * p[2], -p[0] * p[1], -p[0] * p[2]),
            Vec3f(-p[0] * p[1], p[0] * p[0] + p[2] * p[2], -p[1] * p[2]),
            Vec3f(-p[0] * p[2], -p[1] * p[2], p[0] * p[0] + p[1] * p[1]),
        };
    }

    Vec3f computePlacement(const Mesh& mesh,
        const std::set<FaceHandle>& triangles,
        const std::set<VertexHandle>& ring,
        const std::set<HalfEdgeHandle>& boundary,
        const double lengthSqr) const {
        LindstromTurkConstraints constraints;

        addVolumeAndBoundaryOptimizationConstraint(
            mesh, constraints, triangles, boundary, lengthSqr);
        std::cout << "after optimization - " << constraints.count() << " constraints"
                  << std::endl;
        addBoundaryConstraint(mesh, constraints, boundary);
        addVolumeConstraint(mesh, constraints, triangles);
        std::cout << "after volume - " << constraints.count() << " constraints" << std::endl;
        addShapeConstraint(mesh, constraints, ring);
        std::cout << "after shape - " << constraints.count() << " constraints" << std::endl;

        Vec3f p = constraints.getPlacement();
        return p;
    }

    void addVolumeConstraint(const Mesh& mesh,
        LindstromTurkConstraints& constraints,
        const std::set<FaceHandle>& faces) const {
        if (constraints.full()) {
            return;
        }
        Vec3f sumN(0);
        double sumL = 0.;

        for (FaceHandle fh : faces) {
            std::array<Vec3f, 3> tri = mesh.triangle(fh);
            const Vec3f n = cross(tri[1] - tri[0], tri[2] - tri[0]);
            const double l = dot(cross(tri[0], tri[1]), tri[2]);
            sumN += n;
            sumL += l;
        }

        constraints.addConstraint(sumN, sumL);
    }

    void addBoundaryConstraint(const Mesh& mesh,
        LindstromTurkConstraints& constraints,
        const std::set<HalfEdgeHandle>& edges) const {
        if (edges.empty() || constraints.full()) {
            return;
        }

        Vec3f e1(0);
        Vec3f e2(0);

        for (HalfEdgeHandle eh : edges) {
            const Vec3f from = mesh.point(mesh.from(eh));
            const Vec3f to = mesh.point(mesh.to(eh));
            e1 += from - to;
            e2 += cross(from, to);
        }

        const Mat33f H = crossProductMatrixSqr(e1);
        const Vec3f c = cross(e1, e2);

        constraints.addQuadraticConstraint(H, c);
    }

    void addShapeConstraint(const Mesh& mesh,
        LindstromTurkConstraints& constraints,
        const std::set<VertexHandle>& ring) const {
        if (constraints.full()) {
            return;
        }
        Mat33f H = Mat33f::identity() * double(ring.size());
        Vec3f c(0);
        for (VertexHandle v : ring) {
            c -= mesh.point(v);
        }
        constraints.addQuadraticConstraint(H, c);
    }

    void addVolumeAndBoundaryOptimizationConstraint(const Mesh& mesh,
        LindstromTurkConstraints& constraints,
        const std::set<FaceHandle>& faces,
        const std::set<HalfEdgeHandle>& boundary,
        const double lengthSqr) const {
        if (constraints.full()) {
            std::cout << "Skipping optimization contraint" << std::endl;
            return;
        }

        Mat33f H = Mat33f::null();
        Vec3f c(0, 0, 0);

        for (FaceHandle face : faces) {
            std::array<Vec3f, 3> tri = mesh.triangle(face);
            const Vec3f n = cross(tri[1] - tri[0], tri[2] - tri[0]);
            const double l = dot(cross(tri[0], tri[1]), tri[2]);

            H += outerProd(n, n);
            c -= n * l;
        }

        if (!boundary.empty()) {
            Mat33f Hb = Mat33f::null();
            Vec3f cb(0, 0, 0);
            for (HalfEdgeHandle edge : boundary) {
                const Vec3f from = mesh.point(mesh.from(edge));
                const Vec3f to = mesh.point(mesh.to(edge));
                const Vec3f n = cross(from, to);
                const Vec3f dir = from - to;

                Hb += crossProductMatrixSqr(dir);
                cb += cross(dir, n);
            }
            // 9 * boundary weight * homogenizing factor
            const double volumeWeight = 0.5;
            const double boundaryWeight = 0.5;
            const double w = 9 * boundaryWeight * lengthSqr;

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

    double computeVolumeCost(const Mesh& mesh,
        const std::set<FaceHandle>& faces,
        const Vec3f& p) const {
        double cost = 0.;
        for (FaceHandle face : faces) {
            std::array<Vec3f, 3> tri = mesh.triangle(face);

            const Vec3f v01 = tri[1] - tri[0];
            const Vec3f v02 = tri[2] - tri[0];
            const Vec3f n = cross(v01, v02);
            const double l = dot(cross(tri[0], tri[1]), tri[2]);
            cost += sqr(dot(n, p) - l);
        }
        return cost / 36.;
    }

    double computeBoundaryCost(const Mesh& mesh,
        const std::set<HalfEdgeHandle>& edges,
        const Vec3f& p) const {
        double cost = 0.;
        for (HalfEdgeHandle edge : edges) {
            Vec3f from = mesh.point(mesh.from(edge));
            Vec3f to = mesh.point(mesh.to(edge));

            const Vec3f dir = from - to;
            const Vec3f c = cross(dir, from - p);
            cost += dot(c, c);
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
                      const double cosPhi = dot(n1, n2);
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
