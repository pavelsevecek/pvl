#pragma once

#include "KdTree.hpp"
#include "Matrix.hpp"
#include "PlyWriter.hpp"
#include "Svd.hpp"
#include "utility"
#include <iostream>

namespace Pvl {

template <typename Cloud>
Box3f boundingBox(const Cloud& cloud) {
    Box3f box;
    for (const auto& p : cloud) {
        box.extend(p);
    }
    return box;
}

template <typename Cloud>
Vec3f centroid(const Cloud& cloud) {
    Vec3f pos(0.);
    std::size_t count = 0;
    for (const auto& p : cloud) {
        pos += p;
        count++;
    }
    return pos / count;
}

template <typename ConcurrencyTag = SequentialTag, typename Cloud>
std::vector<Vec3f> estimateNormals(Cloud& cloud) {
    return estimateNormals<ConcurrencyTag>(cloud, [](float) { return false; });
}

template <typename ConcurrencyTag = SequentialTag, typename Cloud, typename Progress>
std::vector<Vec3f> estimateNormals(Cloud& cloud, const Progress& progress) {
    KdTree<Vec3f> tree;
    tree.build(cloud);
    std::vector<Vec3f> normals(cloud.size());
    auto meter = makeProgressMeter(cloud.size(), progress);
    std::atomic<bool> wasCancelled{ false };
    ParallelFor<ConcurrencyTag>()(std::size_t(0), cloud.size(), [&](std::size_t i) {
        if (wasCancelled) {
            return;
        }
        const Vec3f& p = cloud[i];
        float radius = 0.01;
        std::vector<int> neighs;
        do {
            neighs.clear();
            tree.rangeQuery(p, radius, std::back_inserter(neighs));
            radius *= 2.f;
        } while (neighs.size() < 20);

        Vec3f centroid(0.);
        for (int j : neighs) {
            centroid += cloud[j];
        }
        centroid /= neighs.size();
        Mat33f cov = Mat33f::null();
        for (int j : neighs) {
            Vec3f diff = cloud[j] - centroid;
            cov += outerProd(diff, diff);
        }
        Svd<float> svd = singularValueDecomposition(cov);

        normals[i] = svd.U.column(argMin(svd.S));

        // initially orient upwards
        normals[i] *= sign(normals[i][2]);
        if (meter.inc()) {
            wasCancelled = true;
            return;
        }
    });
    if (!wasCancelled) {
        return normals;
    } else {
        return {};
    }
}

Vector<float, 6> join(Vec3f v, Vec3f n) {
    Vector<float, 6> w;
    w[0] = v[0];
    w[1] = v[1];
    w[2] = v[2];
    w[3] = n[0];
    w[4] = n[1];
    w[5] = n[2];
    return w;
}

template <typename Cloud, typename Normals>
void orientNormals(Cloud& cloud, Normals& normals) {
    // get consistently oriented patches

    const std::size_t unassigned = std::size_t(-1);
    std::vector<std::size_t> indices(cloud.size(), unassigned);
    std::size_t componentIdx = 0;

    std::vector<std::size_t> stack;
    std::vector<std::size_t> neighs;

    KdTree<Vec3f> tree;
    tree.build(cloud);

    for (std::size_t i = 0; i < cloud.size(); ++i) {
        if (indices[i] == unassigned) {
            indices[i] = componentIdx;
            stack.push_back(i);
            // find new neigbours recursively until we find all particles in the component
            while (!stack.empty()) {
                const std::size_t n1 = stack.back();
                stack.pop_back();

                const Vec3f& p = cloud[n1];
                float radius = 0.001;
                std::vector<int> neighs;
                do {
                    neighs.clear();
                    tree.rangeQuery(p, radius, std::back_inserter(neighs));
                    radius *= 1.5f;
                } while (neighs.size() < 5);

                for (std::size_t n2 : neighs) {
                    if (n1 == n2) {
                        continue;
                    }
                    Vec3f e = normalize(cloud[n2] - cloud[n1]);
                    Vec3f ndash = normals[n2] - 2 * e * dotProd(e, normals[n2]);
                    if (dotProd(normals[n1], ndash) < 0.8) {
                        // do not count as neighbours
                        continue;
                    }
                    if (indices[n2] == unassigned) {
                        indices[n2] = componentIdx;
                        stack.push_back(n2);
                    }
                }
            }
            componentIdx++;
        }
    }

    // declare component 0 as correct
    for (std::size_t c = 1; c < componentIdx; ++c) {
        int votes = 0;
        for (std::size_t i = 0; i < cloud.size(); ++i) {
            if (indices[i] != c) {
                continue;
            }
            const Vec3f& p = cloud[i];
            float radius = 0.001;
            std::vector<int> neighs;
            do {
                neighs.clear();
                tree.rangeQuery(p, radius, std::back_inserter(neighs));
                radius *= 1.5f;
            } while (neighs.size() < 20);
            for (int j : neighs) {
                if (indices[j] == 0) {
                    Vec3f e = normalize(cloud[j] - cloud[i]);
                    Vec3f ndash = normals[j] - 2 * e * dotProd(e, normals[j]);
                    votes += sign(dotProd(normals[i], ndash));
                }
            }
        }
        if (votes >= 0) {
            continue;
        }
        for (std::size_t i = 0; i < cloud.size(); ++i) {
            if (indices[i] != c) {
                continue;
            }
            normals[i] *= -1;
        }
    }

#if 0
    KdTree<Vec3f> tree;
    tree.build(cloud);

    std::set<std::size_t> visited;
    std::queue<std::size_t> queue;
    queue.push(0);
    while (!queue.empty()) {
        std::cout << "queue size: " << queue.size() << std::endl;
        std::cout << "visited size: " << visited.size() << std::endl;
        std::size_t i = queue.front();
        queue.pop();

        const Vec3f& p = cloud[i];
        float radius = 0.003;
        std::vector<int> neighs;
        do {
            neighs.clear();
            tree.rangeQuery(p, radius, std::back_inserter(neighs));
            radius *= 2.f;
        } while (neighs.size() < 20);

        int votes = 0;
        for (int j : neighs) {
            votes += sign(dot(normals[i], normals[j]));
            if (visited.find(j) == visited.end()) {
                queue.push(j);
                visited.insert(j);
            }
        }
        normals[i] *= sign(votes);
    }
#endif

#if 0
    Vec3f center = centroid(cloud);
    std::vector<Vec3f> cameras{
        center + Vec3f(-100., 0, 0),
        center + Vec3f(100., 0, 0),
        center + Vec3f(0., -100, 0),
        center + Vec3f(0., 100, 0),
        center + Vec3f(0., 0., -100),
        center + Vec3f(0., 0., 100),
    };
    for (std::size_t i = 0; i < cloud.size(); ++i) {
        const Vec3f camera =
            *std::min_element(cameras.begin(), cameras.end(), [&](const Vec3f& cam1, const Vec3f& cam2) {
                return normSqr(cam1 - cloud[i]) < normSqr(cam2 - cloud[i]);
            });
        normals[i] *= sign(dot(normals[i], camera - cloud[i]));
    }
#endif
}

} // namespace Pvl
