#pragma once

#include "Box.hpp"
#include "Utils.hpp"
#include <atomic>
#include <cstdint>
#include <set>
#include <shared_mutex>
//#include <tbb/tbb.h>
#include <type_traits>
#include <vector>


namespace Pvl {

/// \brief Base class for nodes of K-d tree.
///
/// Can be derived to include additional user data for each node.
template <typename Float, int Dim>
struct KdNode {
    /// Here X, Y, Z must be 0, 1, 2
    int8_t type;

    /// Bounding box of particles in the node
    BoundingBox<Vector<Float, Dim>> box;

    struct LeafTag {};

    KdNode(const int8_t type)
        : type(type) {}

    KdNode(const LeafTag)
        : type(Dim) {}

    bool isLeaf() const {
        return type == Dim;
    }
};

/// \brief Inner node of K-d tree
template <typename Float, typename Index, int Dim>
struct KdInnerNode : public KdNode<Float, Dim> {
    /// Position where the selected dimension is split
    Float splitPosition;

    /// Index of left child node
    Index left;

    /// Index of right child node
    Index right;

    KdInnerNode()
        : KdNode<Float, Dim>(-1) {}

    KdInnerNode(const int8_t type)
        : KdNode<Float, Dim>(type) {}
};

/// \brief Leaf (bucket) node of K-d tree
template <typename Float, typename Index, int Dim>
struct KdLeafNode : public KdNode<Float, Dim> {
    /// First index of particlse belonging to the leaf
    Index from;

    /// One-past-last index of particles belonging to the leaf
    Index to;

    /// Unused, used so that LeafNode and InnerNode have the same size
    Index padding;

    KdLeafNode()
        : KdNode<Float, Dim>(typename KdNode<Float, Dim>::LeafTag{}) {}

    void setLeaf() {
        this->type = Dim;
    }

    /// Returns the number of points in the leaf. Can be zero.
    Index size() const {
        return to - from;
    }
};

// static_assert(sizeof(Size) == sizeof(float), "Sizes must match to keep this layout");

/// \brief Index iterator with given mapping (index permutation).
///
/// Returns value mapping[index] when dereferenced,
/*class LeafIndexIterator : public IndexIterator {
private:
    ArrayView<const Size> mapping;

public:
     LeafIndexIterator(const Size idx, ArrayView<const Size> mapping)
        : IndexIterator(idx)
        , mapping(mapping) {}

     Size operator*() const {
        return mapping[idx];
    }
};

/// \brief Helper index sequence to iterate over particle indices of a leaf node.
class LeafIndexSequence : public IndexSequence {
private:
    ArrayView<const Size> mapping;

public:
     LeafIndexSequence(const Size from, const Size to, ArrayView<const Size> mapping)
        : IndexSequence(from, to)
        , mapping(mapping) {
        ASSERT(to <= mapping.size());
    }

     LeafIndexIterator begin() const {
        return LeafIndexIterator(from, mapping);
    }

     LeafIndexIterator end() const {
        return LeafIndexIterator(to, mapping);
    }
};
*/

struct EuclideanMetric {
    template <typename T, int Dim>
    T lengthSqr(const Vector<T, Dim>& v) const {
        return normSqr(v);
    }
};


enum class KdChild {
    LEFT = 0,
    RIGHT = 1,
};

template <typename Index, typename Float>
struct Neighbor {
    Index index;
    Float distSqr;

    Neighbor() = default;
    Neighbor(Index index, Float distSqr)
        : index(index)
        , distSqr(distSqr) {}

    operator Index() const {
        return index;
    }
};

/// \brief K-d tree, used for hierarchical clustering of particles and accelerated Kn queries.
///
/// Allows storing arbitrary data at each node of the tree.
///
/// https://www.cs.umd.edu/~mount/Papers/cgc99-smpack.pdf
/// \tparam TNode Nodes of the tree, should always derive from KdNode and should be POD structs.
/// \tparam TMetric Functor returning the squared distance of two vectors.
template <typename Vec, typename Index = uint32_t, typename Metric = EuclideanMetric>
class KdTree : public Noncopyable {
private:
    struct {
        /// Maximal number of particles in the leaf node
        Index leafSize;

        /// Maximal depth for which the build is parallelized
        Index maxParallelDepth;
    } config_;

    using Float = typename Vec::Float;
    static constexpr int Dim = Vec::size();
    using Box = BoundingBox<Vector<Float, Dim>>;
    using InnerNode = KdInnerNode<Float, Index, Dim>;
    using LeafNode = KdLeafNode<Float, Index, Dim>;
    using IndexDiff = std::make_signed_t<Index>;

    /// Holds all nodes, either \ref InnerNode or \ref LeafNode (depending on the value of \ref type).
    Box entireBox_;
    /// \todo optionally weak copy?
    std::vector<Vec> values_;
    std::vector<Index> idxs_;
    std::vector<InnerNode> nodes_;
    std::atomic_int nodeCounter_;
    // std::shared_timed_mutex nodesMutex_;

    static constexpr Index ROOT_PARENT_NODE = Index(-1);

public:
    explicit KdTree(const int leafSize = 25, const int maxParallelDepth = 50) {
        PVL_ASSERT(leafSize >= 1);
        config_.leafSize = leafSize;
        config_.maxParallelDepth = maxParallelDepth;
    }

    template <typename TContainer>
    void build(const TContainer& points) {
        static_assert(sizeof(LeafNode) == sizeof(InnerNode), "Sizes of nodes must match");

        // clean the current tree
        const Index currentCnt = nodes_.size();
        this->init();

        Index index = 0;
        for (const Vec& p : points) {
            entireBox_.extend(p);
            values_.push_back(p);
            idxs_.push_back(index);
            index++;
        }

        if (points.empty()) {
            return;
        }

        const Index nodeCnt = std::max<Index>(2 * points.size() / config_.leafSize + 1, currentCnt);
        nodes_.resize(nodeCnt);

        // tbb::task_group tg;
        buildTree(ROOT_PARENT_NODE, KdChild(-1), 0, points.size(), entireBox_, 0, 0);
        // tg.wait();

        // shrink nodes to only the constructed ones
        // nodes.resize(nodeCounter);

        // ASSERT(this->sanityCheck(), this->sanityCheck().error());
    }

    template <typename TOutIter>
    Index rangeQuery(const Vec& pos, const Float radius, TOutIter neighs) const {
        struct TraversalNode {
            Index idx;
            Vec sizeSqr;
            Float distanceSqr;
        };
        static thread_local std::vector<TraversalNode> nodeStack;

        const Float radiusSqr = sqr(radius);
        const Vec maxDistSqr = sqr(max(Vec(0), max(entireBox_.lower() - pos, pos - entireBox_.upper())));

        // L1 norm
        const Float l1 = normL1(maxDistSqr);
        TraversalNode node{ 0, maxDistSqr, l1 };

        PVL_ASSERT(nodeStack.empty()); // not sure if there can be some nodes from previous search ...

        Index neighCnt = 0;
        Metric metric;
        while (node.distanceSqr < radiusSqr) {
            if (nodes_[node.idx].isLeaf()) {
                // for leaf just add all
                const LeafNode& leaf = (const LeafNode&)nodes_[node.idx];
                if (leaf.size() > 0) {
                    const Float leafDistSqr =
                        metric.lengthSqr(max(Vec(0), max(leaf.box.lower() - pos, pos - leaf.box.upper())));
                    if (leafDistSqr < radiusSqr) {
                        // leaf intersects the sphere
                        for (Index i = leaf.from; i < leaf.to; ++i) {
                            const Index actIndex = idxs_[i];
                            const Float distSqr = metric.lengthSqr(values_[actIndex] - pos);
                            if (distSqr < radiusSqr) {
                                *neighs++ = Neighbor<Index, Float>{ actIndex, distSqr };
                                neighCnt++;
                            }
                        }
                    }
                }
                if (nodeStack.empty()) {
                    break;
                }
                node = nodeStack.back();
                nodeStack.pop_back();
            } else {
                // inner node
                const InnerNode& inner = (InnerNode&)nodes_[node.idx];
                const Index splitDimension = Index(inner.type);
                PVL_ASSERT(splitDimension < Dim);
                const Float splitPosition = inner.splitPosition;
                if (pos[splitDimension] < splitPosition) {
                    // process left subtree, put right on stack
                    TraversalNode right = node;
                    node.idx = inner.left;

                    const Float dx = splitPosition - pos[splitDimension];
                    right.distanceSqr += sqr(dx) - right.sizeSqr[splitDimension];
                    right.sizeSqr[splitDimension] = sqr(dx);
                    if (right.distanceSqr < radiusSqr) {
                        const InnerNode& next = (const InnerNode&)nodes_[right.idx];
                        right.idx = next.right;
                        nodeStack.push_back(right);
                    }
                } else {
                    // process right subtree, put left on stack
                    TraversalNode left = node;
                    node.idx = inner.right;
                    const Float dx = splitPosition - pos[splitDimension];
                    left.distanceSqr += sqr(dx) - left.sizeSqr[splitDimension];
                    left.sizeSqr[splitDimension] = sqr(dx);
                    if (left.distanceSqr < radiusSqr) {
                        const InnerNode& next = (const InnerNode&)nodes_[left.idx];
                        left.idx = next.left;
                        nodeStack.push_back(left);
                    }
                }
            }
        }

        return neighCnt;
    }

    /// \brief Returns the node with given index
    /*TNode& getNode(const Size nodeIdx) {
        return nodes[nodeIdx];
    }

    /// \brief Returns the node with given index
    const TNode& getNode(const Size nodeIdx) const {
        return nodes[nodeIdx];
    }

    /// \brief Returns the number of nodes in the tree
    Size getNodeCnt() const {
        return nodes.size();
    }

    /// \brief Returns the sequence of particles indices belonging to given leaf.
    LeafIndexSequence getLeafIndices(const LeafNode<TNode>& leaf) const {
        return LeafIndexSequence(leaf.from, leaf.to, idxs);
    }*/

    bool sanityCheck() const;

private:
    void init() {
        entireBox_ = Box();
        values_.clear();
        idxs_.clear();
        nodes_.clear();
        nodeCounter_ = 0;
    }

    void buildTree(const Index parent,
        const KdChild child,
        const Index from,
        const Index to,
        const Box& box,
        const Index slidingCnt,
        const Index depth) {
        Box box1, box2;
        Vec boxSize = box.size();

        // split by the dimension of largest extent
        Index splitIdx = argMax(boxSize);

        bool slidingMidpoint = false;
        bool degeneratedBox = false;

        if (to - from <= config_.leafSize) {
            // enough points to fit inside one leaf
            addLeaf(parent, child, from, to);
            return;
        } else {
            // check for singularity of dimensions
            for (int dim = 0; dim < Dim; ++dim) {
                if (isSingular(from, to, splitIdx)) {
                    boxSize[splitIdx] = 0.f;
                    // find new largest dimension
                    splitIdx = argMax(boxSize);

                    if (boxSize == Vec(0)) {
                        // too many overlapping points, just split until they fit within a leaf,
                        // the code can handle this case, but it smells with an error ...
                        PVL_ASSERT(false, "Too many overlapping points, something is probably wrong ...");
                        degeneratedBox = true;
                        break;
                    }
                } else {
                    break;
                }
            }

            // split around center of the box
            Float splitPosition = box.center()[splitIdx];
            IndexDiff n1 = from, n2 = to - 1;

            if (slidingCnt <= 5 && !degeneratedBox) {
                for (;; std::swap(idxs_[n1], idxs_[n2])) {
                    for (; n1 < IndexDiff(to) && values_[idxs_[n1]][splitIdx] <= splitPosition; ++n1)
                        ;
                    for (; n2 >= IndexDiff(from) && values_[idxs_[n2]][splitIdx] >= splitPosition; --n2)
                        ;
                    if (n1 >= n2) {
                        break;
                    }
                }

                if (n1 == IndexDiff(from)) {
                    Index idx = from;
                    splitPosition = values_[idxs_[from]][splitIdx];
                    for (Index i = from + 1; i < to; ++i) {
                        const Float x1 = values_[idxs_[i]][splitIdx];
                        if (x1 < splitPosition) {
                            idx = i;
                            splitPosition = x1;
                        }
                    }
                    std::swap(idxs_[from], idxs_[idx]);
                    n1++;
                    slidingMidpoint = true;
                } else if (n1 == IndexDiff(to)) {
                    Index idx = from;
                    splitPosition = values_[idxs_[from]][splitIdx];
                    for (Index i = from + 1; i < to; ++i) {
                        const Float x2 = values_[idxs_[i]][splitIdx];
                        if (x2 > splitPosition) {
                            idx = i;
                            splitPosition = x2;
                        }
                    }
                    std::swap(idxs_[to - 1], idxs_[idx]);
                    n1--;
                    slidingMidpoint = true;
                }

                std::tie(box1, box2) = splitBox(box, splitIdx, splitPosition);
            } else {
                n1 = (from + to) >> 1;
                // do quick select to sort elements around the midpoint
                typename std::vector<Index>::iterator iter = idxs_.begin();
                if (!degeneratedBox) {
                    std::nth_element(iter + from, iter + n1, iter + to, [this, splitIdx](Index i1, Index i2) {
                        return values_[i1][splitIdx] < values_[i2][splitIdx];
                    });
                }

                std::tie(box1, box2) = splitBox(box, splitIdx, values_[idxs_[n1]][splitIdx]);
            }

            // sanity check
            PVL_ASSERT(checkBoxes(from, to, n1, box1, box2));

            // add inner node and connect it to the parent
            const Index index = addInner(parent, child, splitPosition, splitIdx);

            // recurse to left and right subtree
            const Index nextSlidingCnt = slidingMidpoint ? slidingCnt + 1 : 0;
            // auto processRightSubTree = [this, &scheduler, index, to, n1, box2, nextSlidingCnt, depth] {
            buildTree(index, KdChild::RIGHT, n1, to, box2, nextSlidingCnt, depth + 1);
            //};
            /*if (depth < config.maxParallelDepth) {
                // ad hoc decision - split the build only for few topmost nodes, there is no point in
                // splitting the work for child node in the bottom, it would only overburden the
            ThreadPool. scheduler.submit(processRightSubTree); } else {
                // otherwise simply process both subtrees in the same thread
                processRightSubTree();
            }*/
            buildTree(index, KdChild::LEFT, from, n1, box1, nextSlidingCnt, depth + 1);
        }
    }

    void addLeaf(const Index parent, const KdChild child, const Index from, const Index to) {
        const Index index = nodeCounter_++;
        if (index >= nodes_.size()) {
            // needs more nodes than estimated; allocate up to 2x more than necessary to avoid frequent
            // reallocations
            // nodesMutex.lock();
            nodes_.resize(std::max<Index>(2 * index, nodes_.size()));
            // nodesMutex.unlock();
        }

        // nodesMutex.lock_shared();
        // auto releaseLock = finally([this] { nodesMutex.unlock_shared(); });

        LeafNode& node = (LeafNode&)nodes_[index];
        node.setLeaf();
        PVL_ASSERT(node.isLeaf());

        node.from = node.to = -1;

        node.from = from;
        node.to = to;

        // find the bounding box of the leaf
        Box box;
        for (Index i = from; i < to; ++i) {
            box.extend(values_[idxs_[i]]);
        }
        node.box = box;

        if (parent == ROOT_PARENT_NODE) {
            return;
        }
        InnerNode& parentNode = (InnerNode&)nodes_[parent];
        PVL_ASSERT(!parentNode.isLeaf());
        if (child == KdChild::LEFT) {
            // left child
            parentNode.left = index;
        } else {
            PVL_ASSERT(child == KdChild::RIGHT);
            // right child
            parentNode.right = index;
        }
    }

    Index addInner(const Index parent, const KdChild child, const Float splitPosition, const Index splitIdx) {
        /*static_assert(int(KdNode::Type::X) == 0 && int(KdNode::Type::Y) == 1 && int(KdNode::Type::Z) == 2,
            "Invalid values of KdNode::Type enum");*/

        const Index index = nodeCounter_++;
        if (index >= nodes_.size()) {
            // needs more nodes than estimated; allocate up to 2x more than necessary to avoid frequent
            // reallocations
            // nodesMutex.lock();
            nodes_.resize(std::max<Index>(2 * index, nodes_.size()));
            // nodesMutex.unlock();
        }

        // nodesMutex.lock_shared();
        // auto releaseLock = finally([this] { nodesMutex.unlock_shared(); });
        InnerNode& node = (InnerNode&)nodes_[index];
        node.type = splitIdx;
        PVL_ASSERT(!node.isLeaf());

        node.left = node.right = -1;
        node.box = Box(); // will be computed later

        node.splitPosition = splitPosition;

        if (parent == ROOT_PARENT_NODE) {
            // no need to set up parents
            return index;
        }
        InnerNode& parentNode = (InnerNode&)nodes_[parent];
        if (child == KdChild::LEFT) {
            // left child
            PVL_ASSERT(parentNode.left == Index(-1));
            parentNode.left = index;
        } else {
            PVL_ASSERT(child == KdChild::RIGHT);
            // right child
            PVL_ASSERT(parentNode.right == Index(-1));
            parentNode.right = index;
        }

        return index;
    }

    bool isSingular(const Index from, const Index to, const Index splitIdx) const {
        for (Index i = from; i < to; ++i) {
            if (values_[idxs_[i]][splitIdx] != values_[idxs_[to - 1]][splitIdx]) {
                return false;
            }
        }
        return true;
    }


    bool checkBoxes(const Index from,
        const Index to,
        const Index mid,
        const Box& box1,
        const Box& box2) const {
        for (Index i = from; i < to; ++i) {
            if (i < mid && !box1.contains(values_[idxs_[i]])) {
                return false;
            }
            if (i >= mid && !box2.contains(values_[idxs_[i]])) {
                return false;
            }
        }
        return true;
    }

    /*    bool checkBoxes(const Size from,
            const Size to,
            const Size mid,
            const BoundingBox& box1,
            const BoundingBox& box2) const;*/
};


enum class IterateDirection {
    TOP_DOWN,  ///< From root to leaves
    BOTTOM_UP, ///< From leaves to root
};

/// \brief Calls a functor for every node of a K-d tree tree in specified direction.
///
/// The functor is called with the node as a parameter. For top-down direction, functor may return false
/// to skip all children nodes from processing, otherwise the iteration proceedes through the tree into
/// leaf nodes.
/// \param tree KdTree to iterate.
/// \param scheduler Scheduler used for sequential or parallelized task execution
/// \param functor Functor executed for every node
/// \param nodeIdx Index of the first processed node; use 0 for root node
/// \param depthLimit Maximal depth processed in parallel.
/*template <IterateDirection Dir, typename TNode, typename TMetric, typename TFunctor>
void iterateTree(KdTree<TNode, TMetric>& tree,
    IScheduler& scheduler,
    const TFunctor& functor,
    const Size nodeIdx = 0,
    const Size depthLimit = Size(-1));

/// \copydoc iterateTree
template <IterateDirection Dir, typename TNode, typename TMetric, typename TFunctor>
void iterateTree(const KdTree<TNode, TMetric>& tree,
    IScheduler& scheduler,
    const TFunctor& functor,
    const Size nodeIdx = 0,
    const Size depthLimit = Size(-1));
*/
} // namespace Pvl

#include "KdTree.inl.hpp"
