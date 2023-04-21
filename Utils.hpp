#pragma once
#ifdef emit
#undef emit
#define mpcv_redefine_emit
#endif
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#ifdef mpcv_redefine_emit
#define emit
#undef mpcv_redefine_emit
#endif

namespace Pvl {

struct Noncopyable {
    Noncopyable() = default;
    Noncopyable(const Noncopyable&) = delete;
    Noncopyable(Noncopyable&&) = delete;
    Noncopyable& operator=(const Noncopyable&) = delete;
    Noncopyable& operator=(Noncopyable&&) = delete;
};

struct ParallelTag {};
struct SequentialTag {};

template <typename Tag>
struct ParallelFor;

template <>
struct ParallelFor<SequentialTag> {
    template <typename Index, typename Func>
    void operator()(Index n1, Index n2, const Func& func) {
        for (Index n = n1; n < n2; ++n) {
            func(n);
        }
    }
};

template <>
struct ParallelFor<ParallelTag> {
    template <typename Index, typename Func>
    void operator()(Index n1, Index n2, const Func& func) {
        tbb::parallel_for(n1, n2, func);
    }
};


template <typename Tag>
struct ParallelForEach;

template <>
struct ParallelForEach<SequentialTag> {
    template <typename Iter, typename Func>
    void operator()(Iter from, Iter to, const Func& func) {
        for (Iter iter = from; iter != to; ++iter) {
            func(*iter);
        }
    }
    template <typename Range, typename Func>
    void operator()(const Range& range, const Func& func) {
        for (auto&& value : range) {
            func(value);
        }
    }
};

template <>
struct ParallelForEach<ParallelTag> {
    template <typename Iter, typename Func>
    void operator()(Iter from, Iter to, const Func& func) {
        tbb::parallel_for_each(from, to, func);
    }
    template <typename Range, typename Func>
    void operator()(const Range& range, const Func& func) {
        tbb::parallel_for_each(range, func);
    }
};

template <typename Progress>
class ProgressMeter {
    Progress func_;
    std::atomic<int> counter_;
    int target_;
    int step_;
    int next_;
    std::thread::id callingThreadId_;

public:
    ProgressMeter(int target, Progress&& func)
        : func_(std::move(func))
        , counter_(0) {
        step_ = std::max(target / 100, 10);
        next_ = step_;
        target_ = target;
        callingThreadId_ = std::this_thread::get_id();
    }

    bool inc() {
        counter_++;
        if (std::this_thread::get_id() == callingThreadId_ && counter_ > next_) {
            float value = float(counter_) / target_ * 100;
            next_ += step_;
            return func_(value);
        } else {
            return false;
        }
    }
};

template <typename Progress>
ProgressMeter<Progress> makeProgressMeter(int target, Progress&& func) {
    return ProgressMeter<Progress>(target, std::move(func));
}

} // namespace Pvl
