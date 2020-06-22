#pragma once

namespace Pvl {

template <typename T, typename TDerived>
class IterBase {
protected:
    T iter_;

public:
    IterBase(const T& iter)
        : iter_(iter) {}

    TDerived& operator++() {
        ++iter_;
        return static_cast<TDerived&>(*this);
    }

    bool operator==(const TDerived& other) const {
        return iter_ == other.iter_;
    }

    bool operator!=(const TDerived& other) const {
        return iter_ != other.iter_;
    }
};

template <typename TIter>
class Range {
    TIter begin_;
    TIter end_;

public:
    Range(const TIter& begin, const TIter& end)
        : begin_(begin)
        , end_(end) {}

    TIter begin() const {
        return begin_;
    }

    TIter end() const {
        return end_;
    }
};

} // namespace Pvl
