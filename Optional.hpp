#pragma once
#include "Assert.hpp"
#include <type_traits>

namespace Pvl {

struct OptionalNone {};
const OptionalNone NONE;

template <typename T>
class Optional {
    typename std::aligned_storage<sizeof(T), alignof(T)>::type data_;
    bool used_;

public:
    Optional()
        : used_(false) {}

    Optional(OptionalNone)
        : used_(false) {}

    Optional(const T& value)
        : used_(true) {
        new (&data_) T(value);
    }

    Optional(T&& value)
        : used_(true) {
        new (&data_) T(std::move(value));
    }

    ~Optional() {
        if (used_) {
            value().~T();
        }
    }

    T& value() {
        PVL_ASSERT(used_);
        return reinterpret_cast<T&>(data_);
    }

    const T& value() const {
        PVL_ASSERT(used_);
        return reinterpret_cast<const T&>(data_);
    }

    template <typename Alt>
    T valueOr(Alt&& alt) {
        if (used_) {
            return value();
        } else {
            return std::forward<Alt>(alt);
        }
    }

    template <typename Exception, typename... TArgs>
    T valueOrThrow(TArgs&&... args) {
        if (used_) {
            return value();
        } else {
            throw Exception(std::forward<TArgs>(args)...);
        }
    }

    T* operator->() const {
        PVL_ASSERT(used_);
        return &value();
    }

    explicit operator bool() const {
        return used_;
    }

    bool operator!() const {
        return !used_;
    }
};

template <typename T>
bool hasValue(const T&) {
    return true;
}
template <typename T>
bool hasValue(const Optional<T>& opt) {
    return bool(opt);
}
template <typename T>
const T& value(const T& v) {
    return v;
}
template <typename T>
const T& value(const Optional<T>& opt) {
    return opt.value();
}


} // namespace Pvl
