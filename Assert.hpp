#pragma once

#include <cassert>

#ifndef NDEBUG
#define PVL_ASSERT(x, ...) assert(x)
#else
#define PVL_ASSERT(x, ...) (void)sizeof(x)
#endif
