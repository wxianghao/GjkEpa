#pragma once

#include <gjkepa.hpp>

constexpr Real epsilon = Eigen::NumTraits<Real>::epsilon();
constexpr Real epsRel  = epsilon * 1e4;
constexpr Real epsAbs  = epsilon * 1e2;