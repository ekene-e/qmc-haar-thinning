// SPDX-License-Identifier: MIT
#pragma once

/// @file haar.hpp
/// Umbrella header for the Haar-thinning library.
///
///   #include <haar/haar.hpp>
///
///   haar::Options opt;
///   opt.dim = 3;
///   std::vector<double> pts = haar::thin(100000, opt);   // 100000 x 3, row-major
///
/// See `thinner.hpp` for the online interface and `discrepancy.hpp` for
/// quality measures.

#include "haar/discrepancy.hpp"
#include "haar/options.hpp"
#include "haar/rng.hpp"
#include "haar/thinner.hpp"
#include "haar/version.hpp"
