// SPDX-License-Identifier: MIT
#pragma once

/// @file version.hpp
/// Library version, as a macro triple and as a string.

#define HAAR_VERSION_MAJOR 0
#define HAAR_VERSION_MINOR 1
#define HAAR_VERSION_PATCH 0
#define HAAR_VERSION_STRING "0.1.0"

namespace haar {

/// Semantic version of the library.
struct Version {
    int major;
    int minor;
    int patch;
};

/// Returns the version this header was compiled from.
[[nodiscard]] constexpr Version version() noexcept {
    return {HAAR_VERSION_MAJOR, HAAR_VERSION_MINOR, HAAR_VERSION_PATCH};
}

}  // namespace haar
