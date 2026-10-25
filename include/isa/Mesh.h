#pragma once
#include "isa/Vec3.h"
#include <array>
#include <cstdint>
#include <vector>

namespace isa {
struct Mesh {
    std::vector<Vec3> vertices;
    std::vector<std::array<std::uint32_t,3>> triangles;
    bool empty() const { return triangles.empty() || vertices.empty(); }
};
struct GridSpec {
    int cells{32};
    double extent{2.0}; // symmetric domain [-extent, extent]^3
};
} // namespace isa
