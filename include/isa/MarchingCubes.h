#pragma once
#include "isa/Mesh.h"
#include "isa/Polynomial.h"
#include <functional>

namespace isa {
Mesh extractMarchingCubes(const Polynomial& field, const GridSpec& grid,
                          const std::function<void(int,int)>& progress = {});
} // namespace isa
