#pragma once
#include "isa/Polynomial.h"
#include <array>
#include <string>
#include <vector>

namespace isa {
struct SingularityCandidate {
    Vec3 position;
    double functionResidual{0.0};
    double gradientNorm{0.0};
    std::array<double,3> hessianEigenvalues{};
    std::string classification;
};
struct SingularitySearchResult {
    std::vector<SingularityCandidate> points;
    std::string diagnostic;
};
SingularitySearchResult findSingularityCandidates(const Polynomial& f,double extent,int seedCount=5);
} // namespace isa
