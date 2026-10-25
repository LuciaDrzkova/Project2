#pragma once
#include "isa/Mesh.h"
#include "isa/Polynomial.h"
#include <cstddef>
#include <limits>

namespace isa {
struct MeshMetrics {
    size_t vertices{0}, triangles{0}, edges{0}, connectedComponents{0};
    long long eulerCharacteristic{0};
    size_t boundaryEdges{0}, nonManifoldEdges{0}, degenerateTriangles{0};
    double implicitResidualRms{std::numeric_limits<double>::quiet_NaN()};
    size_t residualSamples{0}, skippedNearCriticalSamples{0};
};
MeshMetrics computeMeshMetrics(const Mesh& mesh,const Polynomial* referenceField=nullptr);
// Approximate symmetric mesh distance using sampled vertices and sampled target triangles.
// This is a convergence diagnostic, not an exact continuous Hausdorff distance.
double approximateSymmetricMeshDistance(const Mesh& a,const Mesh& b,size_t maxSampleVertices=250,size_t maxSampleTriangles=8000);
} // namespace isa
