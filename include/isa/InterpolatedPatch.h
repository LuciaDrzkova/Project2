#pragma once

#include "isa/Mesh.h"

namespace isa {

// Builds a linear interpolant between the negative-k and positive-k surfaces.
// The reference mesh supplies the triangle connectivity/correspondence template;
// each of its vertices is projected to both perturbed meshes by exact
// closest-point queries on triangle faces, then blended linearly.
// alpha=0 returns the exact negativeK mesh; alpha=1 returns the exact positiveK mesh.
// Intermediate alpha values use the reference connectivity and project reference
// vertices to both target meshes before linearly blending positions. This is a
// geometric approximation and does not guarantee topology preservation.
Mesh makeLinearInterpolatedPatch(const Mesh& reference,
                                 const Mesh& negativeK,
                                 const Mesh& positiveK,
                                 double alpha = 0.5);

} // namespace isa
