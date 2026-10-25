#pragma once
#include "isa/Mesh.h"
#include <filesystem>

namespace isa {
// Converts the application's triangle mesh to PMP::SurfaceMesh and saves using PMP I/O.
void exportWithPmp(const Mesh& mesh,const std::filesystem::path& path);
void exportRawObj(const Mesh& mesh,const std::filesystem::path& path);
}
