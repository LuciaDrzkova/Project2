#include "isa/PmpExport.h"
#include <pmp/io/io.h>
#include <pmp/surface_mesh.h>
#include <stdexcept>
#include <fstream>
#include <iomanip>
#include <vector>

namespace isa {
void exportWithPmp(const Mesh& source,const std::filesystem::path& path) {
    if(source.vertices.empty()||source.triangles.empty()) throw std::runtime_error("There is no non-empty mesh to export");
    pmp::SurfaceMesh mesh;
    std::vector<pmp::Vertex> vertices; vertices.reserve(source.vertices.size());
    for(const auto& p:source.vertices) vertices.push_back(mesh.add_vertex(pmp::Point(p.x,p.y,p.z)));
    size_t rejected=0;
    for(const auto& tri:source.triangles) {
        if(tri[0]>=vertices.size()||tri[1]>=vertices.size()||tri[2]>=vertices.size()) { ++rejected; continue; }
        try {
            const auto face=mesh.add_triangle(vertices[tri[0]],vertices[tri[1]],vertices[tri[2]]);
            if(!face.is_valid()) ++rejected;
        } catch(const std::exception&) {
            ++rejected;
        }
    }
    if(rejected>0) throw std::runtime_error("PMP SurfaceMesh rejected "+std::to_string(rejected)+" face(s). This can happen for non-manifold singular meshes. No lossy PMP file was written; use the raw OBJ export to preserve the generated triangle list.");
    if(mesh.n_faces()==0) throw std::runtime_error("PMP produced an empty mesh");
    pmp::write(mesh,path);
}
void exportRawObj(const Mesh& mesh,const std::filesystem::path& path) {
    if(mesh.vertices.empty()||mesh.triangles.empty()) throw std::runtime_error("There is no non-empty mesh to export");
    std::ofstream out(path);
    if(!out) throw std::runtime_error("Could not open output file for writing");
    out<<"# Singularity-aware isosurface prototype raw triangle mesh\n"<<std::setprecision(17);
    for(const auto& p:mesh.vertices)out<<"v "<<p.x<<' '<<p.y<<' '<<p.z<<'\n';
    for(const auto& t:mesh.triangles){
        if(t[0]>=mesh.vertices.size()||t[1]>=mesh.vertices.size()||t[2]>=mesh.vertices.size())continue;
        out<<"f "<<static_cast<unsigned long long>(t[0])+1<<' '<<static_cast<unsigned long long>(t[1])+1<<' '<<static_cast<unsigned long long>(t[2])+1<<'\n';
    }
    if(!out)throw std::runtime_error("An error occurred while writing the OBJ file");
}
} // namespace isa
