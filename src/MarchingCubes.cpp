#include "isa/MarchingCubes.h"
#include "isa/McTables.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace isa {
namespace {
constexpr std::array<std::array<int,3>,8> kCornerOffsets{{
    {{0,0,0}}, {{1,0,0}}, {{1,1,0}}, {{0,1,0}},
    {{0,0,1}}, {{1,0,1}}, {{1,1,1}}, {{0,1,1}}
}};
constexpr std::array<std::array<int,2>,12> kEdgeCorners{{
    {{0,1}},{{1,2}},{{2,3}},{{3,0}},{{4,5}},{{5,6}},{{6,7}},{{7,4}},
    {{0,4}},{{1,5}},{{2,6}},{{3,7}}
}};
size_t gridIndex(int i,int j,int k,int n) {
    const size_t side=static_cast<size_t>(n+1);
    return (static_cast<size_t>(k)*side+static_cast<size_t>(j))*side+static_cast<size_t>(i);
}
Vec3 interpolate(const Vec3& a,const Vec3& b,double fa,double fb) {
    const double d=fa-fb;
    double t=std::abs(d)<1e-30?0.5:fa/d;
    t=std::clamp(t,0.0,1.0);
    return a+(b-a)*t;
}
}
Mesh extractMarchingCubes(const Polynomial& field,const GridSpec& grid,const std::function<void(int,int)>& progress) {
    if (grid.cells<2 || grid.cells>512) throw std::invalid_argument("Grid cell count must be between 2 and 512");
    if (!(grid.extent>0.0) || !std::isfinite(grid.extent)) throw std::invalid_argument("Domain extent must be finite and positive");
    const int n=grid.cells; const size_t side=static_cast<size_t>(n+1);
    if (side*side*side>80'000'000ULL) throw std::runtime_error("Grid is too large for this prototype");
    const double h=(2.0*grid.extent)/static_cast<double>(n);
    std::vector<double> values(side*side*side);
    for (int k=0;k<=n;++k) {
        const double z=-grid.extent+h*k;
        for (int j=0;j<=n;++j) {
            const double y=-grid.extent+h*j;
            for (int i=0;i<=n;++i) {
                const double x=-grid.extent+h*i;
                values[gridIndex(i,j,k,n)]=field.evaluate({x,y,z});
            }
        }
    }
    Mesh mesh;
    // Reuse the vertex generated for each grid edge so adjacent cells share vertices.
    std::map<std::pair<size_t,size_t>,std::uint32_t> edgeVertices;
    for (int k=0;k<n;++k) {
        if (progress && (k%std::max(1,n/20)==0)) progress(k,n);
        for (int j=0;j<n;++j) for (int i=0;i<n;++i) {
            std::array<Vec3,8> p{}; std::array<double,8> v{}; std::array<size_t,8> ids{};
            int cubeCase=0;
            for (int c=0;c<8;++c) {
                const int gi=i+kCornerOffsets[static_cast<size_t>(c)][0];
                const int gj=j+kCornerOffsets[static_cast<size_t>(c)][1];
                const int gk=k+kCornerOffsets[static_cast<size_t>(c)][2];
                ids[static_cast<size_t>(c)]=gridIndex(gi,gj,gk,n);
                p[static_cast<size_t>(c)]={-grid.extent+h*gi,-grid.extent+h*gj,-grid.extent+h*gk};
                v[static_cast<size_t>(c)]=values[ids[static_cast<size_t>(c)]];
                if (v[static_cast<size_t>(c)]<0.0) cubeCase|=(1<<c);
            }
            const auto& row=kMcTriangleTable[static_cast<size_t>(cubeCase)];
            if (row[0]<0) continue;
            std::array<std::uint32_t,12> edgeVertex{};
            std::array<bool,12> wasMade{};
            auto getEdgeVertex=[&](int edge)->std::uint32_t {
                const size_t e=static_cast<size_t>(edge);
                if (wasMade[e]) return edgeVertex[e];
                const int c0=kEdgeCorners[e][0],c1=kEdgeCorners[e][1];
                size_t a=ids[static_cast<size_t>(c0)],b=ids[static_cast<size_t>(c1)];
                const std::pair<size_t,size_t> key{std::min(a,b),std::max(a,b)};
                auto it=edgeVertices.find(key);
                std::uint32_t index;
                if (it!=edgeVertices.end()) index=it->second;
                else {
                    if (mesh.vertices.size()>=static_cast<size_t>(UINT32_MAX)) throw std::runtime_error("Mesh exceeds 32-bit index capacity");
                    index=static_cast<std::uint32_t>(mesh.vertices.size());
                    mesh.vertices.push_back(interpolate(p[static_cast<size_t>(c0)],p[static_cast<size_t>(c1)],v[static_cast<size_t>(c0)],v[static_cast<size_t>(c1)]));
                    edgeVertices.emplace(key,index);
                }
                edgeVertex[e]=index; wasMade[e]=true; return index;
            };
            for (int t=0;t<16 && row[static_cast<size_t>(t)]>=0;t+=3) {
                if (t+2>=16 || row[static_cast<size_t>(t+1)]<0 || row[static_cast<size_t>(t+2)]<0) break;
                const auto a=getEdgeVertex(row[static_cast<size_t>(t)]);
                const auto b=getEdgeVertex(row[static_cast<size_t>(t+1)]);
                const auto c=getEdgeVertex(row[static_cast<size_t>(t+2)]);
                if (a!=b && b!=c && a!=c) mesh.triangles.push_back({a,b,c});
            }
        }
    }
    if (progress) progress(n,n);
    return mesh;
}
} // namespace isa
