#include "isa/Metrics.h"
#include "isa/Polynomial.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <numeric>
#include <set>
#include <utility>
#include <vector>

namespace isa {
namespace {
struct DisjointSet {
    std::vector<size_t> parent,rank;
    explicit DisjointSet(size_t n):parent(n),rank(n,0){std::iota(parent.begin(),parent.end(),0);}
    size_t find(size_t x){while(parent[x]!=x){parent[x]=parent[parent[x]];x=parent[x];}return x;}
    void unite(size_t a,size_t b){a=find(a);b=find(b);if(a==b)return;if(rank[a]<rank[b])std::swap(a,b);parent[b]=a;if(rank[a]==rank[b])++rank[a];}
};
using Edge=std::pair<std::uint32_t,std::uint32_t>;
Edge edgeKey(std::uint32_t a,std::uint32_t b){return a<b?Edge{a,b}:Edge{b,a};}
double pointTriangleDistance2(const Vec3& p,const Vec3& a,const Vec3& b,const Vec3& c) {
    // Closest-point regions from the standard point-triangle distance construction.
    const Vec3 ab=b-a,ac=c-a,ap=p-a; const double d1=dot(ab,ap),d2=dot(ac,ap);
    if(d1<=0.0&&d2<=0.0)return norm2(ap);
    const Vec3 bp=p-b; const double d3=dot(ab,bp),d4=dot(ac,bp);
    if(d3>=0.0&&d4<=d3)return norm2(bp);
    const double vc=d1*d4-d3*d2;
    if(vc<=0.0&&d1>=0.0&&d3<=0.0){const double v=d1/(d1-d3);return norm2(p-(a+ab*v));}
    const Vec3 cp=p-c; const double d5=dot(ab,cp),d6=dot(ac,cp);
    if(d6>=0.0&&d5<=d6)return norm2(cp);
    const double vb=d5*d2-d1*d6;
    if(vb<=0.0&&d2>=0.0&&d6<=0.0){const double w=d2/(d2-d6);return norm2(p-(a+ac*w));}
    const double va=d3*d6-d5*d4;
    if(va<=0.0&&(d4-d3)>=0.0&&(d5-d6)>=0.0){const double w=(d4-d3)/((d4-d3)+(d5-d6));return norm2(p-(b+(c-b)*w));}
    const double denom=va+vb+vc;
    if(std::abs(denom)<1e-30) return std::min({norm2(ap),norm2(bp),norm2(cp)});
    const double inv=1.0/denom,v=vb*inv,w=vc*inv;
    return norm2(p-(a+ab*v+ac*w));
}
std::vector<size_t> evenlySpacedIndices(size_t total,size_t maxCount) {
    std::vector<size_t> out; if(total==0||maxCount==0)return out;
    const size_t count=std::min(total,maxCount); out.reserve(count);
    if(count==1){out.push_back(total/2);return out;}
    for(size_t i=0;i<count;++i)out.push_back(i*(total-1)/(count-1));
    return out;
}
double directedVertexToTriangleRms(const Mesh& source,const Mesh& target,size_t maxVertices,size_t maxTriangles) {
    if(source.vertices.empty()||target.triangles.empty())return std::numeric_limits<double>::infinity();
    const auto vi=evenlySpacedIndices(source.vertices.size(),maxVertices);
    const auto ti=evenlySpacedIndices(target.triangles.size(),maxTriangles);
    double sum=0.0;
    for(size_t id:vi){
        const Vec3& p=source.vertices[id]; double best=std::numeric_limits<double>::infinity();
        for(size_t tid:ti){const auto& t=target.triangles[tid];
            best=std::min(best,pointTriangleDistance2(p,target.vertices[t[0]],target.vertices[t[1]],target.vertices[t[2]]));
            if(best<1e-24)break;
        }
        sum+=best;
    }
    return std::sqrt(sum/static_cast<double>(vi.size()));
}
}
MeshMetrics computeMeshMetrics(const Mesh& mesh,const Polynomial* referenceField) {
    MeshMetrics result;result.vertices=mesh.vertices.size();result.triangles=mesh.triangles.size();
    std::map<Edge,size_t> edgeUses; DisjointSet dsu(mesh.vertices.size()); std::vector<bool> used(mesh.vertices.size(),false);
    const double bboxScale=1.0;
    for(const auto& tri:mesh.triangles){
        const auto a=tri[0],b=tri[1],c=tri[2];
        if(a>=mesh.vertices.size()||b>=mesh.vertices.size()||c>=mesh.vertices.size())continue;
        used[a]=used[b]=used[c]=true;
        dsu.unite(a,b);dsu.unite(b,c);dsu.unite(c,a);
        edgeUses[edgeKey(a,b)]++;edgeUses[edgeKey(b,c)]++;edgeUses[edgeKey(c,a)]++;
        const double area2=norm(cross(mesh.vertices[b]-mesh.vertices[a],mesh.vertices[c]-mesh.vertices[a]));
        if(!std::isfinite(area2)||area2<1e-12*bboxScale*bboxScale)++result.degenerateTriangles;
    }
    result.edges=edgeUses.size();
    for(const auto& [edge,count]:edgeUses){(void)edge;if(count==1)++result.boundaryEdges;else if(count>2)++result.nonManifoldEdges;}
    std::set<size_t> roots;size_t usedCount=0;
    for(size_t i=0;i<used.size();++i)if(used[i]){++usedCount;roots.insert(dsu.find(i));}
    result.connectedComponents=roots.size();
    result.eulerCharacteristic=static_cast<long long>(usedCount)-static_cast<long long>(result.edges)+static_cast<long long>(result.triangles);
    if(referenceField && !mesh.vertices.empty()){
        const PolynomialModel model(*referenceField); const auto samples=evenlySpacedIndices(mesh.vertices.size(),4000);
        double sum=0.0;size_t count=0;
        double maxGrad=1.0;
        for(const auto& p:mesh.vertices)maxGrad=std::max(maxGrad,norm(model.gradient(p)));
        const double gradTol=1e-10*maxGrad;
        for(size_t i:samples){const Vec3& p=mesh.vertices[i];const double g=norm(model.gradient(p));
            if(g<=gradTol){++result.skippedNearCriticalSamples;continue;}
            const double d=std::abs(model.evaluate(p))/g;
            if(std::isfinite(d)){sum+=d*d;++count;}
        }
        result.residualSamples=count;
        if(count>0)result.implicitResidualRms=std::sqrt(sum/static_cast<double>(count));
    }
    return result;
}
double approximateSymmetricMeshDistance(const Mesh& a,const Mesh& b,size_t maxSampleVertices,size_t maxSampleTriangles) {
    if(a.empty()&&b.empty())return 0.0;
    if(a.empty()||b.empty())return std::numeric_limits<double>::infinity();
    const double ab=directedVertexToTriangleRms(a,b,maxSampleVertices,maxSampleTriangles);
    const double ba=directedVertexToTriangleRms(b,a,maxSampleVertices,maxSampleTriangles);
    if(!std::isfinite(ab)||!std::isfinite(ba))return std::numeric_limits<double>::infinity();
    return std::sqrt((ab*ab+ba*ba)/2.0);
}
} // namespace isa
