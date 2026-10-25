#include "isa/Singularity.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <sstream>

namespace isa {
namespace {
bool solve3(std::array<std::array<double,3>,3> a,std::array<double,3> b,std::array<double,3>& x) {
    for (int col=0;col<3;++col) {
        int pivot=col;
        for (int r=col+1;r<3;++r) if (std::abs(a[static_cast<size_t>(r)][static_cast<size_t>(col)])>std::abs(a[static_cast<size_t>(pivot)][static_cast<size_t>(col)])) pivot=r;
        if (std::abs(a[static_cast<size_t>(pivot)][static_cast<size_t>(col)])<1e-13) return false;
        std::swap(a[static_cast<size_t>(pivot)],a[static_cast<size_t>(col)]);
        std::swap(b[static_cast<size_t>(pivot)],b[static_cast<size_t>(col)]);
        const double d=a[static_cast<size_t>(col)][static_cast<size_t>(col)];
        for (int j=col;j<3;++j) a[static_cast<size_t>(col)][static_cast<size_t>(j)]/=d;
        b[static_cast<size_t>(col)]/=d;
        for (int r=0;r<3;++r) if (r!=col) {
            const double factor=a[static_cast<size_t>(r)][static_cast<size_t>(col)];
            for (int j=col;j<3;++j) a[static_cast<size_t>(r)][static_cast<size_t>(j)]-=factor*a[static_cast<size_t>(col)][static_cast<size_t>(j)];
            b[static_cast<size_t>(r)]-=factor*b[static_cast<size_t>(col)];
        }
    }
    x=b; return true;
}
std::array<double,3> symmetricEigenvalues(std::array<std::array<double,3>,3> a) {
    // Jacobi rotations are stable enough for the symmetric 3x3 Hessians used here.
    for (int iter=0;iter<48;++iter) {
        int p=0,q=1; double largest=std::abs(a[0][1]);
        if (std::abs(a[0][2])>largest) { p=0;q=2;largest=std::abs(a[0][2]); }
        if (std::abs(a[1][2])>largest) { p=1;q=2;largest=std::abs(a[1][2]); }
        if (largest<1e-13) break;
        const double app=a[static_cast<size_t>(p)][static_cast<size_t>(p)];
        const double aqq=a[static_cast<size_t>(q)][static_cast<size_t>(q)];
        const double apq=a[static_cast<size_t>(p)][static_cast<size_t>(q)];
        const double phi=0.5*std::atan2(2.0*apq,aqq-app);
        const double c=std::cos(phi),s=std::sin(phi);
        for (int k=0;k<3;++k) if (k!=p && k!=q) {
            const double akp=a[static_cast<size_t>(k)][static_cast<size_t>(p)];
            const double akq=a[static_cast<size_t>(k)][static_cast<size_t>(q)];
            a[static_cast<size_t>(k)][static_cast<size_t>(p)]=a[static_cast<size_t>(p)][static_cast<size_t>(k)]=c*akp-s*akq;
            a[static_cast<size_t>(k)][static_cast<size_t>(q)]=a[static_cast<size_t>(q)][static_cast<size_t>(k)]=s*akp+c*akq;
        }
        a[static_cast<size_t>(p)][static_cast<size_t>(p)]=c*c*app-2*s*c*apq+s*s*aqq;
        a[static_cast<size_t>(q)][static_cast<size_t>(q)]=s*s*app+2*s*c*apq+c*c*aqq;
        a[static_cast<size_t>(p)][static_cast<size_t>(q)]=a[static_cast<size_t>(q)][static_cast<size_t>(p)]=0.0;
    }
    std::array<double,3> values{a[0][0],a[1][1],a[2][2]};
    std::sort(values.begin(),values.end()); return values;
}
std::string classify(const std::array<double,3>& eig) {
    const double largest=std::max({std::abs(eig[0]),std::abs(eig[1]),std::abs(eig[2]),1.0});
    const double tol=largest*1e-7; int pos=0,neg=0,zero=0;
    for (double x:eig) { if (x>tol) ++pos; else if (x< -tol) ++neg; else ++zero; }
    std::ostringstream out;
    if (zero>0) out<<"Degenerate critical point ("<<pos<<"+, "<<neg<<"-, "<<zero<<" near-zero Hessian eigenvalues)";
    else if (pos==3) out<<"Positive-definite Hessian (local minimum of F)";
    else if (neg==3) out<<"Negative-definite Hessian (local maximum of F)";
    else out<<"Indefinite Hessian ("<<pos<<"+, "<<neg<<"-)";
    if (zero==0 && pos>0 && neg>0) out<<"; ordinary quadratic-type candidate";
    return out.str();
}
bool inside(const Vec3& p,double extent,double margin=1.02) {
    const double b=extent*margin; return std::abs(p.x)<=b && std::abs(p.y)<=b && std::abs(p.z)<=b;
}
}
SingularitySearchResult findSingularityCandidates(const Polynomial& f,double extent,int seedCount) {
    SingularitySearchResult result;
    if (!(extent>0.0) || !std::isfinite(extent)) { result.diagnostic="Invalid domain extent."; return result; }
    seedCount=std::clamp(seedCount,3,9);
    const PolynomialModel model(f);
    double maxAbsF=1.0,maxGrad=1.0;
    for (int ix : {-1,1}) for (int iy : {-1,1}) for (int iz : {-1,1}) {
        const Vec3 p{ix*extent,iy*extent,iz*extent};
        maxAbsF=std::max(maxAbsF,std::abs(model.evaluate(p)));
        maxGrad=std::max(maxGrad,norm(model.gradient(p)));
    }
    const double fTol=1e-8*maxAbsF, gTol=1e-8*maxGrad;
    int convergedCriticalPoints=0;
    for (int iz=0;iz<seedCount;++iz) for (int iy=0;iy<seedCount;++iy) for (int ix=0;ix<seedCount;++ix) {
        Vec3 p{-extent+2.0*extent*ix/(seedCount-1),-extent+2.0*extent*iy/(seedCount-1),-extent+2.0*extent*iz/(seedCount-1)};
        bool converged=false;
        for (int it=0;it<40;++it) {
            const Vec3 g=model.gradient(p); const double gn=norm(g);
            if (gn<=gTol) { converged=true; break; }
            std::array<double,3> delta{};
            if (!solve3(model.hessianAt(p),{-g.x,-g.y,-g.z},delta)) break;
            Vec3 step{delta[0],delta[1],delta[2]};
            const double stepNorm=norm(step);
            if (!std::isfinite(stepNorm) || stepNorm>4.0*extent) break;
            // Mild damping makes Newton iterations less likely to jump out of the search domain.
            if (stepNorm>extent) step=step*(extent/stepNorm);
            p+=step;
            if (!inside(p,extent,1.5)) break;
        }
        if (!converged || !inside(p,extent) || norm(model.gradient(p))>gTol*10.0) continue;
        ++convergedCriticalPoints;
        const double residual=std::abs(model.evaluate(p));
        if (residual>fTol) continue; // critical point, but not on the requested isosurface
        bool duplicate=false;
        const double dedup=std::max(1e-7,extent*1e-5);
        for (const auto& old:result.points) if (norm(old.position-p)<dedup) { duplicate=true; break; }
        if (duplicate) continue;
        const auto eig=symmetricEigenvalues(model.hessianAt(p));
        result.points.push_back({p,residual,norm(model.gradient(p)),eig,classify(eig)});
    }
    std::ostringstream message;
    message<<"Grid-seeded Newton search: "<<convergedCriticalPoints<<" critical-point convergence(s), "<<result.points.size()<<" candidate singularity/singularities on F=0. This is a numerical search, not a completeness proof.";
    result.diagnostic=message.str(); return result;
}
} // namespace isa
