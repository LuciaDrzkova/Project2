#include "isa/MarchingCubes.h"
#include "isa/Metrics.h"
#include "isa/Polynomial.h"
#include "isa/Singularity.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
void require(bool condition,const std::string& message){if(!condition)throw std::runtime_error(message);}
void close(double actual,double expected,double tolerance,const std::string& message){if(std::abs(actual-expected)>tolerance)throw std::runtime_error(message+": expected "+std::to_string(expected)+", got "+std::to_string(actual));}
}
int main(){
    try {
        const auto sphere=isa::Polynomial::parse("x^2+y^2+z^2=1");
        close(sphere.evaluate({1,0,0}),0.0,1e-12,"Equation parser should move the right-hand side to F=0");
        const isa::PolynomialModel sm(sphere);
        const auto g=sm.gradient({1,2,3});close(g.x,2,1e-12,"Sphere F_x");close(g.y,4,1e-12,"Sphere F_y");close(g.z,6,1e-12,"Sphere F_z");
        const auto product=isa::Polynomial::parse("2x(y+z)-x*y");
        close(product.evaluate({1,2,3}),8.0,1e-12,"Parser should support implicit multiplication and parentheses");
        close(product.derivative(0).evaluate({1,2,3}),8.0,1e-12,"Symbolic derivative of parsed polynomial");
        const auto sphereSingularities=isa::findSingularityCandidates(sphere,1.5,5);
        require(sphereSingularities.points.empty(),"Sphere should have no singular point on its zero set");
        const auto cone=isa::Polynomial::parse("x^2+y^2-z^2");
        const auto coneSingularities=isa::findSingularityCandidates(cone,1.5,5);
        require(!coneSingularities.points.empty(),"Cone origin should be detected as a singularity candidate");
        const auto found=coneSingularities.points.front();
        close(found.position.x,0,1e-7,"Cone singular x");close(found.position.y,0,1e-7,"Cone singular y");close(found.position.z,0,1e-7,"Cone singular z");
        require(found.classification.find("Indefinite")!=std::string::npos,"Cone Hessian should be classified as indefinite");
        const auto hyper=isa::Polynomial::parse("x^2+y^2-z^2-0.2*x");
        const auto hyperSingularities=isa::findSingularityCandidates(hyper,1.5,5);
        require(hyperSingularities.points.empty(),"Shifted hyperboloid should be smooth on its zero set");
        const auto mesh=isa::extractMarchingCubes(sphere,{24,1.5});
        require(!mesh.empty(),"Marching Cubes should create a sphere mesh");
        const auto metrics=isa::computeMeshMetrics(mesh,&sphere);
        require(metrics.vertices>0 && metrics.triangles>0,"Mesh metrics should count vertices and triangles");
        require(metrics.degenerateTriangles==0,"Sphere mesh should not contain degenerate triangles at this resolution");
        require(metrics.connectedComponents==1,"Sphere mesh should have one connected component");
        require(std::isfinite(metrics.implicitResidualRms)&&metrics.implicitResidualRms<0.05,"Sphere residual should be small at this resolution");
        const double same=isa::approximateSymmetricMeshDistance(mesh,mesh);
        close(same,0.0,1e-12,"Mesh distance to itself");
        const auto perturbed=isa::Polynomial::parse("x^2+y^2-z^2-0.02*x");
        const auto perturbedMesh=isa::extractMarchingCubes(perturbed,{24,1.5});
        require(!perturbedMesh.empty(),"Perturbed cone should generate an isosurface mesh");
        std::cout<<"All core tests passed.\n";
        return 0;
    }catch(const std::exception& e){std::cerr<<"Test failed: "<<e.what()<<'\n';return 1;}
}
