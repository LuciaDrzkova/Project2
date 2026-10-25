#pragma once
#include "isa/Vec3.h"
#include <array>
#include <map>
#include <string>
#include <vector>

namespace isa {
using Exponents = std::array<int,3>;

class Polynomial {
public:
    using Terms = std::map<Exponents,double>;
    Polynomial() = default;
    explicit Polynomial(double constant);
    static Polynomial variable(int axis);
    static Polynomial parse(const std::string& expression);

    double evaluate(const Vec3& p) const;
    Polynomial derivative(int axis) const;
    int degree() const;
    bool isZero(double epsilon=1e-14) const;
    const Terms& terms() const { return terms_; }

    Polynomial operator+(const Polynomial& rhs) const;
    Polynomial operator-(const Polynomial& rhs) const;
    Polynomial operator*(const Polynomial& rhs) const;
    Polynomial operator*(double scalar) const;
    Polynomial pow(unsigned exponent) const;

private:
    Terms terms_;
    explicit Polynomial(Terms terms);
    void clean(double epsilon=1e-14);
    friend class PolynomialParser;
};

struct PolynomialModel {
    explicit PolynomialModel(Polynomial polynomial);
    Polynomial f;
    Polynomial fx, fy, fz;
    std::array<std::array<Polynomial,3>,3> hessian;
    double evaluate(const Vec3& p) const { return f.evaluate(p); }
    Vec3 gradient(const Vec3& p) const { return {fx.evaluate(p), fy.evaluate(p), fz.evaluate(p)}; }
    std::array<std::array<double,3>,3> hessianAt(const Vec3& p) const;
};
} // namespace isa
