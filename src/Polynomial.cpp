#include "isa/Polynomial.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string_view>

namespace isa {
Polynomial::Polynomial(double constant) { if (std::abs(constant)>1e-14) terms_[{0,0,0}]=constant; }
Polynomial::Polynomial(Terms terms) : terms_(std::move(terms)) { clean(); }
Polynomial Polynomial::variable(int axis) {
    if (axis<0 || axis>2) throw std::invalid_argument("Variable axis must be 0, 1, or 2");
    Polynomial p; Exponents e{0,0,0}; e[static_cast<size_t>(axis)]=1; p.terms_[e]=1.0; return p;
}
void Polynomial::clean(double epsilon) {
    for (auto it=terms_.begin(); it!=terms_.end();) {
        if (std::abs(it->second)<epsilon) it=terms_.erase(it); else ++it;
    }
}
double Polynomial::evaluate(const Vec3& p) const {
    double result=0.0;
    for (const auto& [e,c] : terms_) {
        result += c*std::pow(p.x,e[0])*std::pow(p.y,e[1])*std::pow(p.z,e[2]);
    }
    return result;
}
Polynomial Polynomial::derivative(int axis) const {
    if (axis<0 || axis>2) throw std::invalid_argument("Derivative axis must be 0, 1, or 2");
    Terms out;
    for (const auto& [e,c] : terms_) {
        if (e[static_cast<size_t>(axis)]==0) continue;
        Exponents d=e; const int power=d[static_cast<size_t>(axis)]--;
        out[d] += c*static_cast<double>(power);
    }
    return Polynomial(std::move(out));
}
int Polynomial::degree() const {
    int d=-1;
    for (const auto& [e,c] : terms_) { (void)c; d=std::max(d,e[0]+e[1]+e[2]); }
    return d;
}
bool Polynomial::isZero(double epsilon) const {
    for (const auto& [e,c] : terms_) { (void)e; if (std::abs(c)>=epsilon) return false; }
    return true;
}
Polynomial Polynomial::operator+(const Polynomial& rhs) const {
    Terms out=terms_; for (const auto& [e,c]:rhs.terms_) out[e]+=c; return Polynomial(std::move(out));
}
Polynomial Polynomial::operator-(const Polynomial& rhs) const { return *this + rhs*(-1.0); }
Polynomial Polynomial::operator*(double scalar) const {
    Terms out=terms_; for (auto& [e,c]:out) { (void)e; c*=scalar; } return Polynomial(std::move(out));
}
Polynomial Polynomial::operator*(const Polynomial& rhs) const {
    Terms out;
    for (const auto& [a,ca]:terms_) for (const auto& [b,cb]:rhs.terms_) {
        Exponents e{a[0]+b[0],a[1]+b[1],a[2]+b[2]};
        if (e[0]+e[1]+e[2]>12) throw std::runtime_error("Polynomial degree would exceed the supported maximum (12)");
        out[e]+=ca*cb;
        if (out.size()>20000) throw std::runtime_error("Expression expands to too many polynomial terms");
    }
    return Polynomial(std::move(out));
}
Polynomial Polynomial::pow(unsigned exponent) const {
    if (exponent>12) throw std::runtime_error("Powers above 12 are not supported");
    Polynomial result(1.0), base=*this; unsigned n=exponent;
    while (n>0) { if (n&1U) result=result*base; n>>=1U; if (n) base=base*base; }
    return result;
}

class PolynomialParser {
public:
    explicit PolynomialParser(std::string_view text) : text_(text) {}
    Polynomial parseAll() {
        Polynomial value=parseExpression(); skipSpace();
        if (pos_!=text_.size()) fail("Unexpected character");
        return value;
    }
private:
    std::string_view text_; size_t pos_{0};
    [[noreturn]] void fail(const std::string& what) const {
        throw std::runtime_error(what+" at character "+std::to_string(pos_+1));
    }
    void skipSpace() { while (pos_<text_.size() && std::isspace(static_cast<unsigned char>(text_[pos_]))) ++pos_; }
    bool consume(char c) { skipSpace(); if (pos_<text_.size() && text_[pos_]==c) { ++pos_; return true; } return false; }
    bool nextStartsImplicitFactor() {
        skipSpace(); if (pos_>=text_.size()) return false;
        const char c=text_[pos_]; return c=='x'||c=='y'||c=='z'||c=='(';
    }
    Polynomial parseExpression() {
        Polynomial value=parseTerm();
        while (true) {
            if (consume('+')) value=value+parseTerm();
            else if (consume('-')) value=value-parseTerm();
            else break;
        }
        return value;
    }
    Polynomial parseTerm() {
        Polynomial value=parseUnary();
        while (true) {
            if (consume('*')) value=value*parseUnary();
            else if (nextStartsImplicitFactor()) value=value*parseUnary();
            else break;
        }
        return value;
    }
    Polynomial parseUnary() {
        if (consume('+')) return parseUnary();
        if (consume('-')) return parseUnary()*(-1.0);
        return parsePower();
    }
    Polynomial parsePower() {
        Polynomial base=parsePrimary();
        if (consume('^')) {
            skipSpace(); const size_t start=pos_;
            while (pos_<text_.size() && std::isdigit(static_cast<unsigned char>(text_[pos_]))) ++pos_;
            if (start==pos_) fail("Exponent must be a non-negative integer");
            unsigned n=0;
            try { n=static_cast<unsigned>(std::stoul(std::string(text_.substr(start,pos_-start)))); }
            catch (...) { fail("Invalid exponent"); }
            if (n>12) fail("Exponent above 12 is not supported");
            base=base.pow(n);
        }
        return base;
    }
    Polynomial parsePrimary() {
        skipSpace(); if (pos_>=text_.size()) fail("Expected a number, variable, or parenthesis");
        const char c=text_[pos_];
        if (c=='(') {
            ++pos_; Polynomial value=parseExpression();
            if (!consume(')')) fail("Missing closing parenthesis");
            return value;
        }
        if (c=='x'||c=='y'||c=='z') { ++pos_; return Polynomial::variable(c=='x'?0:(c=='y'?1:2)); }
        if (std::isdigit(static_cast<unsigned char>(c)) || c=='.') {
            const std::string remainder(text_.substr(pos_)); char* end=nullptr;
            const double v=std::strtod(remainder.c_str(),&end);
            if (end==remainder.c_str() || !std::isfinite(v)) fail("Invalid number");
            pos_+=static_cast<size_t>(end-remainder.c_str()); return Polynomial(v);
        }
        fail("Expected a number, x, y, z, or parenthesis");
    }
};

Polynomial Polynomial::parse(const std::string& expression) {
    if (expression.empty()) throw std::runtime_error("Enter a polynomial expression");
    // Accept either F=0 or an implicit equation F_left=F_right (exactly one '=').
    size_t eq=expression.find('=');
    if (eq!=std::string::npos) {
        if (expression.find('=',eq+1)!=std::string::npos) throw std::runtime_error("Only one '=' is allowed in an equation");
        Polynomial left=PolynomialParser(std::string_view(expression).substr(0,eq)).parseAll();
        Polynomial right=PolynomialParser(std::string_view(expression).substr(eq+1)).parseAll();
        return left-right;
    }
    return PolynomialParser(expression).parseAll();
}

PolynomialModel::PolynomialModel(Polynomial polynomial)
    : f(std::move(polynomial)), fx(f.derivative(0)), fy(f.derivative(1)), fz(f.derivative(2)) {
    const std::array<Polynomial,3> first{fx,fy,fz};
    for (size_t i=0;i<3;++i) for (size_t j=0;j<3;++j) hessian[i][j]=first[i].derivative(static_cast<int>(j));
}
std::array<std::array<double,3>,3> PolynomialModel::hessianAt(const Vec3& p) const {
    std::array<std::array<double,3>,3> out{};
    for (size_t i=0;i<3;++i) for (size_t j=0;j<3;++j) out[i][j]=hessian[i][j].evaluate(p);
    return out;
}
} // namespace isa
