#pragma once
#include <cmath>

namespace isa {
struct Vec3 {
    double x{0.0}, y{0.0}, z{0.0};
    Vec3() = default;
    Vec3(double px, double py, double pz) : x(px), y(py), z(pz) {}
    Vec3 operator+(const Vec3& b) const { return {x+b.x, y+b.y, z+b.z}; }
    Vec3 operator-(const Vec3& b) const { return {x-b.x, y-b.y, z-b.z}; }
    Vec3 operator*(double s) const { return {x*s, y*s, z*s}; }
    Vec3 operator/(double s) const { return {x/s, y/s, z/s}; }
    Vec3& operator+=(const Vec3& b) { x+=b.x; y+=b.y; z+=b.z; return *this; }
};
inline double dot(const Vec3& a, const Vec3& b) { return a.x*b.x+a.y*b.y+a.z*b.z; }
inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};
}
inline double norm2(const Vec3& v) { return dot(v,v); }
inline double norm(const Vec3& v) { return std::sqrt(norm2(v)); }
inline Vec3 normalized(const Vec3& v) { const double n=norm(v); return n>1e-15 ? v/n : Vec3{}; }
} // namespace isa
