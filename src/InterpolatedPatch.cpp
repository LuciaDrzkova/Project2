#include "isa/InterpolatedPatch.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace isa {
namespace {

struct Bounds {
    Vec3 lo{ std::numeric_limits<double>::infinity(),
             std::numeric_limits<double>::infinity(),
             std::numeric_limits<double>::infinity() };
    Vec3 hi{-std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity() };
};

void expand(Bounds& b, const Vec3& p) {
    b.lo.x = std::min(b.lo.x, p.x); b.hi.x = std::max(b.hi.x, p.x);
    b.lo.y = std::min(b.lo.y, p.y); b.hi.y = std::max(b.hi.y, p.y);
    b.lo.z = std::min(b.lo.z, p.z); b.hi.z = std::max(b.hi.z, p.z);
}

Vec3 componentMin(const Vec3& a, const Vec3& b) {
    return {std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z)};
}
Vec3 componentMax(const Vec3& a, const Vec3& b) {
    return {std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)};
}
double component(const Vec3& p, int axis) {
    return axis == 0 ? p.x : (axis == 1 ? p.y : p.z);
}
double distanceToBoundsSquared(const Vec3& p, const Bounds& b) {
    const double dx = p.x < b.lo.x ? b.lo.x - p.x : (p.x > b.hi.x ? p.x - b.hi.x : 0.0);
    const double dy = p.y < b.lo.y ? b.lo.y - p.y : (p.y > b.hi.y ? p.y - b.hi.y : 0.0);
    const double dz = p.z < b.lo.z ? b.lo.z - p.z : (p.z > b.hi.z ? p.z - b.hi.z : 0.0);
    return dx * dx + dy * dy + dz * dz;
}

Vec3 closestPointOnSegment(const Vec3& p, const Vec3& a, const Vec3& b) {
    const Vec3 ab = b - a;
    const double denom = norm2(ab);
    if (denom <= 1e-30) return a;
    const double t = std::clamp(dot(p - a, ab) / denom, 0.0, 1.0);
    return a + ab * t;
}

// Closest point on triangle using Voronoi regions. Degenerate triangles fall
// back to their three segments so the query remains well-defined.
Vec3 closestPointOnTriangle(const Vec3& p, const Vec3& a, const Vec3& b, const Vec3& c) {
    const Vec3 ab = b - a;
    const Vec3 ac = c - a;
    if (norm2(cross(ab, ac)) <= 1e-28) {
        const Vec3 q0 = closestPointOnSegment(p, a, b);
        const Vec3 q1 = closestPointOnSegment(p, b, c);
        const Vec3 q2 = closestPointOnSegment(p, c, a);
        const double d0 = norm2(p - q0), d1 = norm2(p - q1), d2 = norm2(p - q2);
        return d0 <= d1 && d0 <= d2 ? q0 : (d1 <= d2 ? q1 : q2);
    }

    const Vec3 ap = p - a;
    const double d1 = dot(ab, ap), d2 = dot(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0) return a;

    const Vec3 bp = p - b;
    const double d3 = dot(ab, bp), d4 = dot(ac, bp);
    if (d3 >= 0.0 && d4 <= d3) return b;

    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        const double v = d1 / (d1 - d3);
        return a + ab * v;
    }

    const Vec3 cp = p - c;
    const double d5 = dot(ab, cp), d6 = dot(ac, cp);
    if (d6 >= 0.0 && d5 <= d6) return c;

    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        const double w = d2 / (d2 - d6);
        return a + ac * w;
    }

    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return b + (c - b) * w;
    }

    const double denom = va + vb + vc;
    if (std::abs(denom) <= 1e-30) {
        const Vec3 q0 = closestPointOnSegment(p, a, b);
        const Vec3 q1 = closestPointOnSegment(p, b, c);
        const Vec3 q2 = closestPointOnSegment(p, c, a);
        const double e0 = norm2(p - q0), e1 = norm2(p - q1), e2 = norm2(p - q2);
        return e0 <= e1 && e0 <= e2 ? q0 : (e1 <= e2 ? q1 : q2);
    }
    const double inv = 1.0 / denom;
    const double v = vb * inv;
    const double w = vc * inv;
    return a + ab * v + ac * w;
}

class TriangleBvh {
public:
    explicit TriangleBvh(const Mesh& mesh) : mesh_(mesh) {
        records_.reserve(mesh.triangles.size());
        for (const auto& ids : mesh.triangles) {
            if (ids[0] >= mesh.vertices.size() || ids[1] >= mesh.vertices.size() || ids[2] >= mesh.vertices.size())
                continue;
            const Vec3& a = mesh.vertices[ids[0]];
            const Vec3& b = mesh.vertices[ids[1]];
            const Vec3& c = mesh.vertices[ids[2]];
            TriangleRecord rec;
            rec.ids = ids;
            rec.bounds.lo = componentMin(a, componentMin(b, c));
            rec.bounds.hi = componentMax(a, componentMax(b, c));
            rec.centroid = (a + b + c) / 3.0;
            records_.push_back(rec);
        }
        if (records_.empty()) throw std::invalid_argument("Cannot interpolate against an empty/invalid target mesh.");
        order_.resize(records_.size());
        std::iota(order_.begin(), order_.end(), std::size_t{0});
        nodes_.reserve(records_.size() * 2);
        build(0, order_.size());
    }

    Vec3 closestPoint(const Vec3& query) const {
        double bestDistance2 = std::numeric_limits<double>::infinity();
        Vec3 bestPoint{};
        queryNode(0, query, bestDistance2, bestPoint);
        if (!std::isfinite(bestDistance2) || !std::isfinite(bestPoint.x) ||
            !std::isfinite(bestPoint.y) || !std::isfinite(bestPoint.z)) {
            throw std::runtime_error("Closest-point projection failed while constructing the interpolated patch.");
        }
        return bestPoint;
    }

private:
    struct TriangleRecord {
        std::array<std::uint32_t, 3> ids{};
        Bounds bounds;
        Vec3 centroid;
    };
    struct Node {
        Bounds bounds;
        std::size_t begin{0}, end{0};
        int left{-1}, right{-1};
        bool leaf() const { return left < 0 && right < 0; }
    };

    const Mesh& mesh_;
    std::vector<TriangleRecord> records_;
    std::vector<std::size_t> order_;
    std::vector<Node> nodes_;

    int build(std::size_t begin, std::size_t end) {
        Node node;
        node.begin = begin;
        node.end = end;
        Bounds centroidBounds;
        for (std::size_t i = begin; i < end; ++i) {
            const TriangleRecord& tri = records_[order_[i]];
            node.bounds.lo = componentMin(node.bounds.lo, tri.bounds.lo);
            node.bounds.hi = componentMax(node.bounds.hi, tri.bounds.hi);
            expand(centroidBounds, tri.centroid);
        }
        const int index = static_cast<int>(nodes_.size());
        nodes_.push_back(node);
        if (end - begin <= 8) return index;

        const Vec3 span = centroidBounds.hi - centroidBounds.lo;
        int axis = span.x >= span.y && span.x >= span.z ? 0 : (span.y >= span.z ? 1 : 2);
        if (component(span, axis) <= 1e-15) return index;
        const std::size_t mid = begin + (end - begin) / 2;
        std::nth_element(order_.begin() + static_cast<std::ptrdiff_t>(begin),
                         order_.begin() + static_cast<std::ptrdiff_t>(mid),
                         order_.begin() + static_cast<std::ptrdiff_t>(end),
                         [&](std::size_t lhs, std::size_t rhs) {
                             return component(records_[lhs].centroid, axis) < component(records_[rhs].centroid, axis);
                         });
        const int left = build(begin, mid);
        const int right = build(mid, end);
        nodes_[static_cast<std::size_t>(index)].left = left;
        nodes_[static_cast<std::size_t>(index)].right = right;
        return index;
    }

    void queryNode(int nodeIndex, const Vec3& query,
                   double& bestDistance2, Vec3& bestPoint) const {
        const Node& node = nodes_[static_cast<std::size_t>(nodeIndex)];
        if (distanceToBoundsSquared(query, node.bounds) > bestDistance2) return;
        if (node.leaf()) {
            for (std::size_t i = node.begin; i < node.end; ++i) {
                const TriangleRecord& tri = records_[order_[i]];
                const Vec3& a = mesh_.vertices[tri.ids[0]];
                const Vec3& b = mesh_.vertices[tri.ids[1]];
                const Vec3& c = mesh_.vertices[tri.ids[2]];
                const Vec3 candidate = closestPointOnTriangle(query, a, b, c);
                const double distance2 = norm2(query - candidate);
                if (distance2 < bestDistance2) {
                    bestDistance2 = distance2;
                    bestPoint = candidate;
                }
            }
            return;
        }

        const Node& left = nodes_[static_cast<std::size_t>(node.left)];
        const Node& right = nodes_[static_cast<std::size_t>(node.right)];
        const bool leftFirst = distanceToBoundsSquared(query, left.bounds) <= distanceToBoundsSquared(query, right.bounds);
        const int first = leftFirst ? node.left : node.right;
        const int second = leftFirst ? node.right : node.left;
        queryNode(first, query, bestDistance2, bestPoint);
        queryNode(second, query, bestDistance2, bestPoint);
    }
};

} // namespace

Mesh makeLinearInterpolatedPatch(const Mesh& reference,
                                 const Mesh& negativeK,
                                 const Mesh& positiveK,
                                 double alpha) {
    if (reference.empty()) throw std::invalid_argument("The original F=0 reference mesh is empty.");
    if (negativeK.empty()) throw std::invalid_argument("The negative-k isosurface is empty.");
    if (positiveK.empty()) throw std::invalid_argument("The positive-k isosurface is empty.");
    if (!std::isfinite(alpha)) throw std::invalid_argument("Patch interpolation alpha must be finite.");
    alpha = std::clamp(alpha, 0.0, 1.0);

    // Endpoint invariants: alpha=0/1 must be the exact meshes selected by the
    // k sweep, not a projection of the original F mesh onto those surfaces.
    // The projection-based common-connectivity approximation is used only for
    // intermediate alpha values, where the input meshes may have different
    // triangulations or topology.
    if (alpha == 0.0) return negativeK;
    if (alpha == 1.0) return positiveK;

    const TriangleBvh negativeBvh(negativeK);
    const TriangleBvh positiveBvh(positiveK);

    Mesh patch;
    patch.vertices.reserve(reference.vertices.size());
    patch.triangles = reference.triangles;
    for (const Vec3& point : reference.vertices) {
        const Vec3 onNegative = negativeBvh.closestPoint(point);
        const Vec3 onPositive = positiveBvh.closestPoint(point);
        patch.vertices.push_back(onNegative * (1.0 - alpha) + onPositive * alpha);
    }

    // Preserve the reference mesh's local orientation after the vertex motion.
    // Closest-point interpolation can occasionally move a face through a folded
    // configuration; reverse such a face if its normal opposes the reference.
    for (auto& triangle : patch.triangles) {
        const auto a = triangle[0], b = triangle[1], c = triangle[2];
        if (a >= reference.vertices.size() || b >= reference.vertices.size() || c >= reference.vertices.size() ||
            a >= patch.vertices.size() || b >= patch.vertices.size() || c >= patch.vertices.size()) continue;
        const Vec3 referenceNormal = cross(reference.vertices[b] - reference.vertices[a],
                                           reference.vertices[c] - reference.vertices[a]);
        const Vec3 patchNormal = cross(patch.vertices[b] - patch.vertices[a],
                                       patch.vertices[c] - patch.vertices[a]);
        if (norm2(referenceNormal) > 1e-28 && norm2(patchNormal) > 1e-28 &&
            dot(referenceNormal, patchNormal) < 0.0) {
            std::swap(triangle[1], triangle[2]);
        }
    }

    // This patch is parameterized by the original mesh so every vertex and face
    // has a stable correspondence even if the +/-k Marching Cubes meshes differ.
    // The result is a linear geometric interpolant, not a topology guarantee.
    return patch;
}

} // namespace isa
