#include "isa/LineOverlay.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <utility>
#include <iostream>
#include <iomanip>
#include <limits>

namespace isa {
namespace {
constexpr const char* kLineVertexShader = R"GLSL(#version 330 core
layout(location = 0) in vec3 position;
uniform mat4 mvp;
void main() { gl_Position = mvp * vec4(position, 1.0); }
)GLSL";
constexpr const char* kLineFragmentShader = R"GLSL(#version 330 core
uniform vec3 line_color;
out vec4 frag_color;
void main() { frag_color = vec4(line_color, 1.0); }
)GLSL";

const char* glErrorName(GLenum error) {
    switch (error) {
        case GL_INVALID_ENUM: return "GL_INVALID_ENUM";
        case GL_INVALID_VALUE: return "GL_INVALID_VALUE";
        case GL_INVALID_OPERATION: return "GL_INVALID_OPERATION";
        case GL_INVALID_FRAMEBUFFER_OPERATION: return "GL_INVALID_FRAMEBUFFER_OPERATION";
        case GL_OUT_OF_MEMORY: return "GL_OUT_OF_MEMORY";
        default: return "unknown OpenGL error";
    }
}

// PMP checks the GL error queue immediately before it draws its overlays. An
// overlay error can otherwise be reported on the following frame, far from its
// source. Drain and label errors here while the overlay is the active operation.
void reportOverlayGlErrors(const char* stage) {
    bool found = false;
    for (GLenum error = glGetError(); error != GL_NO_ERROR; error = glGetError()) {
        if (!found) {
            std::cerr << "[LineOverlay] OpenGL error after " << stage << ": ";
            found = true;
        } else {
            std::cerr << "[LineOverlay] additional OpenGL error after " << stage << ": ";
        }
        std::cerr << glErrorName(error) << " (0x" << std::hex << error << std::dec << ")\n";
    }
}
}

LineOverlay::LineOverlay(const pmp::vec3& color, float width)
    : color_(color), lineWidth_(std::max(1.0f, width)) {}

LineOverlay::~LineOverlay() {
    if (vbo_ != 0) glDeleteBuffers(1, &vbo_);
    if (vao_ != 0) glDeleteVertexArrays(1, &vao_);
}

void LineOverlay::appendLine(double ax, double ay, double az,
                             double bx, double by, double bz) {
    positions_.push_back(static_cast<float>(ax));
    positions_.push_back(static_cast<float>(ay));
    positions_.push_back(static_cast<float>(az));
    positions_.push_back(static_cast<float>(bx));
    positions_.push_back(static_cast<float>(by));
    positions_.push_back(static_cast<float>(bz));
}

void LineOverlay::setGrid(int cells, double extent) {
    positions_.clear();
    if (cells < 1 || !std::isfinite(extent) || extent <= 0.0) return;
    const double step = (2.0 * extent) / static_cast<double>(cells);

    // Each full-length line corresponds to one family of parallel grid edges.
    // This yields the same axis-aligned lattice as the Marching Cubes sampler.
    for (int j = 0; j <= cells; ++j) {
        const double y = -extent + j * step;
        for (int k = 0; k <= cells; ++k) {
            const double z = -extent + k * step;
            appendLine(-extent, y, z, extent, y, z);
        }
    }
    for (int i = 0; i <= cells; ++i) {
        const double x = -extent + i * step;
        for (int k = 0; k <= cells; ++k) {
            const double z = -extent + k * step;
            appendLine(x, -extent, z, x, extent, z);
        }
    }
    for (int i = 0; i <= cells; ++i) {
        const double x = -extent + i * step;
        for (int j = 0; j <= cells; ++j) {
            const double y = -extent + j * step;
            appendLine(x, y, -extent, x, y, extent);
        }
    }
}

void LineOverlay::setMeshEdges(const Mesh& mesh) {
    positions_.clear();
    using Edge = std::pair<std::uint32_t, std::uint32_t>;
    std::set<Edge> edges;
    for (const auto& triangle : mesh.triangles) {
        for (int i = 0; i < 3; ++i) {
            auto a = triangle[static_cast<std::size_t>(i)];
            auto b = triangle[static_cast<std::size_t>((i + 1) % 3)];
            if (a >= mesh.vertices.size() || b >= mesh.vertices.size() || a == b) continue;
            if (b < a) std::swap(a, b);
            edges.emplace(a, b);
        }
    }
    for (const auto& edge : edges) {
        const auto& a = mesh.vertices[edge.first];
        const auto& b = mesh.vertices[edge.second];
        appendLine(a.x, a.y, a.z, b.x, b.y, b.z);
    }
}

void LineOverlay::clear() {
    positions_.clear();
}

void LineOverlay::update_buffers() {
    if (vao_ == 0) {
        glGenVertexArrays(1, &vao_);
        glGenBuffers(1, &vbo_);
    }
    glBindVertexArray(vao_);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_);
    glBufferData(GL_ARRAY_BUFFER,
                 static_cast<GLsizeiptr>(positions_.size() * sizeof(float)),
                 positions_.empty() ? nullptr : positions_.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);
}

void LineOverlay::draw(const pmp::mat4& projection, const pmp::mat4& modelview) {
    if (!enabled_ || positions_.empty()) return;
    if (!shader_.is_valid()) {
        shader_.source(kLineVertexShader, kLineFragmentShader);
        reportOverlayGlErrors("line shader creation");
    }

    GLint previousDepthFunc = GL_LESS;
    GLfloat previousLineWidth = 1.0f;
    const GLboolean depthWasEnabled = glIsEnabled(GL_DEPTH_TEST);
    glGetIntegerv(GL_DEPTH_FUNC, &previousDepthFunc);
    glGetFloatv(GL_LINE_WIDTH, &previousLineWidth);

    // Wide OpenGL lines are not portable. In particular, macOS core-profile
    // contexts commonly support only width 1.0. Requesting 1.8 can generate
    // GL_INVALID_VALUE on some drivers; clamp to the driver's actual range.
    GLfloat supportedLineWidths[2] = {1.0f, 1.0f};
    glGetFloatv(GL_ALIASED_LINE_WIDTH_RANGE, supportedLineWidths);
    reportOverlayGlErrors("querying supported line-width range");

    GLfloat safeLineWidth = 1.0f;
    GLfloat safePreviousLineWidth = 1.0f;
    if (std::isfinite(supportedLineWidths[0]) &&
        std::isfinite(supportedLineWidths[1]) &&
        supportedLineWidths[0] > 0.0f &&
        supportedLineWidths[1] >= supportedLineWidths[0]) {
        safeLineWidth = std::clamp(lineWidth_, supportedLineWidths[0], supportedLineWidths[1]);
        safePreviousLineWidth = std::clamp(previousLineWidth, supportedLineWidths[0], supportedLineWidths[1]);
    }

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glLineWidth(safeLineWidth);
    reportOverlayGlErrors("setting depth state and supported line width");

    shader_.use();
    shader_.set_uniform("mvp", projection * modelview);
    shader_.set_uniform("line_color", color_);
    reportOverlayGlErrors("line shader uniforms");

    glBindVertexArray(vao_);
    const std::size_t vertexCount = positions_.size() / 3;
    if (vertexCount <= static_cast<std::size_t>(std::numeric_limits<GLsizei>::max())) {
        glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(vertexCount));
    } else {
        std::cerr << "[LineOverlay] Skipping draw: vertex count exceeds GLsizei range ("
                  << vertexCount << ").\n";
    }
    glBindVertexArray(0);
    shader_.disable();
    reportOverlayGlErrors("line draw call");

    glLineWidth(safePreviousLineWidth);
    glDepthFunc(static_cast<GLenum>(previousDepthFunc));
    if (!depthWasEnabled) glDisable(GL_DEPTH_TEST);
    reportOverlayGlErrors("OpenGL state restoration");
}

} // namespace isa
