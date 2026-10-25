#pragma once

#include "isa/Mesh.h"
#include <pmp/mat_vec.h>
#include <pmp/viewers/drawable.h>
#include <pmp/viewers/shader.h>
#include <pmp/viewers/gl.h>

#include <cstdint>
#include <vector>

namespace isa {

// OpenGL line overlay used to render the sampling lattice and an original-mesh
// wireframe without adding those lines/faces to PMP's mesh statistics.
class LineOverlay final : public pmp::Drawable {
public:
    explicit LineOverlay(const pmp::vec3& color, float width = 1.0f);
    ~LineOverlay() override;

    void setEnabled(bool enabled) { enabled_ = enabled; }
    bool enabled() const { return enabled_; }
    void setGrid(int cells, double extent);
    void setMeshEdges(const Mesh& mesh);
    void clear();

    void update_buffers() override;
    void draw(const pmp::mat4& projection, const pmp::mat4& modelview) override;

private:
    void appendLine(double ax, double ay, double az,
                    double bx, double by, double bz);

    std::vector<float> positions_;
    pmp::vec3 color_;
    float lineWidth_;
    bool enabled_{false};
    GLuint vao_{0};
    GLuint vbo_{0};
    pmp::Shader shader_;
};

} // namespace isa
