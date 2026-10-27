#pragma once

#include "isa/Mesh.h"
#include "isa/LineOverlay.h"
#include "isa/Singularity.h"
#include <pmp/viewers/mesh_viewer.h>

#include <memory>
#include <string>
#include <vector>

namespace isa {

struct KSweepRow {
    double k{0.0};
    std::size_t vertices{0};
    std::size_t triangles{0};
    std::size_t components{0};
    double normalizedRmsDistance{0.0};
};

// A desktop application built on PMP's GLFW/OpenGL viewer and Dear ImGui.
// It deliberately has no Qt dependency.
class SingularityViewer final : public pmp::MeshViewer {
public:
    SingularityViewer();

protected:
    void process_imgui() override;

private:
    char functionText_[512] = "x^2+y^2-z^2";
    char perturbationText_[256] = "-x";
    char pmpExportPath_[512] = "isosurface_pmp.obj";
    char rawExportPath_[512] = "isosurface_raw.obj";
    char csvExportPath_[512] = "k_sweep.csv";

    int presetIndex_{0};
    int cells_{32};
    int sampleCount_{9};
    double extent_{1.5};
    double convergenceTolerance_{0.002};
    double k_{0.02};
    double kMin_{-0.02};
    double kMax_{0.02};
    float patchAlpha_{0.5f};
    std::string kSweepFunction_;
    std::string kSweepG_;
    int kSweepCells_{0};
    double kSweepExtent_{0.0};
    bool usePerturbation_{false};
    bool automaticResolution_{false};
    bool showGrid_{false};
    bool showOriginalWireframe_{true};
    bool meshWasExtractedWithPerturbation_{false};

    std::shared_ptr<LineOverlay> gridOverlay_;
    std::shared_ptr<LineOverlay> originalWireOverlay_;
    std::string viewerRepresentation_ = "not built";
    std::size_t rejectedDisplayFaces_{0};
    Mesh currentMesh_;
    Polynomial currentActiveField_;
    std::string currentFieldDescription_ = "No extracted field yet";
    std::vector<SingularityCandidate> singularities_;
    std::vector<KSweepRow> kRows_;
    std::string status_ = "Ready. Select a preset and extract the zero isosurface.";
    std::string metricsText_ = "The current mesh metrics will appear here.";

    Polynomial basePolynomialFromUi() const;
    Polynomial activePolynomialFromUi(double overrideK = 0.0, bool useOverrideK = false) const;
    GridSpec gridFromUi(int cells = -1) const;
    void setPreset(int index);
    void extractCurrent();
    void selectAutomaticResolution();
    void detectSingularities();
    void runKSweep();
    void createInterpolatedPatchFromSweep();
    void exportCurrentPmpMesh();
    void exportCurrentRawObj();
    void exportSweepCsv();
    void rebuildViewerMesh();
    void updateMetrics(const Polynomial& activeField, const std::string& note = {});
    void setStatus(std::string message);
};

} // namespace isa
