#include "isa/SingularityViewer.h"

#include "isa/MarchingCubes.h"
#include "isa/Metrics.h"
#include "isa/PmpExport.h"
#include "isa/Polynomial.h"

#include <imgui.h>
#include <pmp/mat_vec.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace isa {
namespace {

constexpr const char* kPresetNames[] = {
    "Double cone (singular)",
    "Sphere (smooth)",
    "One-sheet hyperboloid",
    "Crossing planes",
    "Cusp-like surface",
    "Ellipsoid"
};
constexpr const char* kPresetExpressions[] = {
    "x^2+y^2-z^2",
    "x^2+y^2+z^2-1",
    "x^2+y^2-z^2-0.2*x",
    "x*y",
    "z^2-x^3-y^2",
    "x^2+2*y^2+0.5*z^2-1"
};
constexpr int kNumPresets = static_cast<int>(sizeof(kPresetNames) / sizeof(kPresetNames[0]));

std::string formatNumber(double value, int precision = 6) {
    if (!std::isfinite(value)) return "n/a";
    std::ostringstream out;
    out << std::setprecision(precision) << std::defaultfloat << value;
    return out.str();
}

std::string activeDescription(const char* f, bool perturbed, const char* g, double k) {
    if (!perturbed) return std::string(f);
    return "(" + std::string(f) + ") + (" + formatNumber(k) + ") * (" + std::string(g) + ")";
}

void copyToBuffer(char* destination, std::size_t capacity, const char* text) {
    if (capacity == 0) return;
    std::size_t i = 0;
    if (text) {
        for (; i + 1 < capacity && text[i] != '\0'; ++i) destination[i] = text[i];
    }
    destination[i] = '\0';
}

} // namespace

SingularityViewer::SingularityViewer()
    : pmp::MeshViewer("Singularity-aware Isosurface Explorer", 1480, 960, true) {
    gridOverlay_ = std::make_shared<LineOverlay>(pmp::vec3(0.22, 0.38, 0.48), 1.0f);
    originalWireOverlay_ = std::make_shared<LineOverlay>(pmp::vec3(1.0, 0.34, 0.08), 1.8f);
    gridOverlay_->setEnabled(showGrid_);
    originalWireOverlay_->setEnabled(false);
    renderer_.add_overlay(gridOverlay_);
    renderer_.add_overlay(originalWireOverlay_);
    set_scene(pmp::vec3(0.0, 0.0, 0.0), static_cast<float>(extent_ * 1.8));
    setStatus("Ready. The GUI is provided by PMP's viewer and Dear ImGui; Qt is not used.");
}

Polynomial SingularityViewer::basePolynomialFromUi() const {
    return Polynomial::parse(std::string(functionText_));
}

Polynomial SingularityViewer::activePolynomialFromUi(double overrideK, bool useOverrideK) const {
    Polynomial base = basePolynomialFromUi();
    if (!usePerturbation_) return base;
    const Polynomial g = Polynomial::parse(std::string(perturbationText_));
    return base + g * (useOverrideK ? overrideK : k_);
}

GridSpec SingularityViewer::gridFromUi(int cells) const {
    return {cells > 0 ? cells : cells_, extent_};
}

void SingularityViewer::setPreset(int index) {
    if (index < 0 || index >= kNumPresets) return;
    presetIndex_ = index;
    copyToBuffer(functionText_, sizeof(functionText_), kPresetExpressions[index]);
}

void SingularityViewer::setStatus(std::string message) {
    status_ = std::move(message);
}

void SingularityViewer::process_imgui() {
    pmp::MeshViewer::process_imgui();

    if (ImGui::Begin("Singularity-aware Isosurface Explorer")) {
        ImGui::TextWrapped("Polynomial zero-set extraction with singularity diagnostics. The original F is retained separately from each F + kG experiment.");
        ImGui::Separator();

        if (ImGui::BeginCombo("Preset", kPresetNames[presetIndex_])) {
            for (int i = 0; i < kNumPresets; ++i) {
                const bool selected = (presetIndex_ == i);
                if (ImGui::Selectable(kPresetNames[i], selected)) setPreset(i);
                if (selected) ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
        }
        ImGui::InputText("F(x,y,z) = 0", functionText_, sizeof(functionText_));
        ImGui::TextWrapped("Polynomial syntax: +, -, *, parentheses, implicit multiplication (2x), and nonnegative integer powers. Use expressions such as x^2+y^2-z^2.");
        ImGui::InputDouble("Domain half-width", &extent_, 0.1, 0.25, "%.4f");
        extent_ = std::clamp(extent_, 0.1, 100.0);

        ImGui::Separator();
        ImGui::TextUnformatted("Perturbation experiment");
        ImGui::Checkbox("Extract F + kG", &usePerturbation_);
        ImGui::InputText("G(x,y,z)", perturbationText_, sizeof(perturbationText_));
        ImGui::InputDouble("k", &k_, 0.001, 0.01, "%.8g");
        ImGui::TextWrapped("A nonzero k can remove a singularity or change topology. A small k may produce a feature smaller than one grid cell.");
        ImGui::TextWrapped("The displayed surface changes only after you click Extract / update mesh.");
        if (std::string(functionText_) == "x^2+y^2-z^2" && std::string(perturbationText_) == "-x") {
            const double cellWidth = (2.0 * extent_) / static_cast<double>(std::max(1, cells_));
            ImGui::TextWrapped("For this cone with G=-x, the perturbed waist radius is |k|/2 = %.6g; current grid step is h = %.6g. If the radius is much smaller than h, Marching Cubes cannot reliably resolve the visible difference.", std::abs(k_) * 0.5, cellWidth);
        }
        if (ImGui::Checkbox("Overlay original F as orange wireframe", &showOriginalWireframe_)) {
            originalWireOverlay_->setEnabled(showOriginalWireframe_ && meshWasExtractedWithPerturbation_);
        }

        ImGui::Separator();
        ImGui::TextUnformatted("Extraction and resolution");
        ImGui::Checkbox("Automatic resolution (geometry + topology heuristic)", &automaticResolution_);
        ImGui::InputInt("Grid cells per axis", &cells_);
        cells_ = std::clamp(cells_, 8, 180);
        ImGui::InputDouble("Normalized convergence tolerance", &convergenceTolerance_, 0.0005, 0.001, "%.6f");
        convergenceTolerance_ = std::clamp(convergenceTolerance_, 0.00001, 0.1);
        if (ImGui::Checkbox("Show 3D sampling grid", &showGrid_)) gridOverlay_->setEnabled(showGrid_);
        ImGui::TextWrapped("Grid lines show the uniform cubic sampling lattice used by Marching Cubes; they update on the next extraction.");
        if (ImGui::Button("Extract / update mesh")) extractCurrent();
        ImGui::SameLine();
        if (ImGui::Button("Find singularities on original F")) detectSingularities();

        ImGui::Separator();
        ImGui::TextUnformatted("Automatic k sweep");
        ImGui::InputDouble("k minimum", &kMin_, 0.005, 0.01, "%.8g");
        ImGui::InputDouble("k maximum", &kMax_, 0.005, 0.01, "%.8g");
        ImGui::InputInt("k sample count", &sampleCount_);
        sampleCount_ = std::clamp(sampleCount_, 3, 41);
        if (ImGui::Button("Run k sweep")) runKSweep();

        ImGui::Separator();
        ImGui::TextUnformatted("Exports (paths are relative to the app's working directory unless absolute)");
        ImGui::InputText("PMP mesh path", pmpExportPath_, sizeof(pmpExportPath_));
        if (ImGui::Button("Export with PMP")) exportCurrentPmpMesh();
        ImGui::InputText("Raw OBJ path", rawExportPath_, sizeof(rawExportPath_));
        if (ImGui::Button("Export raw OBJ (all triangles)")) exportCurrentRawObj();
        ImGui::InputText("Sweep CSV path", csvExportPath_, sizeof(csvExportPath_));
        if (ImGui::Button("Export k sweep CSV")) exportSweepCsv();

        ImGui::Separator();
        ImGui::TextUnformatted("Status");
        ImGui::TextWrapped("%s", status_.c_str());
        ImGui::Separator();
        ImGui::TextUnformatted("Mesh diagnostics");
        ImGui::TextWrapped("%s", metricsText_.c_str());

        if (!singularities_.empty()) {
            ImGui::Separator();
            ImGui::Text("Singularity candidates: %llu", static_cast<unsigned long long>(singularities_.size()));
            for (std::size_t i = 0; i < singularities_.size(); ++i) {
                const auto& s = singularities_[i];
                ImGui::Text("%llu: (%.7g, %.7g, %.7g), |F|=%.3g, |grad F|=%.3g",
                    static_cast<unsigned long long>(i + 1), s.position.x, s.position.y, s.position.z,
                    s.functionResidual, s.gradientNorm);
                ImGui::TextWrapped("   Hessian diagnostic: %s", s.classification.c_str());
            }
        }

        if (!kRows_.empty() && ImGui::BeginTable("k sweep results", 5, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp)) {
            ImGui::TableSetupColumn("k");
            ImGui::TableSetupColumn("Vertices");
            ImGui::TableSetupColumn("Triangles");
            ImGui::TableSetupColumn("Components");
            ImGui::TableSetupColumn("Norm. symmetric RMS distance");
            ImGui::TableHeadersRow();
            for (const auto& row : kRows_) {
                ImGui::TableNextRow();
                ImGui::TableNextColumn(); ImGui::Text("%.6g", row.k);
                ImGui::TableNextColumn(); ImGui::Text("%llu", static_cast<unsigned long long>(row.vertices));
                ImGui::TableNextColumn(); ImGui::Text("%llu", static_cast<unsigned long long>(row.triangles));
                ImGui::TableNextColumn(); ImGui::Text("%llu", static_cast<unsigned long long>(row.components));
                ImGui::TableNextColumn(); ImGui::Text("%.6g", row.normalizedRmsDistance);
            }
            ImGui::EndTable();
        }
        ImGui::TextWrapped("k-sweep rule: test equally spaced nonzero k values and choose the one with the lowest sampled symmetric RMS mesh distance to F=0, normalized by the domain-cube diagonal. This usually favors k values closest to zero; it does not enforce regularity or topology preservation.");
        ImGui::TextWrapped("Research limitation: this is a classic Marching Cubes baseline with numerical singularity candidates. It does not yet certify local topology or implement a proven singularity-aware patch.");
    }
    ImGui::End();
}

void SingularityViewer::extractCurrent() {
    try {
        const Polynomial active = activePolynomialFromUi();
        if (automaticResolution_) {
            selectAutomaticResolution();
            return;
        }
        const GridSpec grid = gridFromUi();
        const auto start = std::chrono::steady_clock::now();
        Mesh result = extractMarchingCubes(active, grid);
        const auto end = std::chrono::steady_clock::now();
        const double seconds = std::chrono::duration<double>(end - start).count();
        if (result.empty()) {
            throw std::runtime_error("No isosurface was found in this domain. Try a larger domain or a different polynomial.");
        }
        currentMesh_ = std::move(result);
        currentActiveField_ = active;
        currentFieldDescription_ = activeDescription(functionText_, usePerturbation_, perturbationText_, k_);
        meshWasExtractedWithPerturbation_ = usePerturbation_;
        singularities_ = findSingularityCandidates(active, extent_).points;
        rebuildViewerMesh();
        updateMetrics(active, "Extraction finished in " + formatNumber(seconds, 4) + " seconds.");
    } catch (const std::exception& e) {
        setStatus(std::string("Extraction failed: ") + e.what());
    }
}

void SingularityViewer::selectAutomaticResolution() {
    try {
        const Polynomial active = activePolynomialFromUi();
        const double domainDiagonal = 2.0 * extent_ * std::sqrt(3.0);
        const std::vector<int> candidates{16, 32, 64, 128};
        Mesh previous;
        Mesh selectedMesh;
        int selectedCells = 0;
        bool converged = false;
        bool topologyStable = false;
        double measuredChange = std::numeric_limits<double>::quiet_NaN();
        const auto start = std::chrono::steady_clock::now();

        for (const int n : candidates) {
            setStatus("Automatic resolution search: extracting " + std::to_string(n) + " cells per axis...");
            Mesh candidate = extractMarchingCubes(active, {n, extent_});
            if (candidate.empty()) continue;

            selectedCells = n;
            selectedMesh = candidate;
            if (!previous.empty()) {
                measuredChange = approximateSymmetricMeshDistance(candidate, previous) / domainDiagonal;
                const MeshMetrics previousMetrics = computeMeshMetrics(previous);
                const MeshMetrics candidateMetrics = computeMeshMetrics(candidate);
                topologyStable =
                    previousMetrics.connectedComponents == candidateMetrics.connectedComponents &&
                    previousMetrics.eulerCharacteristic == candidateMetrics.eulerCharacteristic &&
                    previousMetrics.boundaryEdges == candidateMetrics.boundaryEdges &&
                    previousMetrics.nonManifoldEdges == candidateMetrics.nonManifoldEdges;
                if (std::isfinite(measuredChange) && measuredChange <= convergenceTolerance_ && topologyStable) {
                    converged = true;
                    selectedMesh = std::move(candidate);
                    break;
                }
            }
            previous = std::move(candidate);
        }
        if (selectedMesh.empty()) {
            throw std::runtime_error("No isosurface was found at any tested resolution in this domain.");
        }

        currentMesh_ = std::move(selectedMesh);
        cells_ = selectedCells;
        currentActiveField_ = active;
        currentFieldDescription_ = activeDescription(functionText_, usePerturbation_, perturbationText_, k_);
        meshWasExtractedWithPerturbation_ = usePerturbation_;
        singularities_ = findSingularityCandidates(active, extent_).points;
        rebuildViewerMesh();
        const auto end = std::chrono::steady_clock::now();
        const double seconds = std::chrono::duration<double>(end - start).count();
        std::string note;
        if (converged) {
            note = "Convergence heuristic passed at " + std::to_string(selectedCells) +
                " cells per axis; normalized sampled symmetric RMS distance=" + formatNumber(measuredChange) +
                "; topology signature stable (components, Euler characteristic, boundary/non-manifold edges). Runtime=" +
                formatNumber(seconds, 4) + " s.";
        } else {
            note = "Maximum non-empty tested resolution (" + std::to_string(selectedCells) +
                ") reached without satisfying both geometry tolerance and topology-signature stability. Last normalized sampled symmetric RMS distance=" +
                formatNumber(measuredChange) + "; topology signature stable at last comparison=" +
                (topologyStable ? "yes" : "no") + ". This remains a heuristic, not a correctness guarantee. Runtime=" +
                formatNumber(seconds, 4) + " s.";
        }
        updateMetrics(active, note);
    } catch (const std::exception& e) {
        setStatus(std::string("Automatic resolution search failed: ") + e.what());
    }
}

void SingularityViewer::detectSingularities() {
    try {
        const Polynomial original = basePolynomialFromUi();
        const auto result = findSingularityCandidates(original, extent_);
        singularities_ = result.points;
        std::ostringstream out;
        out << result.diagnostic << "\n"
            << "Search is numerical and does not guarantee completeness.\n"
            << "Candidate count: " << singularities_.size() << '\n';
        for (std::size_t i = 0; i < singularities_.size(); ++i) {
            const auto& s = singularities_[i];
            out << i + 1 << ": (" << s.position.x << ", " << s.position.y << ", " << s.position.z
                << "), |F|=" << s.functionResidual << ", |grad F|=" << s.gradientNorm
                << ", Hessian: " << s.classification << '\n';
        }
        status_ = out.str();
        if (!currentMesh_.empty()) updateMetrics(currentActiveField_, status_);
    } catch (const std::exception& e) {
        setStatus(std::string("Singularity search failed: ") + e.what());
    }
}

void SingularityViewer::runKSweep() {
    try {
        const Polynomial base = basePolynomialFromUi();
        const Polynomial g = Polynomial::parse(std::string(perturbationText_));
        if (kMax_ < kMin_) throw std::runtime_error("k maximum must be greater than or equal to k minimum.");
        const GridSpec grid = gridFromUi();
        const Mesh reference = extractMarchingCubes(base, grid);
        if (reference.empty()) throw std::runtime_error("The original F=0 surface is empty in this domain, so k candidates cannot be compared.");
        kRows_.clear();
        double bestK = 0.0;
        double bestDeviation = std::numeric_limits<double>::infinity();
        bool hasBest = false;
        const double diameter = 2.0 * extent_ * std::sqrt(3.0);
        const int sampleTotal = std::clamp(sampleCount_, 3, 41);

        for (int i = 0; i < sampleTotal; ++i) {
            const double fraction = static_cast<double>(i) / static_cast<double>(sampleTotal - 1);
            const double candidateK = kMin_ + (kMax_ - kMin_) * fraction;
            if (std::abs(candidateK) < 1e-12) continue;
            setStatus("k sweep: evaluating sample " + std::to_string(i + 1) + " of " + std::to_string(sampleTotal) + "...");
            const Polynomial candidateField = base + g * candidateK;
            const Mesh candidateMesh = extractMarchingCubes(candidateField, grid);
            if (candidateMesh.empty()) continue;
            const MeshMetrics metrics = computeMeshMetrics(candidateMesh);
            const double deviation = approximateSymmetricMeshDistance(candidateMesh, reference) / diameter;
            kRows_.push_back({candidateK, metrics.vertices, metrics.triangles, metrics.connectedComponents, deviation});
            if (std::isfinite(deviation) && deviation < bestDeviation) {
                bestK = candidateK;
                bestDeviation = deviation;
                hasBest = true;
            }
        }
        if (!hasBest) throw std::runtime_error("No nonzero k candidate produced a comparable mesh.");
        usePerturbation_ = true;
        k_ = bestK;
        const Polynomial selectedField = base + g * bestK;
        currentMesh_ = extractMarchingCubes(selectedField, grid);
        currentActiveField_ = selectedField;
        currentFieldDescription_ = activeDescription(functionText_, true, perturbationText_, k_);
        meshWasExtractedWithPerturbation_ = true;
        singularities_ = findSingularityCandidates(selectedField, extent_).points;
        rebuildViewerMesh();
        updateMetrics(selectedField, "k sweep complete. Selected the tested nonzero k with the lowest sampled normalized symmetric RMS mesh distance: k=" + formatNumber(bestK, 8) + ", distance=" + formatNumber(bestDeviation) + ". This objective favors closeness to the original geometry and does not prove regularity or topology preservation.");
    } catch (const std::exception& e) {
        setStatus(std::string("k sweep failed: ") + e.what());
    }
}

void SingularityViewer::exportCurrentPmpMesh() {
    try {
        exportWithPmp(currentMesh_, std::string(pmpExportPath_));
        setStatus("Mesh exported through PMP to: " + std::string(pmpExportPath_));
    } catch (const std::exception& e) {
        setStatus(std::string("PMP export failed: ") + e.what() + " If the mesh is singular/non-manifold, try raw OBJ export.");
    }
}

void SingularityViewer::exportCurrentRawObj() {
    try {
        exportRawObj(currentMesh_, std::string(rawExportPath_));
        setStatus("Raw OBJ exported to: " + std::string(rawExportPath_));
    } catch (const std::exception& e) {
        setStatus(std::string("Raw OBJ export failed: ") + e.what());
    }
}

void SingularityViewer::exportSweepCsv() {
    try {
        if (kRows_.empty()) throw std::runtime_error("Run a k sweep before exporting its table.");
        std::ofstream out(csvExportPath_);
        if (!out) throw std::runtime_error("Could not open the requested CSV file for writing.");
        out << "k,vertices,triangles,connected_components,normalized_symmetric_rms_mesh_distance_to_original\n";
        out << std::setprecision(17);
        for (const auto& row : kRows_) {
            out << row.k << ',' << row.vertices << ',' << row.triangles << ',' << row.components << ',' << row.normalizedRmsDistance << '\n';
        }
        if (!out) throw std::runtime_error("An error occurred while writing the CSV file.");
        setStatus("k sweep exported to: " + std::string(csvExportPath_));
    } catch (const std::exception& e) {
        setStatus(std::string("CSV export failed: ") + e.what());
    }
}

void SingularityViewer::rebuildViewerMesh() {
    mesh_.clear();
    rejectedDisplayFaces_ = 0;
    bool sharedIndexedRepresentation = true;
    std::vector<pmp::Vertex> vertices;
    vertices.reserve(currentMesh_.vertices.size());
    for (const auto& point : currentMesh_.vertices) {
        vertices.push_back(mesh_.add_vertex(pmp::Point(point.x, point.y, point.z)));
    }

    // First try a proper indexed SurfaceMesh. When the extracted triangle list
    // is manifold-compatible, PMP and the extractor now share the same V/E/F
    // connectivity and their statistics agree.
    try {
        for (const auto& triangle : currentMesh_.triangles) {
            if (triangle[0] >= vertices.size() || triangle[1] >= vertices.size() || triangle[2] >= vertices.size()) {
                sharedIndexedRepresentation = false;
                break;
            }
            const pmp::Face face = mesh_.add_triangle(vertices[triangle[0]], vertices[triangle[1]], vertices[triangle[2]]);
            if (!face.is_valid()) {
                sharedIndexedRepresentation = false;
                break;
            }
        }
    } catch (const std::exception&) {
        sharedIndexedRepresentation = false;
    }

    if (!sharedIndexedRepresentation) {
        // A manifold-oriented SurfaceMesh cannot encode every singular or
        // non-manifold triangle soup. Clear the partial attempt and make a
        // display-only mesh with one independent triangle per face. The source
        // isa::Mesh remains untouched and is still the research/export source.
        mesh_.clear();
        rejectedDisplayFaces_ = 0;
        for (const auto& triangle : currentMesh_.triangles) {
            if (triangle[0] >= currentMesh_.vertices.size() ||
                triangle[1] >= currentMesh_.vertices.size() ||
                triangle[2] >= currentMesh_.vertices.size()) {
                ++rejectedDisplayFaces_;
                continue;
            }
            const auto& p0 = currentMesh_.vertices[triangle[0]];
            const auto& p1 = currentMesh_.vertices[triangle[1]];
            const auto& p2 = currentMesh_.vertices[triangle[2]];
            const pmp::Vertex v0 = mesh_.add_vertex(pmp::Point(p0.x, p0.y, p0.z));
            const pmp::Vertex v1 = mesh_.add_vertex(pmp::Point(p1.x, p1.y, p1.z));
            const pmp::Vertex v2 = mesh_.add_vertex(pmp::Point(p2.x, p2.y, p2.z));
            try {
                const pmp::Face face = mesh_.add_triangle(v0, v1, v2);
                if (!face.is_valid()) ++rejectedDisplayFaces_;
            } catch (const std::exception&) {
                ++rejectedDisplayFaces_;
            }
        }
        viewerRepresentation_ = "triangle-soup display fallback (PMP duplicates triangle corners; topology counts intentionally differ)";
    } else {
        viewerRepresentation_ = "shared indexed PMP SurfaceMesh (source and viewer V/E/F should agree)";
    }

    gridOverlay_->setGrid(cells_, extent_);
    gridOverlay_->setEnabled(showGrid_);
    if (meshWasExtractedWithPerturbation_) {
        const Mesh originalMesh = extractMarchingCubes(basePolynomialFromUi(), gridFromUi());
        originalWireOverlay_->setMeshEdges(originalMesh);
        originalWireOverlay_->setEnabled(showOriginalWireframe_);
    } else {
        originalWireOverlay_->clear();
        originalWireOverlay_->setEnabled(false);
    }

    update_mesh();
    if (!currentMesh_.empty()) {
        set_scene(pmp::vec3(0.0, 0.0, 0.0), static_cast<float>(extent_ * 1.8));
        view_all();
    }
}

void SingularityViewer::updateMetrics(const Polynomial& activeField, const std::string& note) {
    const MeshMetrics m = computeMeshMetrics(currentMesh_, &activeField);
    std::ostringstream out;
    out << "Active field used for the displayed mesh: " << currentFieldDescription_ << '\n'
        << "Grid: " << cells_ << " cells per axis; domain [-" << extent_ << ", " << extent_ << "]^3\n"
        << "\nExtractor triangle-list statistics: vertices=" << m.vertices << " | edges=" << m.edges
        << " | faces/triangles=" << m.triangles << '\n'
        << "PMP Mesh Info (actual viewer SurfaceMesh): vertices=" << mesh_.n_vertices()
        << " | edges=" << mesh_.n_edges() << " | faces=" << mesh_.n_faces() << '\n'
        << "PMP display representation: " << viewerRepresentation_ << '\n';
    if (rejectedDisplayFaces_ > 0) {
        out << "WARNING: PMP viewer could not represent " << rejectedDisplayFaces_
            << " source triangles; raw OBJ retains the complete source triangle list.\n";
    }
    out << "Components: " << m.connectedComponents << " | Euler characteristic: " << m.eulerCharacteristic << '\n'
        << "Boundary edges: " << m.boundaryEdges << " | Non-manifold edges: " << m.nonManifoldEdges
        << " | Degenerate triangles: " << m.degenerateTriangles << '\n'
        << "Approx. RMS first-order residual distance to active field: " << formatNumber(m.implicitResidualRms)
        << " (samples=" << m.residualSamples << ", skipped near-critical=" << m.skippedNearCriticalSamples << ")\n"
        << "Numerical singularity candidates currently shown: " << singularities_.size() << '\n';
    if (!note.empty()) out << "\n" << note << '\n';
    out << "\nResearch limitation: classic Marching Cubes and numerical diagnostics do not certify the exact topology of the polynomial zero set.";
    metricsText_ = out.str();
    if (!note.empty()) status_ = note;
}

} // namespace isa
