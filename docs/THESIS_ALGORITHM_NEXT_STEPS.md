# Next steps for the thesis-specific extraction method

The current app is a baseline and experiment harness. The main thesis contribution—an extraction method with a demonstrated local-topology guarantee—still needs to be designed and validated. Keep the following stages separate from the existing classical Marching Cubes baseline so all comparisons remain reproducible.

## Research target

Start with isolated nondegenerate quadratic-type singularities of polynomial surfaces. At a candidate point `p`, require `F(p)=0` and `grad(F)(p)=0`; classify the Hessian signature. For an indefinite nondegenerate Hessian in three variables, the quadratic leading term has a cone-like zero set. Do not extend any claim to degenerate or non-isolated singularities unless separately tested.

## Recommended next implementation increments

1. **Improve singularity validation.** Check the candidate residual and gradient using scale-aware tolerances. Repeat the search with several seed grids and report agreement/disagreement. For thesis-grade completeness, compare against symbolic or certified root isolation on small test polynomials.
2. **Add a local link test.** Intersect a small neighborhood of the candidate with a sphere centered at the singularity. For the model `x^2+y^2-z^2=0`, the link is two circles. Sample multiple radii to check stability and do not rely on a single mesh's Euler characteristic alone.
3. **Add local refinement first.** Mark cells within a configurable multiple of the local grid spacing from a candidate and refine them. Record the number of cells and triangles saved versus uniform refinement. Ensure neighboring levels are connected without cracks.
4. **Design a local patch separately.** Use the quadratic Taylor model to characterize local branches. Define the boundary of the patch region, compute how the reference surface intersects that boundary, and construct a patch that matches those boundary curves. Test it first on the exact quadratic cone before higher-degree examples.
5. **Do not silently smooth the singular point.** If a patch regularizes the surface, mark that as a different mode and measure the topological change relative to the original zero set.
6. **Validate grid-position robustness.** Translate the grid relative to the same algebraic surface, rotate/scale the test surfaces, and compare the singular-point neighborhood, not just the global triangle count.
7. **Compare to multiple baselines.** At minimum use classic Marching Cubes and an ambiguity-resolving variant such as MC33/Lewiner. If practical, add Marching Tetrahedra or Dual Contouring.

## Minimum validation suite

- Smooth controls: sphere and ellipsoid.
- Target: `x^2+y^2-z^2=0`.
- Coordinate transforms: scaled and rotated copies of the target.
- Higher-degree local-equivalent model: `x^2+y^2-z^2+a*z^3=0` for a small, documented `a` and domain that keeps the tested local feature visible.
- Stress cases: crossing planes `x*y=0` and cusp-like `z^2-x^3-y^2=0`.

## Acceptance criteria to define before running experiments

For every run, record geometric deviation from the original polynomial surface, local link connectivity, connected components, boundary/non-manifold diagnostics, triangle count, grid-cell count and runtime. Define numerical tolerances before comparing algorithms. A successful run on one alignment or one resolution is not enough to support a general correctness claim.
