# Validation status — PMP viewer edition

- The dependency-light mathematical core (`Polynomial`, classic Marching Cubes, numerical singularity search and mesh metrics) configured, compiled and passed its regression tests in the available build environment with `ISA_BUILD_APP=OFF`.
- The revised GUI code adds an indexed-mesh attempt with triangle-soup fallback, native PMP statistics in the diagnostics, OpenGL grid/original-surface overlays, clearer k-sweep metric labeling and a topology-signature check in automatic resolution selection.
- The full GUI was **not compiled or launched in this environment** because PMP's viewer dependencies are not available locally and network access from the build container is unavailable. In particular, the new OpenGL overlay should be confirmed by building and running the app on the target Mac.
- The application uses no Qt headers or CMake `find_package(Qt6)` calls.
- This remains a classic Marching Cubes experimental baseline; it does not claim a proven singularity-preserving local extraction method.

On macOS, run a full CMake configure/build as described in the README. The first build is the definitive check of the local Xcode/OpenGL/PMP toolchain.
