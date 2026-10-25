# Validation status — PMP viewer edition

- The dependency-light mathematical core (`Polynomial`, classic Marching Cubes, numerical singularity search and mesh metrics) has been configured, compiled and tested with `ISA_BUILD_APP=OFF`.
- The automated core regression suite is expected to run with the commands in the README; the original core suite passed in the preceding Qt-based package before replacing the application layer.
- The GUI has been migrated from Qt Widgets to PMP's `MeshViewer` and ImGui API based on PMP 3.0.0. The GUI and PMP adapter could not be compiled in this execution environment because PMP's viewer dependencies are not installed and network access from the build container is unavailable.
- The application uses no Qt headers or CMake `find_package(Qt6)` calls.
- This remains a classic Marching Cubes experimental baseline; it does not claim a proven singularity-preserving local extraction method.

On the target Mac, run a full CMake configure/build as described in the README. The first build is the definitive check of the local Xcode/OpenGL/PMP toolchain.
