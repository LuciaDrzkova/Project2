cmake_minimum_required(VERSION 3.24)
project(SingularityAwareIsosurfaces VERSION 0.2.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(ISA_BUILD_APP "Build the PMP-native desktop application (no Qt)" ON)
option(ISA_FETCH_PMP "Fetch the pinned PMP 3.0.0 dependency for the app" ON)
option(ISA_PMP_USE_DOUBLE "Set ON when the selected PMP build uses 64-bit scalar precision" OFF)
option(ISA_BUILD_TESTS "Build mathematical core tests" ON)

add_library(isa_core
    src/Polynomial.cpp
    src/MarchingCubes.cpp
    src/Singularity.cpp
    src/Metrics.cpp)
target_include_directories(isa_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(isa_core PUBLIC cxx_std_20)
if(MSVC)
    target_compile_options(isa_core PRIVATE /W4)
else()
    target_compile_options(isa_core PRIVATE -Wall -Wextra -Wpedantic)
endif()

if(ISA_BUILD_APP)
    # PMP's own viewer uses GLFW + OpenGL + Dear ImGui. No Qt installation is needed.
    set(ISA_PMP_SOURCE_DIR "" CACHE PATH "Optional path to an existing PMP source checkout (for example external/pmp-library)")
    message(STATUS "ISA_PMP_SOURCE_DIR = '${ISA_PMP_SOURCE_DIR}'")
    message(STATUS "ISA_FETCH_PMP = ${ISA_FETCH_PMP}")  
    if(ISA_PMP_SOURCE_DIR)
        set(PMP_SCALAR_TYPE 64 CACHE STRING "Use double precision in PMP" FORCE)
        set(PMP_BUILD_VIEWERS ON CACHE BOOL "Build PMP's GLFW/ImGui viewers" FORCE)
        set(PMP_BUILD_EXAMPLES OFF CACHE BOOL "Do not build PMP examples" FORCE)
        set(PMP_BUILD_TESTS OFF CACHE BOOL "Do not build PMP tests" FORCE)
        set(PMP_BUILD_DOCS OFF CACHE BOOL "Do not build PMP docs" FORCE)
        set(PMP_BUILD_REGRESSIONS OFF CACHE BOOL "Do not build PMP regressions" FORCE)
        set(PMP_INSTALL OFF CACHE BOOL "Consume PMP as a subproject" FORCE)
        set(PMP_STRICT_COMPILATION OFF CACHE BOOL "Do not treat dependency warnings as errors" FORCE)
        add_subdirectory("${ISA_PMP_SOURCE_DIR}" "${CMAKE_BINARY_DIR}/pmp-build")
        set(ISA_PMP_USE_DOUBLE ON CACHE BOOL "The local PMP source is configured for double precision" FORCE)
    elseif(ISA_FETCH_PMP)
        include(FetchContent)
        set(PMP_SCALAR_TYPE 64 CACHE STRING "Use double precision in PMP" FORCE)
        set(PMP_BUILD_VIEWERS ON CACHE BOOL "Build PMP's GLFW/ImGui viewers" FORCE)
        set(PMP_BUILD_EXAMPLES OFF CACHE BOOL "Do not build PMP examples" FORCE)
        set(PMP_BUILD_TESTS OFF CACHE BOOL "Do not build PMP tests" FORCE)
        set(PMP_BUILD_DOCS OFF CACHE BOOL "Do not build PMP docs" FORCE)
        set(PMP_BUILD_REGRESSIONS OFF CACHE BOOL "Do not build PMP regressions" FORCE)
        set(PMP_INSTALL OFF CACHE BOOL "Consume PMP as a subproject" FORCE)
        set(PMP_STRICT_COMPILATION OFF CACHE BOOL "Do not treat dependency warnings as errors" FORCE)
        FetchContent_Declare(pmp
            GIT_REPOSITORY https://github.com/pmp-library/pmp-library.git
            GIT_TAG main
            GIT_SHALLOW TRUE)
        FetchContent_MakeAvailable(pmp)
        set(ISA_PMP_USE_DOUBLE ON CACHE BOOL "The fetched PMP dependency uses double precision" FORCE)
    else()
        find_package(pmp CONFIG REQUIRED)
    endif()

    if(NOT TARGET pmp)
        message(FATAL_ERROR "PMP was not found. Set ISA_PMP_SOURCE_DIR, leave ISA_FETCH_PMP=ON, or install a compatible PMP package.")
    endif()
    if(NOT TARGET pmp_viewers)
        message(FATAL_ERROR "The selected PMP build does not include pmp_viewers. Rebuild/install PMP with -DPMP_BUILD_VIEWERS=ON, or leave ISA_FETCH_PMP=ON so this project fetches a compatible PMP build.")
    endif()

    add_executable(singularity_aware_isosurfaces
        src/main.cpp
        src/SingularityViewer.cpp
        src/PmpExport.cpp
        include/isa/SingularityViewer.h)
    target_link_libraries(singularity_aware_isosurfaces PRIVATE isa_core pmp pmp_viewers)
    if(ISA_PMP_USE_DOUBLE)
        target_compile_definitions(singularity_aware_isosurfaces PRIVATE PMP_SCALAR_TYPE_64)
    endif()
    if(MSVC)
        target_compile_options(singularity_aware_isosurfaces PRIVATE /W4)
    else()
        target_compile_options(singularity_aware_isosurfaces PRIVATE -Wall -Wextra -Wpedantic)
    endif()
endif()

if(ISA_BUILD_TESTS)
    include(CTest)
    add_executable(isa_core_tests tests/core_tests.cpp)
    target_link_libraries(isa_core_tests PRIVATE isa_core)
    add_test(NAME isa_core_tests COMMAND isa_core_tests)
endif()
