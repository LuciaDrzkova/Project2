#include "isa/SingularityViewer.h"

#include <exception>
#include <iostream>

int main() {
    try {
        isa::SingularityViewer viewer;
        return viewer.run();
    } catch (const std::exception& e) {
        std::cerr << "Singularity-aware isosurface application failed: " << e.what() << '\n';
        return 1;
    }
}
