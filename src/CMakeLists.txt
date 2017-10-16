# src-level CMakeLists

add_library(hartree-fock
        SCFSolver.cpp
        SCF_functions.cpp)


# Include this project's headers
target_include_directories(hartree-fock PRIVATE ../include)

# Add my own libint wrapper
target_include_directories(hartree-fock PRIVATE /usr/local/include/libint-wrapper)
target_link_libraries(hartree-fock PRIVATE /usr/local/lib/liblibint-wrapper.a)

# Include the Eigen headers
target_link_libraries(hartree-fock PRIVATE Eigen3::Eigen)

# Add libint
target_include_directories(hartree-fock PRIVATE /usr/local/libint/2.2.0/include)
target_include_directories(hartree-fock PRIVATE /usr/local/libint/2.2.0/include/libint2)
target_link_libraries(hartree-fock PRIVATE /usr/local/libint/2.2.0/lib/libint2.a)
