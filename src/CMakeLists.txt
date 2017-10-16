# src-level CMakeLists


# The goal of this project is to create a library called 'hartree-fock'
add_library(hartree-fock
        SCFSolver.cpp)


# Include this project's headers
target_include_directories(hartree-fock PRIVATE ../include)

# add my own libint-wrapper headers
target_include_directories(hartree-fock PRIVATE /opt/local/include/libint-wrapper)

# add the libint headers
target_include_directories(hartree-fock PRIVATE /usr/local/libint/2.2.0/include)
target_include_directories(hartree-fock PRIVATE /usr/local/libint/2.2.0/include/libint2)

# add the eigen headers (for libint)
target_link_libraries(hartree-fock PRIVATE Eigen3::Eigen)


# link to libint
target_link_libraries(hartree-fock PRIVATE /usr/local/libint/2.2.0/lib/libint2.a)

# link to my own libint-wrapper
target_link_libraries(hartree-fock PRIVATE /opt/local/lib/liblibint-wrapper.a)
