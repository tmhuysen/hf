# top-level CMakeLists

cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(hartree-fock LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)



add_executable(hartree-fock
        SCF_Solver.cpp)


# add the libint-eigen headers
target_include_directories(hartree-fock PRIVATE /opt/local/include/libint-eigen)

# add the libint headers
target_include_directories(hartree-fock PRIVATE /usr/local/libint/2.2.0/include)
target_include_directories(hartree-fock PRIVATE /usr/local/libint/2.2.0/include/libint2)

# add the eigen headers
target_include_directories(hartree-fock PRIVATE /opt/local/include/eigen3)


# link to libint
target_link_libraries(hartree-fock /usr/local/libint/2.2.0/lib/libint2.a)

# link to libint-eigen-wrapper
target_link_libraries(hartree-fock /opt/local/lib/liblibint-eigen.a)
