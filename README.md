# hf v3.0.0

[![Build Status](https://travis-ci.org/GQCG/hf.svg?branch=develop)](https://travis-ci.org/GQCG/hf)

A C++ library that performs Hartree-Fock calculations using Eigen3 for matrix manipulations, based on a wrapper around libint2.


## Dependencies
[![Eigen3 Dependency](https://img.shields.io/badge/Eigen-3.3.4+-000000.svg)](http://eigen.tuxfamily.org/index.php?title=Main_Page)
[![Boost Dependency](https://img.shields.io/badge/Boost-1.65.1+-000000.svg)](www.boost.org)
[![libint2 Dependency](https://img.shields.io/badge/libint-2.3.1+-000000.svg)](https://github.com/evaleev/libint)

[![cpputil Dependency](https://img.shields.io/badge/cpputil-1.2.1+-blue.svg)](https://github.com/GQCG/cpputil)
[![libwint Dependency](https://img.shields.io/badge/libwint-3.0.0+-blue.svg)](https://github.com/GQCG/libwint)

You might have to copy the sto-6g.g94 (provided in `docs/sto-6g.g94`) file to the relevant directory. In a default installation, the file should be copied to `/usr/local/libint/2.3.1/share/libint/2.3.1/basis/sto-6g.g94`.


## Installation
To install this library:
1. clone the master branch, which contains the latest release

        git clone https://github.com/GQCG/hf.git --branch master --single-branch
        cd hf

2. perform an out-of-source cmake build:

        mkdir build && cd build
        cmake -DINSTALLATION_PREFIX=prefix ..
        make && make test && sudo make install

    where
    * `prefix` is the installation prefix (defaulted to `/usr/local`) you want the library to be installed at:
        * the library `libhf.a` will be installed in `prefix/hf/lib`
        * the header files (and cmake files, see Usage) will be installed in `prefix/hf/include`


## Usage
Basic usage of this library can be found in the `tests` directory. If you use CMake in other projects, you can add the following CMake command to the CMakeLists.txt-file:

    find_package(hf 3.0.0)

where `x.y.z` is the version number. CMake then provides the commands `hf_INCLUDE_DIRS` to be used in your `target_include_directories` and the library `hf` to be used in your `target_link_libraries`.
