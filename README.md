# hf

[![Build Status](https://travis-ci.org/lelemmen/hf.svg?branch=master)](https://travis-ci.org/lelemmen/hf)

A C++ library that performs Hartree-Fock calculations using Eigen3 for matrix manipulations, based on a wrapper around libint2.


## Dependencies
[![libint2 Dependency](https://img.shields.io/badge/libint-2.3.1+-blue.svg)](https://github.com/evaleev/libint)
[![Eigen3 Dependency](https://img.shields.io/badge/Eigen-3+-blue.svg)](http://eigen.tuxfamily.org/index.php?title=Main_Page)
[![libwrp Dependency](https://img.shields.io/badge/libwrp-2.1.1+-blue.svg)](https://github.com/lelemmen/libwrp)
You might have to copy the sto-6g.g94 (provided in `tests/reference_data/sto-6g.g94`) file to the relevant directory. In a default installation, the file should be copied to `/usr/local/libint/2.3.1/share/libint/2.3.1/basis/sto-6g.g94`.


## Installation
To install this library:
1. clone the master branch

        git clone https://github.com/lelemmen/hf.git
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

    find_package(hf x.y.z)

where `x.y.z` is the version number. CMake then provides the commands `hf_INCLUDE_DIRS` to be used in your `target_include_directories` and the library `hf` to be used in your `target_link_libraries`.
