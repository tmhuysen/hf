# top-level CMakeLists

cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(hartree-fock LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)


# Find the Eigen3 package (http://eigen.tuxfamily.org/dox/TopicCMakeGuide.html)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)


# source files for the library
add_subdirectory(src)


# Set up testing
#enable_testing()
#find_package(Boost REQUIRED)
#add_subdirectory(tests)
