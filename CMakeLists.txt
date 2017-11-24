# top-level CMakeLists
#   CMakeLists.txt files modified from (https://github.com/pablospe/cmake-example-library)
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)

# Specify the project, its version and the used language
project(hf
        VERSION 1.0.1
        LANGUAGES CXX)

# To enable the C++11 standard, we have to add flags like -std=gnu+11 or so. CMake can take care of this:
set(CMAKE_CXX_STANDARD 11)



# Set CMake variables
include(${CMAKE_SOURCE_DIR}/cmake/SetCMakeEnvironment.cmake)

# Install Findlibint2.cmake. This is done before FindPackages.cmake, since we are using find_package(libint2)
# CMake looks for Find<package>.cmake files in its ${CMAKE_MODULE_PATH}. By default, this is empty. Since we don't want to bother with installing our custom-made Findlibint2.cmake, we will set the ${CMAKE_MODULE_PATH} accordingly
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_SOURCE_DIR}/cmake")

# Find all required packages
include(${CMAKE_SOURCE_DIR}/cmake/FindPackages.cmake)



# Add the sources for the library
add_subdirectory(${PROJECT_SOURCE_FOLDER})



# Set up testing
enable_testing()
add_subdirectory(${PROJECT_TESTS_FOLDER})
