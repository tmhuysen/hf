# In this CMake file, we will find all required packages

# Find the Eigen3 package - needed for this project, and for libint2
# Options to find this package from (http://eigen.tuxfamily.org/dox/TopicCMakeGuide.html)
find_package(Eigen3 3.3.4 REQUIRED)

# Find the boost package - needed for unittests
find_package(Boost REQUIRED)

# Find my own libint wrapper
find_package(libwint 3.0.0 REQUIRED)

# Find cpputil
find_package(cpputil 1.2.1 REQUIRED)
