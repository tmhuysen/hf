# In this CMake file, we will find all required packages

# Find the Eigen3 package - needed for this project, and for libint2
# Options to find this package from (http://eigen.tuxfamily.org/dox/TopicCMakeGuide.html)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Find the boost package - needed for unittests
find_package(Boost REQUIRED)

# Find my own libint wrapper
find_package(libwrp 2.1.1 REQUIRED)
