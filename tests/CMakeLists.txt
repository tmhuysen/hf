# tests-level CMakeLists

# For every *.cpp-file (every file corresponding to one set of tests) ...
file(GLOB test_sources ./*.cpp)

foreach(test_source ${test_sources})
    # ... extract the filename without extension as a name for our executable ...
    get_filename_component(test_name ${test_source} NAME_WE)

    # ... add an executable based on the test source ...
    add_executable(${test_name} ${test_source})

    # ... add the boost headers (for testing) ...
    target_include_directories(${test_name} PRIVATE ${Boost_INCLUDE_DIRS})

    # ... add this project's headers ...
    target_include_directories(${test_name} PRIVATE ../include)

    # ... add the eigen3 headers (for libint) ...
    target_link_libraries(${test_name} PUBLIC Eigen3::Eigen)

    # ... add the libint2 headers ...
    target_include_directories(${test_name} PRIVATE /usr/local/libint/2.2.0/include)
    target_include_directories(${test_name} PRIVATE /usr/local/libint/2.2.0/include/libint2)

    # ... link to this project's library ...
    target_link_libraries(${test_name} PRIVATE hartree-fock)

    # ... add my own libint-wrapper headers ...
    target_include_directories(${test_name} PRIVATE /opt/local/include/libint-wrapper)


    # ... link to libint ...
    target_link_libraries(${test_name} /usr/local/libint/2.2.0/lib/libint2.a)

    # ... link to my own libint-wrapper ...
    target_link_libraries(${test_name} /opt/local/lib/liblibint-wrapper.a)


    # ... and finally make a test ...
    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()