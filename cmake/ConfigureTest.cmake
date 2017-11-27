# In this CMake file, we will include header files and link to libraries for a given test source

# ... add the boost headers (for testing) ...
target_include_directories(${TEST_NAME} PRIVATE ${Boost_INCLUDE_DIRS})

# ... add this project's headers ...
target_include_directories(${TEST_NAME} PRIVATE ${PROJECT_INCLUDE_FOLDER})

# ... add the Eigen3 headers ...
target_link_libraries(${TEST_NAME} PUBLIC Eigen3::Eigen)

# ... add the libint wrapper ...
target_include_directories(${TEST_NAME} PUBLIC ${libwrp_INCLUDE_DIRS})
target_link_libraries(${TEST_NAME} PUBLIC libwrp)

# ... link to this project's library ...
target_link_libraries(${TEST_NAME} PRIVATE ${LIBRARY_NAME})
