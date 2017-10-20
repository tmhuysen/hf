# In this CMake file, we will include the headers and link to the necessary libraries


# Include this project's headers
target_include_directories(${LIBRARY_NAME} PRIVATE ${PROJECT_INCLUDE_FOLDER})

# Include the Eigen3 headers
target_link_libraries(${LIBRARY_NAME} PUBLIC Eigen3::Eigen)

# Add libint2
target_include_directories(${LIBRARY_NAME} PRIVATE ${libint2_INCLUDE_DIRS})
target_link_libraries(${LIBRARY_NAME} PUBLIC ${libint2_LIBRARIES})


# Include the libint wrapper
target_include_directories(${LIBRARY_NAME} PRIVATE ${libwrp_INCLUDE_DIRS})
target_link_libraries(${LIBRARY_NAME} PRIVATE libwrp)
