#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "hf" for configuration "Release"
set_property(TARGET hf APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hf PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "/usr/local/hf/lib/libhf.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS hf )
list(APPEND _IMPORT_CHECK_FILES_FOR_hf "/usr/local/hf/lib/libhf.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
