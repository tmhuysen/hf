# Install script for directory: /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/hf/lib/libhf.a")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/hf/lib" TYPE STATIC_LIBRARY FILES "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/in/src/libhf.a")
  if(EXISTS "$ENV{DESTDIR}/usr/local/hf/lib/libhf.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/hf/lib/libhf.a")
    execute_process(COMMAND "/opt/local/bin/ranlib" "$ENV{DESTDIR}/usr/local/hf/lib/libhf.a")
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/hf/include/BaseSCFSolver.hpp;/usr/local/hf/include/DIISSCFSolver.hpp;/usr/local/hf/include/PlainSCFSolver.hpp;/usr/local/hf/include/RHF.hpp;/usr/local/hf/include/RHFC.hpp;/usr/local/hf/include/SCFSolverType.hpp;/usr/local/hf/include/common.hpp;/usr/local/hf/include/hf.hpp;/usr/local/hf/include/version.hpp")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/hf/include" TYPE FILE FILES
    "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/include/BaseSCFSolver.hpp"
    "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/include/DIISSCFSolver.hpp"
    "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/include/PlainSCFSolver.hpp"
    "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/include/RHF.hpp"
    "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/include/RHFC.hpp"
    "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/include/SCFSolverType.hpp"
    "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/include/common.hpp"
    "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/include/hf.hpp"
    "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/include/version.hpp"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/usr/local/hf/cmake/hfTargets.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}/usr/local/hf/cmake/hfTargets.cmake"
         "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/in/src/CMakeFiles/Export/_usr/local/hf/cmake/hfTargets.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}/usr/local/hf/cmake/hfTargets-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}/usr/local/hf/cmake/hfTargets.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/hf/cmake/hfTargets.cmake")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/hf/cmake" TYPE FILE FILES "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/in/src/CMakeFiles/Export/_usr/local/hf/cmake/hfTargets.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
     "/usr/local/hf/cmake/hfTargets-release.cmake")
    if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
        message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
    endif()
    if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
        message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
    endif()
file(INSTALL DESTINATION "/usr/local/hf/cmake" TYPE FILE FILES "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/in/src/CMakeFiles/Export/_usr/local/hf/cmake/hfTargets-release.cmake")
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/hf/cmake/hfConfig.cmake;/usr/local/hf/cmake/hfConfigVersion.cmake")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/hf/cmake" TYPE FILE FILES
    "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/cmake/hfConfig.cmake"
    "/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/cmake/hfConfigVersion.cmake"
    )
endif()

