# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.11

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/local/bin/cmake

# The command to remove a file.
RM = /opt/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/b3

# Include any dependencies generated for this target.
include tests/CMakeFiles/RHFC_DIIS_test.dir/depend.make

# Include the progress variables for this target.
include tests/CMakeFiles/RHFC_DIIS_test.dir/progress.make

# Include the compile flags for this target's objects.
include tests/CMakeFiles/RHFC_DIIS_test.dir/flags.make

tests/CMakeFiles/RHFC_DIIS_test.dir/RHFC_DIIS_test.cpp.o: tests/CMakeFiles/RHFC_DIIS_test.dir/flags.make
tests/CMakeFiles/RHFC_DIIS_test.dir/RHFC_DIIS_test.cpp.o: ../tests/RHFC_DIIS_test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/b3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests/CMakeFiles/RHFC_DIIS_test.dir/RHFC_DIIS_test.cpp.o"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/b3/tests && /opt/local/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/RHFC_DIIS_test.dir/RHFC_DIIS_test.cpp.o -c /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/tests/RHFC_DIIS_test.cpp

tests/CMakeFiles/RHFC_DIIS_test.dir/RHFC_DIIS_test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/RHFC_DIIS_test.dir/RHFC_DIIS_test.cpp.i"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/b3/tests && /opt/local/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/tests/RHFC_DIIS_test.cpp > CMakeFiles/RHFC_DIIS_test.dir/RHFC_DIIS_test.cpp.i

tests/CMakeFiles/RHFC_DIIS_test.dir/RHFC_DIIS_test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/RHFC_DIIS_test.dir/RHFC_DIIS_test.cpp.s"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/b3/tests && /opt/local/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/tests/RHFC_DIIS_test.cpp -o CMakeFiles/RHFC_DIIS_test.dir/RHFC_DIIS_test.cpp.s

# Object files for target RHFC_DIIS_test
RHFC_DIIS_test_OBJECTS = \
"CMakeFiles/RHFC_DIIS_test.dir/RHFC_DIIS_test.cpp.o"

# External object files for target RHFC_DIIS_test
RHFC_DIIS_test_EXTERNAL_OBJECTS =

tests/RHFC_DIIS_test: tests/CMakeFiles/RHFC_DIIS_test.dir/RHFC_DIIS_test.cpp.o
tests/RHFC_DIIS_test: tests/CMakeFiles/RHFC_DIIS_test.dir/build.make
tests/RHFC_DIIS_test: src/libhf.a
tests/RHFC_DIIS_test: /usr/local/cpputil/lib/libcpputil.a
tests/RHFC_DIIS_test: /usr/local/libwint/lib/libwint.a
tests/RHFC_DIIS_test: /usr/local/libint/2.3.1/lib/libint2.a
tests/RHFC_DIIS_test: tests/CMakeFiles/RHFC_DIIS_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/b3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable RHFC_DIIS_test"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/b3/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/RHFC_DIIS_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/CMakeFiles/RHFC_DIIS_test.dir/build: tests/RHFC_DIIS_test

.PHONY : tests/CMakeFiles/RHFC_DIIS_test.dir/build

tests/CMakeFiles/RHFC_DIIS_test.dir/clean:
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/b3/tests && $(CMAKE_COMMAND) -P CMakeFiles/RHFC_DIIS_test.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/RHFC_DIIS_test.dir/clean

tests/CMakeFiles/RHFC_DIIS_test.dir/depend:
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/b3 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/tests /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/b3 /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/b3/tests /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/b3/tests/CMakeFiles/RHFC_DIIS_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/RHFC_DIIS_test.dir/depend
