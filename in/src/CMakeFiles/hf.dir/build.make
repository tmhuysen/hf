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
CMAKE_BINARY_DIR = /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/in

# Include any dependencies generated for this target.
include src/CMakeFiles/hf.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/hf.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/hf.dir/flags.make

src/CMakeFiles/hf.dir/BaseSCFSolver.cpp.o: src/CMakeFiles/hf.dir/flags.make
src/CMakeFiles/hf.dir/BaseSCFSolver.cpp.o: ../src/BaseSCFSolver.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/in/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/hf.dir/BaseSCFSolver.cpp.o"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/in/src && /opt/local/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/hf.dir/BaseSCFSolver.cpp.o -c /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/src/BaseSCFSolver.cpp

src/CMakeFiles/hf.dir/BaseSCFSolver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/hf.dir/BaseSCFSolver.cpp.i"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/in/src && /opt/local/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/src/BaseSCFSolver.cpp > CMakeFiles/hf.dir/BaseSCFSolver.cpp.i

src/CMakeFiles/hf.dir/BaseSCFSolver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/hf.dir/BaseSCFSolver.cpp.s"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/in/src && /opt/local/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/src/BaseSCFSolver.cpp -o CMakeFiles/hf.dir/BaseSCFSolver.cpp.s

src/CMakeFiles/hf.dir/DIISSCFSolver.cpp.o: src/CMakeFiles/hf.dir/flags.make
src/CMakeFiles/hf.dir/DIISSCFSolver.cpp.o: ../src/DIISSCFSolver.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/in/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/CMakeFiles/hf.dir/DIISSCFSolver.cpp.o"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/in/src && /opt/local/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/hf.dir/DIISSCFSolver.cpp.o -c /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/src/DIISSCFSolver.cpp

src/CMakeFiles/hf.dir/DIISSCFSolver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/hf.dir/DIISSCFSolver.cpp.i"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/in/src && /opt/local/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/src/DIISSCFSolver.cpp > CMakeFiles/hf.dir/DIISSCFSolver.cpp.i

src/CMakeFiles/hf.dir/DIISSCFSolver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/hf.dir/DIISSCFSolver.cpp.s"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/in/src && /opt/local/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/src/DIISSCFSolver.cpp -o CMakeFiles/hf.dir/DIISSCFSolver.cpp.s

src/CMakeFiles/hf.dir/PlainSCFSolver.cpp.o: src/CMakeFiles/hf.dir/flags.make
src/CMakeFiles/hf.dir/PlainSCFSolver.cpp.o: ../src/PlainSCFSolver.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/in/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object src/CMakeFiles/hf.dir/PlainSCFSolver.cpp.o"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/in/src && /opt/local/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/hf.dir/PlainSCFSolver.cpp.o -c /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/src/PlainSCFSolver.cpp

src/CMakeFiles/hf.dir/PlainSCFSolver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/hf.dir/PlainSCFSolver.cpp.i"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/in/src && /opt/local/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/src/PlainSCFSolver.cpp > CMakeFiles/hf.dir/PlainSCFSolver.cpp.i

src/CMakeFiles/hf.dir/PlainSCFSolver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/hf.dir/PlainSCFSolver.cpp.s"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/in/src && /opt/local/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/src/PlainSCFSolver.cpp -o CMakeFiles/hf.dir/PlainSCFSolver.cpp.s

src/CMakeFiles/hf.dir/RHF.cpp.o: src/CMakeFiles/hf.dir/flags.make
src/CMakeFiles/hf.dir/RHF.cpp.o: ../src/RHF.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/in/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object src/CMakeFiles/hf.dir/RHF.cpp.o"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/in/src && /opt/local/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/hf.dir/RHF.cpp.o -c /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/src/RHF.cpp

src/CMakeFiles/hf.dir/RHF.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/hf.dir/RHF.cpp.i"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/in/src && /opt/local/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/src/RHF.cpp > CMakeFiles/hf.dir/RHF.cpp.i

src/CMakeFiles/hf.dir/RHF.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/hf.dir/RHF.cpp.s"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/in/src && /opt/local/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/src/RHF.cpp -o CMakeFiles/hf.dir/RHF.cpp.s

src/CMakeFiles/hf.dir/RHFC.cpp.o: src/CMakeFiles/hf.dir/flags.make
src/CMakeFiles/hf.dir/RHFC.cpp.o: ../src/RHFC.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/in/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object src/CMakeFiles/hf.dir/RHFC.cpp.o"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/in/src && /opt/local/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/hf.dir/RHFC.cpp.o -c /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/src/RHFC.cpp

src/CMakeFiles/hf.dir/RHFC.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/hf.dir/RHFC.cpp.i"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/in/src && /opt/local/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/src/RHFC.cpp > CMakeFiles/hf.dir/RHFC.cpp.i

src/CMakeFiles/hf.dir/RHFC.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/hf.dir/RHFC.cpp.s"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/in/src && /opt/local/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/src/RHFC.cpp -o CMakeFiles/hf.dir/RHFC.cpp.s

# Object files for target hf
hf_OBJECTS = \
"CMakeFiles/hf.dir/BaseSCFSolver.cpp.o" \
"CMakeFiles/hf.dir/DIISSCFSolver.cpp.o" \
"CMakeFiles/hf.dir/PlainSCFSolver.cpp.o" \
"CMakeFiles/hf.dir/RHF.cpp.o" \
"CMakeFiles/hf.dir/RHFC.cpp.o"

# External object files for target hf
hf_EXTERNAL_OBJECTS =

src/libhf.a: src/CMakeFiles/hf.dir/BaseSCFSolver.cpp.o
src/libhf.a: src/CMakeFiles/hf.dir/DIISSCFSolver.cpp.o
src/libhf.a: src/CMakeFiles/hf.dir/PlainSCFSolver.cpp.o
src/libhf.a: src/CMakeFiles/hf.dir/RHF.cpp.o
src/libhf.a: src/CMakeFiles/hf.dir/RHFC.cpp.o
src/libhf.a: src/CMakeFiles/hf.dir/build.make
src/libhf.a: src/CMakeFiles/hf.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/in/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX static library libhf.a"
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/in/src && $(CMAKE_COMMAND) -P CMakeFiles/hf.dir/cmake_clean_target.cmake
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/in/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/hf.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/hf.dir/build: src/libhf.a

.PHONY : src/CMakeFiles/hf.dir/build

src/CMakeFiles/hf.dir/clean:
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/in/src && $(CMAKE_COMMAND) -P CMakeFiles/hf.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/hf.dir/clean

src/CMakeFiles/hf.dir/depend:
	cd /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/in && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/src /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/in /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/in/src /Users/wulfix/Desktop/Cursussen_Gent/ThesisDir/Libraries/GQCG/hf/in/src/CMakeFiles/hf.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/hf.dir/depend
