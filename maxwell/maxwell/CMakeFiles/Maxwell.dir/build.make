# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/rhatch/OOPS

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/rhatch/OOPS/maxwell

# Include any dependencies generated for this target.
include maxwell/CMakeFiles/Maxwell.dir/depend.make

# Include the progress variables for this target.
include maxwell/CMakeFiles/Maxwell.dir/progress.make

# Include the compile flags for this target's objects.
include maxwell/CMakeFiles/Maxwell.dir/flags.make

maxwell/CMakeFiles/Maxwell.dir/src/maxwell.cpp.o: maxwell/CMakeFiles/Maxwell.dir/flags.make
maxwell/CMakeFiles/Maxwell.dir/src/maxwell.cpp.o: src/maxwell.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/rhatch/OOPS/maxwell/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object maxwell/CMakeFiles/Maxwell.dir/src/maxwell.cpp.o"
	cd /home/rhatch/OOPS/maxwell/maxwell && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Maxwell.dir/src/maxwell.cpp.o -c /home/rhatch/OOPS/maxwell/src/maxwell.cpp

maxwell/CMakeFiles/Maxwell.dir/src/maxwell.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Maxwell.dir/src/maxwell.cpp.i"
	cd /home/rhatch/OOPS/maxwell/maxwell && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/rhatch/OOPS/maxwell/src/maxwell.cpp > CMakeFiles/Maxwell.dir/src/maxwell.cpp.i

maxwell/CMakeFiles/Maxwell.dir/src/maxwell.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Maxwell.dir/src/maxwell.cpp.s"
	cd /home/rhatch/OOPS/maxwell/maxwell && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/rhatch/OOPS/maxwell/src/maxwell.cpp -o CMakeFiles/Maxwell.dir/src/maxwell.cpp.s

maxwell/CMakeFiles/Maxwell.dir/src/maxwell.cpp.o.requires:

.PHONY : maxwell/CMakeFiles/Maxwell.dir/src/maxwell.cpp.o.requires

maxwell/CMakeFiles/Maxwell.dir/src/maxwell.cpp.o.provides: maxwell/CMakeFiles/Maxwell.dir/src/maxwell.cpp.o.requires
	$(MAKE) -f maxwell/CMakeFiles/Maxwell.dir/build.make maxwell/CMakeFiles/Maxwell.dir/src/maxwell.cpp.o.provides.build
.PHONY : maxwell/CMakeFiles/Maxwell.dir/src/maxwell.cpp.o.provides

maxwell/CMakeFiles/Maxwell.dir/src/maxwell.cpp.o.provides.build: maxwell/CMakeFiles/Maxwell.dir/src/maxwell.cpp.o


maxwell/CMakeFiles/Maxwell.dir/src/main.cpp.o: maxwell/CMakeFiles/Maxwell.dir/flags.make
maxwell/CMakeFiles/Maxwell.dir/src/main.cpp.o: src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/rhatch/OOPS/maxwell/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object maxwell/CMakeFiles/Maxwell.dir/src/main.cpp.o"
	cd /home/rhatch/OOPS/maxwell/maxwell && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Maxwell.dir/src/main.cpp.o -c /home/rhatch/OOPS/maxwell/src/main.cpp

maxwell/CMakeFiles/Maxwell.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Maxwell.dir/src/main.cpp.i"
	cd /home/rhatch/OOPS/maxwell/maxwell && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/rhatch/OOPS/maxwell/src/main.cpp > CMakeFiles/Maxwell.dir/src/main.cpp.i

maxwell/CMakeFiles/Maxwell.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Maxwell.dir/src/main.cpp.s"
	cd /home/rhatch/OOPS/maxwell/maxwell && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/rhatch/OOPS/maxwell/src/main.cpp -o CMakeFiles/Maxwell.dir/src/main.cpp.s

maxwell/CMakeFiles/Maxwell.dir/src/main.cpp.o.requires:

.PHONY : maxwell/CMakeFiles/Maxwell.dir/src/main.cpp.o.requires

maxwell/CMakeFiles/Maxwell.dir/src/main.cpp.o.provides: maxwell/CMakeFiles/Maxwell.dir/src/main.cpp.o.requires
	$(MAKE) -f maxwell/CMakeFiles/Maxwell.dir/build.make maxwell/CMakeFiles/Maxwell.dir/src/main.cpp.o.provides.build
.PHONY : maxwell/CMakeFiles/Maxwell.dir/src/main.cpp.o.provides

maxwell/CMakeFiles/Maxwell.dir/src/main.cpp.o.provides.build: maxwell/CMakeFiles/Maxwell.dir/src/main.cpp.o


# Object files for target Maxwell
Maxwell_OBJECTS = \
"CMakeFiles/Maxwell.dir/src/maxwell.cpp.o" \
"CMakeFiles/Maxwell.dir/src/main.cpp.o"

# External object files for target Maxwell
Maxwell_EXTERNAL_OBJECTS =

maxwell/Maxwell: maxwell/CMakeFiles/Maxwell.dir/src/maxwell.cpp.o
maxwell/Maxwell: maxwell/CMakeFiles/Maxwell.dir/src/main.cpp.o
maxwell/Maxwell: maxwell/CMakeFiles/Maxwell.dir/build.make
maxwell/Maxwell: liboops.a
maxwell/Maxwell: maxwell/CMakeFiles/Maxwell.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/rhatch/OOPS/maxwell/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable Maxwell"
	cd /home/rhatch/OOPS/maxwell/maxwell && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Maxwell.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
maxwell/CMakeFiles/Maxwell.dir/build: maxwell/Maxwell

.PHONY : maxwell/CMakeFiles/Maxwell.dir/build

maxwell/CMakeFiles/Maxwell.dir/requires: maxwell/CMakeFiles/Maxwell.dir/src/maxwell.cpp.o.requires
maxwell/CMakeFiles/Maxwell.dir/requires: maxwell/CMakeFiles/Maxwell.dir/src/main.cpp.o.requires

.PHONY : maxwell/CMakeFiles/Maxwell.dir/requires

maxwell/CMakeFiles/Maxwell.dir/clean:
	cd /home/rhatch/OOPS/maxwell/maxwell && $(CMAKE_COMMAND) -P CMakeFiles/Maxwell.dir/cmake_clean.cmake
.PHONY : maxwell/CMakeFiles/Maxwell.dir/clean

maxwell/CMakeFiles/Maxwell.dir/depend:
	cd /home/rhatch/OOPS/maxwell && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/rhatch/OOPS /home/rhatch/OOPS/maxwell /home/rhatch/OOPS/maxwell /home/rhatch/OOPS/maxwell/maxwell /home/rhatch/OOPS/maxwell/maxwell/CMakeFiles/Maxwell.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : maxwell/CMakeFiles/Maxwell.dir/depend

