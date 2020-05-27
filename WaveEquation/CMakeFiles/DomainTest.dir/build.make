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
CMAKE_BINARY_DIR = /home/rhatch/OOPS/WaveEquation

# Include any dependencies generated for this target.
include CMakeFiles/DomainTest.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/DomainTest.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/DomainTest.dir/flags.make

CMakeFiles/DomainTest.dir/src/grid.cpp.o: CMakeFiles/DomainTest.dir/flags.make
CMakeFiles/DomainTest.dir/src/grid.cpp.o: ../src/grid.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/rhatch/OOPS/WaveEquation/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/DomainTest.dir/src/grid.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/DomainTest.dir/src/grid.cpp.o -c /home/rhatch/OOPS/src/grid.cpp

CMakeFiles/DomainTest.dir/src/grid.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DomainTest.dir/src/grid.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/rhatch/OOPS/src/grid.cpp > CMakeFiles/DomainTest.dir/src/grid.cpp.i

CMakeFiles/DomainTest.dir/src/grid.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DomainTest.dir/src/grid.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/rhatch/OOPS/src/grid.cpp -o CMakeFiles/DomainTest.dir/src/grid.cpp.s

CMakeFiles/DomainTest.dir/src/grid.cpp.o.requires:

.PHONY : CMakeFiles/DomainTest.dir/src/grid.cpp.o.requires

CMakeFiles/DomainTest.dir/src/grid.cpp.o.provides: CMakeFiles/DomainTest.dir/src/grid.cpp.o.requires
	$(MAKE) -f CMakeFiles/DomainTest.dir/build.make CMakeFiles/DomainTest.dir/src/grid.cpp.o.provides.build
.PHONY : CMakeFiles/DomainTest.dir/src/grid.cpp.o.provides

CMakeFiles/DomainTest.dir/src/grid.cpp.o.provides.build: CMakeFiles/DomainTest.dir/src/grid.cpp.o


CMakeFiles/DomainTest.dir/src/domain.cpp.o: CMakeFiles/DomainTest.dir/flags.make
CMakeFiles/DomainTest.dir/src/domain.cpp.o: ../src/domain.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/rhatch/OOPS/WaveEquation/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/DomainTest.dir/src/domain.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/DomainTest.dir/src/domain.cpp.o -c /home/rhatch/OOPS/src/domain.cpp

CMakeFiles/DomainTest.dir/src/domain.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DomainTest.dir/src/domain.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/rhatch/OOPS/src/domain.cpp > CMakeFiles/DomainTest.dir/src/domain.cpp.i

CMakeFiles/DomainTest.dir/src/domain.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DomainTest.dir/src/domain.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/rhatch/OOPS/src/domain.cpp -o CMakeFiles/DomainTest.dir/src/domain.cpp.s

CMakeFiles/DomainTest.dir/src/domain.cpp.o.requires:

.PHONY : CMakeFiles/DomainTest.dir/src/domain.cpp.o.requires

CMakeFiles/DomainTest.dir/src/domain.cpp.o.provides: CMakeFiles/DomainTest.dir/src/domain.cpp.o.requires
	$(MAKE) -f CMakeFiles/DomainTest.dir/build.make CMakeFiles/DomainTest.dir/src/domain.cpp.o.provides.build
.PHONY : CMakeFiles/DomainTest.dir/src/domain.cpp.o.provides

CMakeFiles/DomainTest.dir/src/domain.cpp.o.provides.build: CMakeFiles/DomainTest.dir/src/domain.cpp.o


CMakeFiles/DomainTest.dir/src/rk4.cpp.o: CMakeFiles/DomainTest.dir/flags.make
CMakeFiles/DomainTest.dir/src/rk4.cpp.o: ../src/rk4.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/rhatch/OOPS/WaveEquation/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/DomainTest.dir/src/rk4.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/DomainTest.dir/src/rk4.cpp.o -c /home/rhatch/OOPS/src/rk4.cpp

CMakeFiles/DomainTest.dir/src/rk4.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DomainTest.dir/src/rk4.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/rhatch/OOPS/src/rk4.cpp > CMakeFiles/DomainTest.dir/src/rk4.cpp.i

CMakeFiles/DomainTest.dir/src/rk4.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DomainTest.dir/src/rk4.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/rhatch/OOPS/src/rk4.cpp -o CMakeFiles/DomainTest.dir/src/rk4.cpp.s

CMakeFiles/DomainTest.dir/src/rk4.cpp.o.requires:

.PHONY : CMakeFiles/DomainTest.dir/src/rk4.cpp.o.requires

CMakeFiles/DomainTest.dir/src/rk4.cpp.o.provides: CMakeFiles/DomainTest.dir/src/rk4.cpp.o.requires
	$(MAKE) -f CMakeFiles/DomainTest.dir/build.make CMakeFiles/DomainTest.dir/src/rk4.cpp.o.provides.build
.PHONY : CMakeFiles/DomainTest.dir/src/rk4.cpp.o.provides

CMakeFiles/DomainTest.dir/src/rk4.cpp.o.provides.build: CMakeFiles/DomainTest.dir/src/rk4.cpp.o


CMakeFiles/DomainTest.dir/src/ode.cpp.o: CMakeFiles/DomainTest.dir/flags.make
CMakeFiles/DomainTest.dir/src/ode.cpp.o: ../src/ode.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/rhatch/OOPS/WaveEquation/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/DomainTest.dir/src/ode.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/DomainTest.dir/src/ode.cpp.o -c /home/rhatch/OOPS/src/ode.cpp

CMakeFiles/DomainTest.dir/src/ode.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DomainTest.dir/src/ode.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/rhatch/OOPS/src/ode.cpp > CMakeFiles/DomainTest.dir/src/ode.cpp.i

CMakeFiles/DomainTest.dir/src/ode.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DomainTest.dir/src/ode.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/rhatch/OOPS/src/ode.cpp -o CMakeFiles/DomainTest.dir/src/ode.cpp.s

CMakeFiles/DomainTest.dir/src/ode.cpp.o.requires:

.PHONY : CMakeFiles/DomainTest.dir/src/ode.cpp.o.requires

CMakeFiles/DomainTest.dir/src/ode.cpp.o.provides: CMakeFiles/DomainTest.dir/src/ode.cpp.o.requires
	$(MAKE) -f CMakeFiles/DomainTest.dir/build.make CMakeFiles/DomainTest.dir/src/ode.cpp.o.provides.build
.PHONY : CMakeFiles/DomainTest.dir/src/ode.cpp.o.provides

CMakeFiles/DomainTest.dir/src/ode.cpp.o.provides.build: CMakeFiles/DomainTest.dir/src/ode.cpp.o


CMakeFiles/DomainTest.dir/src/solverdata.cpp.o: CMakeFiles/DomainTest.dir/flags.make
CMakeFiles/DomainTest.dir/src/solverdata.cpp.o: ../src/solverdata.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/rhatch/OOPS/WaveEquation/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/DomainTest.dir/src/solverdata.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/DomainTest.dir/src/solverdata.cpp.o -c /home/rhatch/OOPS/src/solverdata.cpp

CMakeFiles/DomainTest.dir/src/solverdata.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DomainTest.dir/src/solverdata.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/rhatch/OOPS/src/solverdata.cpp > CMakeFiles/DomainTest.dir/src/solverdata.cpp.i

CMakeFiles/DomainTest.dir/src/solverdata.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DomainTest.dir/src/solverdata.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/rhatch/OOPS/src/solverdata.cpp -o CMakeFiles/DomainTest.dir/src/solverdata.cpp.s

CMakeFiles/DomainTest.dir/src/solverdata.cpp.o.requires:

.PHONY : CMakeFiles/DomainTest.dir/src/solverdata.cpp.o.requires

CMakeFiles/DomainTest.dir/src/solverdata.cpp.o.provides: CMakeFiles/DomainTest.dir/src/solverdata.cpp.o.requires
	$(MAKE) -f CMakeFiles/DomainTest.dir/build.make CMakeFiles/DomainTest.dir/src/solverdata.cpp.o.provides.build
.PHONY : CMakeFiles/DomainTest.dir/src/solverdata.cpp.o.provides

CMakeFiles/DomainTest.dir/src/solverdata.cpp.o.provides.build: CMakeFiles/DomainTest.dir/src/solverdata.cpp.o


CMakeFiles/DomainTest.dir/src/cubic.cpp.o: CMakeFiles/DomainTest.dir/flags.make
CMakeFiles/DomainTest.dir/src/cubic.cpp.o: ../src/cubic.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/rhatch/OOPS/WaveEquation/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/DomainTest.dir/src/cubic.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/DomainTest.dir/src/cubic.cpp.o -c /home/rhatch/OOPS/src/cubic.cpp

CMakeFiles/DomainTest.dir/src/cubic.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DomainTest.dir/src/cubic.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/rhatch/OOPS/src/cubic.cpp > CMakeFiles/DomainTest.dir/src/cubic.cpp.i

CMakeFiles/DomainTest.dir/src/cubic.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DomainTest.dir/src/cubic.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/rhatch/OOPS/src/cubic.cpp -o CMakeFiles/DomainTest.dir/src/cubic.cpp.s

CMakeFiles/DomainTest.dir/src/cubic.cpp.o.requires:

.PHONY : CMakeFiles/DomainTest.dir/src/cubic.cpp.o.requires

CMakeFiles/DomainTest.dir/src/cubic.cpp.o.provides: CMakeFiles/DomainTest.dir/src/cubic.cpp.o.requires
	$(MAKE) -f CMakeFiles/DomainTest.dir/build.make CMakeFiles/DomainTest.dir/src/cubic.cpp.o.provides.build
.PHONY : CMakeFiles/DomainTest.dir/src/cubic.cpp.o.provides

CMakeFiles/DomainTest.dir/src/cubic.cpp.o.provides.build: CMakeFiles/DomainTest.dir/src/cubic.cpp.o


CMakeFiles/DomainTest.dir/src/interpolator.cpp.o: CMakeFiles/DomainTest.dir/flags.make
CMakeFiles/DomainTest.dir/src/interpolator.cpp.o: ../src/interpolator.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/rhatch/OOPS/WaveEquation/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/DomainTest.dir/src/interpolator.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/DomainTest.dir/src/interpolator.cpp.o -c /home/rhatch/OOPS/src/interpolator.cpp

CMakeFiles/DomainTest.dir/src/interpolator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DomainTest.dir/src/interpolator.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/rhatch/OOPS/src/interpolator.cpp > CMakeFiles/DomainTest.dir/src/interpolator.cpp.i

CMakeFiles/DomainTest.dir/src/interpolator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DomainTest.dir/src/interpolator.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/rhatch/OOPS/src/interpolator.cpp -o CMakeFiles/DomainTest.dir/src/interpolator.cpp.s

CMakeFiles/DomainTest.dir/src/interpolator.cpp.o.requires:

.PHONY : CMakeFiles/DomainTest.dir/src/interpolator.cpp.o.requires

CMakeFiles/DomainTest.dir/src/interpolator.cpp.o.provides: CMakeFiles/DomainTest.dir/src/interpolator.cpp.o.requires
	$(MAKE) -f CMakeFiles/DomainTest.dir/build.make CMakeFiles/DomainTest.dir/src/interpolator.cpp.o.provides.build
.PHONY : CMakeFiles/DomainTest.dir/src/interpolator.cpp.o.provides

CMakeFiles/DomainTest.dir/src/interpolator.cpp.o.provides.build: CMakeFiles/DomainTest.dir/src/interpolator.cpp.o


CMakeFiles/DomainTest.dir/src/cubicinterpolator.cpp.o: CMakeFiles/DomainTest.dir/flags.make
CMakeFiles/DomainTest.dir/src/cubicinterpolator.cpp.o: ../src/cubicinterpolator.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/rhatch/OOPS/WaveEquation/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/DomainTest.dir/src/cubicinterpolator.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/DomainTest.dir/src/cubicinterpolator.cpp.o -c /home/rhatch/OOPS/src/cubicinterpolator.cpp

CMakeFiles/DomainTest.dir/src/cubicinterpolator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DomainTest.dir/src/cubicinterpolator.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/rhatch/OOPS/src/cubicinterpolator.cpp > CMakeFiles/DomainTest.dir/src/cubicinterpolator.cpp.i

CMakeFiles/DomainTest.dir/src/cubicinterpolator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DomainTest.dir/src/cubicinterpolator.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/rhatch/OOPS/src/cubicinterpolator.cpp -o CMakeFiles/DomainTest.dir/src/cubicinterpolator.cpp.s

CMakeFiles/DomainTest.dir/src/cubicinterpolator.cpp.o.requires:

.PHONY : CMakeFiles/DomainTest.dir/src/cubicinterpolator.cpp.o.requires

CMakeFiles/DomainTest.dir/src/cubicinterpolator.cpp.o.provides: CMakeFiles/DomainTest.dir/src/cubicinterpolator.cpp.o.requires
	$(MAKE) -f CMakeFiles/DomainTest.dir/build.make CMakeFiles/DomainTest.dir/src/cubicinterpolator.cpp.o.provides.build
.PHONY : CMakeFiles/DomainTest.dir/src/cubicinterpolator.cpp.o.provides

CMakeFiles/DomainTest.dir/src/cubicinterpolator.cpp.o.provides.build: CMakeFiles/DomainTest.dir/src/cubicinterpolator.cpp.o


CMakeFiles/DomainTest.dir/src/polynomialinterpolator.cpp.o: CMakeFiles/DomainTest.dir/flags.make
CMakeFiles/DomainTest.dir/src/polynomialinterpolator.cpp.o: ../src/polynomialinterpolator.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/rhatch/OOPS/WaveEquation/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/DomainTest.dir/src/polynomialinterpolator.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/DomainTest.dir/src/polynomialinterpolator.cpp.o -c /home/rhatch/OOPS/src/polynomialinterpolator.cpp

CMakeFiles/DomainTest.dir/src/polynomialinterpolator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DomainTest.dir/src/polynomialinterpolator.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/rhatch/OOPS/src/polynomialinterpolator.cpp > CMakeFiles/DomainTest.dir/src/polynomialinterpolator.cpp.i

CMakeFiles/DomainTest.dir/src/polynomialinterpolator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DomainTest.dir/src/polynomialinterpolator.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/rhatch/OOPS/src/polynomialinterpolator.cpp -o CMakeFiles/DomainTest.dir/src/polynomialinterpolator.cpp.s

CMakeFiles/DomainTest.dir/src/polynomialinterpolator.cpp.o.requires:

.PHONY : CMakeFiles/DomainTest.dir/src/polynomialinterpolator.cpp.o.requires

CMakeFiles/DomainTest.dir/src/polynomialinterpolator.cpp.o.provides: CMakeFiles/DomainTest.dir/src/polynomialinterpolator.cpp.o.requires
	$(MAKE) -f CMakeFiles/DomainTest.dir/build.make CMakeFiles/DomainTest.dir/src/polynomialinterpolator.cpp.o.provides.build
.PHONY : CMakeFiles/DomainTest.dir/src/polynomialinterpolator.cpp.o.provides

CMakeFiles/DomainTest.dir/src/polynomialinterpolator.cpp.o.provides.build: CMakeFiles/DomainTest.dir/src/polynomialinterpolator.cpp.o


CMakeFiles/DomainTest.dir/src/paramreader.cpp.o: CMakeFiles/DomainTest.dir/flags.make
CMakeFiles/DomainTest.dir/src/paramreader.cpp.o: ../src/paramreader.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/rhatch/OOPS/WaveEquation/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/DomainTest.dir/src/paramreader.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/DomainTest.dir/src/paramreader.cpp.o -c /home/rhatch/OOPS/src/paramreader.cpp

CMakeFiles/DomainTest.dir/src/paramreader.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DomainTest.dir/src/paramreader.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/rhatch/OOPS/src/paramreader.cpp > CMakeFiles/DomainTest.dir/src/paramreader.cpp.i

CMakeFiles/DomainTest.dir/src/paramreader.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DomainTest.dir/src/paramreader.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/rhatch/OOPS/src/paramreader.cpp -o CMakeFiles/DomainTest.dir/src/paramreader.cpp.s

CMakeFiles/DomainTest.dir/src/paramreader.cpp.o.requires:

.PHONY : CMakeFiles/DomainTest.dir/src/paramreader.cpp.o.requires

CMakeFiles/DomainTest.dir/src/paramreader.cpp.o.provides: CMakeFiles/DomainTest.dir/src/paramreader.cpp.o.requires
	$(MAKE) -f CMakeFiles/DomainTest.dir/build.make CMakeFiles/DomainTest.dir/src/paramreader.cpp.o.provides.build
.PHONY : CMakeFiles/DomainTest.dir/src/paramreader.cpp.o.provides

CMakeFiles/DomainTest.dir/src/paramreader.cpp.o.provides.build: CMakeFiles/DomainTest.dir/src/paramreader.cpp.o


CMakeFiles/DomainTest.dir/tests/domainTest.cpp.o: CMakeFiles/DomainTest.dir/flags.make
CMakeFiles/DomainTest.dir/tests/domainTest.cpp.o: ../tests/domainTest.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/rhatch/OOPS/WaveEquation/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object CMakeFiles/DomainTest.dir/tests/domainTest.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/DomainTest.dir/tests/domainTest.cpp.o -c /home/rhatch/OOPS/tests/domainTest.cpp

CMakeFiles/DomainTest.dir/tests/domainTest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DomainTest.dir/tests/domainTest.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/rhatch/OOPS/tests/domainTest.cpp > CMakeFiles/DomainTest.dir/tests/domainTest.cpp.i

CMakeFiles/DomainTest.dir/tests/domainTest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DomainTest.dir/tests/domainTest.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/rhatch/OOPS/tests/domainTest.cpp -o CMakeFiles/DomainTest.dir/tests/domainTest.cpp.s

CMakeFiles/DomainTest.dir/tests/domainTest.cpp.o.requires:

.PHONY : CMakeFiles/DomainTest.dir/tests/domainTest.cpp.o.requires

CMakeFiles/DomainTest.dir/tests/domainTest.cpp.o.provides: CMakeFiles/DomainTest.dir/tests/domainTest.cpp.o.requires
	$(MAKE) -f CMakeFiles/DomainTest.dir/build.make CMakeFiles/DomainTest.dir/tests/domainTest.cpp.o.provides.build
.PHONY : CMakeFiles/DomainTest.dir/tests/domainTest.cpp.o.provides

CMakeFiles/DomainTest.dir/tests/domainTest.cpp.o.provides.build: CMakeFiles/DomainTest.dir/tests/domainTest.cpp.o


# Object files for target DomainTest
DomainTest_OBJECTS = \
"CMakeFiles/DomainTest.dir/src/grid.cpp.o" \
"CMakeFiles/DomainTest.dir/src/domain.cpp.o" \
"CMakeFiles/DomainTest.dir/src/rk4.cpp.o" \
"CMakeFiles/DomainTest.dir/src/ode.cpp.o" \
"CMakeFiles/DomainTest.dir/src/solverdata.cpp.o" \
"CMakeFiles/DomainTest.dir/src/cubic.cpp.o" \
"CMakeFiles/DomainTest.dir/src/interpolator.cpp.o" \
"CMakeFiles/DomainTest.dir/src/cubicinterpolator.cpp.o" \
"CMakeFiles/DomainTest.dir/src/polynomialinterpolator.cpp.o" \
"CMakeFiles/DomainTest.dir/src/paramreader.cpp.o" \
"CMakeFiles/DomainTest.dir/tests/domainTest.cpp.o"

# External object files for target DomainTest
DomainTest_EXTERNAL_OBJECTS =

DomainTest: CMakeFiles/DomainTest.dir/src/grid.cpp.o
DomainTest: CMakeFiles/DomainTest.dir/src/domain.cpp.o
DomainTest: CMakeFiles/DomainTest.dir/src/rk4.cpp.o
DomainTest: CMakeFiles/DomainTest.dir/src/ode.cpp.o
DomainTest: CMakeFiles/DomainTest.dir/src/solverdata.cpp.o
DomainTest: CMakeFiles/DomainTest.dir/src/cubic.cpp.o
DomainTest: CMakeFiles/DomainTest.dir/src/interpolator.cpp.o
DomainTest: CMakeFiles/DomainTest.dir/src/cubicinterpolator.cpp.o
DomainTest: CMakeFiles/DomainTest.dir/src/polynomialinterpolator.cpp.o
DomainTest: CMakeFiles/DomainTest.dir/src/paramreader.cpp.o
DomainTest: CMakeFiles/DomainTest.dir/tests/domainTest.cpp.o
DomainTest: CMakeFiles/DomainTest.dir/build.make
DomainTest: CMakeFiles/DomainTest.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/rhatch/OOPS/WaveEquation/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Linking CXX executable DomainTest"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/DomainTest.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/DomainTest.dir/build: DomainTest

.PHONY : CMakeFiles/DomainTest.dir/build

CMakeFiles/DomainTest.dir/requires: CMakeFiles/DomainTest.dir/src/grid.cpp.o.requires
CMakeFiles/DomainTest.dir/requires: CMakeFiles/DomainTest.dir/src/domain.cpp.o.requires
CMakeFiles/DomainTest.dir/requires: CMakeFiles/DomainTest.dir/src/rk4.cpp.o.requires
CMakeFiles/DomainTest.dir/requires: CMakeFiles/DomainTest.dir/src/ode.cpp.o.requires
CMakeFiles/DomainTest.dir/requires: CMakeFiles/DomainTest.dir/src/solverdata.cpp.o.requires
CMakeFiles/DomainTest.dir/requires: CMakeFiles/DomainTest.dir/src/cubic.cpp.o.requires
CMakeFiles/DomainTest.dir/requires: CMakeFiles/DomainTest.dir/src/interpolator.cpp.o.requires
CMakeFiles/DomainTest.dir/requires: CMakeFiles/DomainTest.dir/src/cubicinterpolator.cpp.o.requires
CMakeFiles/DomainTest.dir/requires: CMakeFiles/DomainTest.dir/src/polynomialinterpolator.cpp.o.requires
CMakeFiles/DomainTest.dir/requires: CMakeFiles/DomainTest.dir/src/paramreader.cpp.o.requires
CMakeFiles/DomainTest.dir/requires: CMakeFiles/DomainTest.dir/tests/domainTest.cpp.o.requires

.PHONY : CMakeFiles/DomainTest.dir/requires

CMakeFiles/DomainTest.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/DomainTest.dir/cmake_clean.cmake
.PHONY : CMakeFiles/DomainTest.dir/clean

CMakeFiles/DomainTest.dir/depend:
	cd /home/rhatch/OOPS/WaveEquation && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/rhatch/OOPS /home/rhatch/OOPS /home/rhatch/OOPS/WaveEquation /home/rhatch/OOPS/WaveEquation /home/rhatch/OOPS/WaveEquation/CMakeFiles/DomainTest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/DomainTest.dir/depend

