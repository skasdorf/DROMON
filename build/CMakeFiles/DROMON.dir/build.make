# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/skasdorf/Documents/DROMON-main

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/skasdorf/Documents/DROMON-main/build

# Include any dependencies generated for this target.
include CMakeFiles/DROMON.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/DROMON.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/DROMON.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/DROMON.dir/flags.make

CMakeFiles/DROMON.dir/source/mesh.cpp.o: CMakeFiles/DROMON.dir/flags.make
CMakeFiles/DROMON.dir/source/mesh.cpp.o: ../source/mesh.cpp
CMakeFiles/DROMON.dir/source/mesh.cpp.o: CMakeFiles/DROMON.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/skasdorf/Documents/DROMON-main/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/DROMON.dir/source/mesh.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/DROMON.dir/source/mesh.cpp.o -MF CMakeFiles/DROMON.dir/source/mesh.cpp.o.d -o CMakeFiles/DROMON.dir/source/mesh.cpp.o -c /home/skasdorf/Documents/DROMON-main/source/mesh.cpp

CMakeFiles/DROMON.dir/source/mesh.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DROMON.dir/source/mesh.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/skasdorf/Documents/DROMON-main/source/mesh.cpp > CMakeFiles/DROMON.dir/source/mesh.cpp.i

CMakeFiles/DROMON.dir/source/mesh.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DROMON.dir/source/mesh.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/skasdorf/Documents/DROMON-main/source/mesh.cpp -o CMakeFiles/DROMON.dir/source/mesh.cpp.s

CMakeFiles/DROMON.dir/source/GeomBase.cpp.o: CMakeFiles/DROMON.dir/flags.make
CMakeFiles/DROMON.dir/source/GeomBase.cpp.o: ../source/GeomBase.cpp
CMakeFiles/DROMON.dir/source/GeomBase.cpp.o: CMakeFiles/DROMON.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/skasdorf/Documents/DROMON-main/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/DROMON.dir/source/GeomBase.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/DROMON.dir/source/GeomBase.cpp.o -MF CMakeFiles/DROMON.dir/source/GeomBase.cpp.o.d -o CMakeFiles/DROMON.dir/source/GeomBase.cpp.o -c /home/skasdorf/Documents/DROMON-main/source/GeomBase.cpp

CMakeFiles/DROMON.dir/source/GeomBase.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DROMON.dir/source/GeomBase.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/skasdorf/Documents/DROMON-main/source/GeomBase.cpp > CMakeFiles/DROMON.dir/source/GeomBase.cpp.i

CMakeFiles/DROMON.dir/source/GeomBase.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DROMON.dir/source/GeomBase.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/skasdorf/Documents/DROMON-main/source/GeomBase.cpp -o CMakeFiles/DROMON.dir/source/GeomBase.cpp.s

CMakeFiles/DROMON.dir/source/MeshGenerator.cpp.o: CMakeFiles/DROMON.dir/flags.make
CMakeFiles/DROMON.dir/source/MeshGenerator.cpp.o: ../source/MeshGenerator.cpp
CMakeFiles/DROMON.dir/source/MeshGenerator.cpp.o: CMakeFiles/DROMON.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/skasdorf/Documents/DROMON-main/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/DROMON.dir/source/MeshGenerator.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/DROMON.dir/source/MeshGenerator.cpp.o -MF CMakeFiles/DROMON.dir/source/MeshGenerator.cpp.o.d -o CMakeFiles/DROMON.dir/source/MeshGenerator.cpp.o -c /home/skasdorf/Documents/DROMON-main/source/MeshGenerator.cpp

CMakeFiles/DROMON.dir/source/MeshGenerator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DROMON.dir/source/MeshGenerator.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/skasdorf/Documents/DROMON-main/source/MeshGenerator.cpp > CMakeFiles/DROMON.dir/source/MeshGenerator.cpp.i

CMakeFiles/DROMON.dir/source/MeshGenerator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DROMON.dir/source/MeshGenerator.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/skasdorf/Documents/DROMON-main/source/MeshGenerator.cpp -o CMakeFiles/DROMON.dir/source/MeshGenerator.cpp.s

CMakeFiles/DROMON.dir/source/Point.cpp.o: CMakeFiles/DROMON.dir/flags.make
CMakeFiles/DROMON.dir/source/Point.cpp.o: ../source/Point.cpp
CMakeFiles/DROMON.dir/source/Point.cpp.o: CMakeFiles/DROMON.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/skasdorf/Documents/DROMON-main/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/DROMON.dir/source/Point.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/DROMON.dir/source/Point.cpp.o -MF CMakeFiles/DROMON.dir/source/Point.cpp.o.d -o CMakeFiles/DROMON.dir/source/Point.cpp.o -c /home/skasdorf/Documents/DROMON-main/source/Point.cpp

CMakeFiles/DROMON.dir/source/Point.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DROMON.dir/source/Point.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/skasdorf/Documents/DROMON-main/source/Point.cpp > CMakeFiles/DROMON.dir/source/Point.cpp.i

CMakeFiles/DROMON.dir/source/Point.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DROMON.dir/source/Point.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/skasdorf/Documents/DROMON-main/source/Point.cpp -o CMakeFiles/DROMON.dir/source/Point.cpp.s

CMakeFiles/DROMON.dir/source/DataOut.cpp.o: CMakeFiles/DROMON.dir/flags.make
CMakeFiles/DROMON.dir/source/DataOut.cpp.o: ../source/DataOut.cpp
CMakeFiles/DROMON.dir/source/DataOut.cpp.o: CMakeFiles/DROMON.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/skasdorf/Documents/DROMON-main/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/DROMON.dir/source/DataOut.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/DROMON.dir/source/DataOut.cpp.o -MF CMakeFiles/DROMON.dir/source/DataOut.cpp.o.d -o CMakeFiles/DROMON.dir/source/DataOut.cpp.o -c /home/skasdorf/Documents/DROMON-main/source/DataOut.cpp

CMakeFiles/DROMON.dir/source/DataOut.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DROMON.dir/source/DataOut.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/skasdorf/Documents/DROMON-main/source/DataOut.cpp > CMakeFiles/DROMON.dir/source/DataOut.cpp.i

CMakeFiles/DROMON.dir/source/DataOut.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DROMON.dir/source/DataOut.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/skasdorf/Documents/DROMON-main/source/DataOut.cpp -o CMakeFiles/DROMON.dir/source/DataOut.cpp.s

CMakeFiles/DROMON.dir/source/MeshBase.cpp.o: CMakeFiles/DROMON.dir/flags.make
CMakeFiles/DROMON.dir/source/MeshBase.cpp.o: ../source/MeshBase.cpp
CMakeFiles/DROMON.dir/source/MeshBase.cpp.o: CMakeFiles/DROMON.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/skasdorf/Documents/DROMON-main/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/DROMON.dir/source/MeshBase.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/DROMON.dir/source/MeshBase.cpp.o -MF CMakeFiles/DROMON.dir/source/MeshBase.cpp.o.d -o CMakeFiles/DROMON.dir/source/MeshBase.cpp.o -c /home/skasdorf/Documents/DROMON-main/source/MeshBase.cpp

CMakeFiles/DROMON.dir/source/MeshBase.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DROMON.dir/source/MeshBase.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/skasdorf/Documents/DROMON-main/source/MeshBase.cpp > CMakeFiles/DROMON.dir/source/MeshBase.cpp.i

CMakeFiles/DROMON.dir/source/MeshBase.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DROMON.dir/source/MeshBase.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/skasdorf/Documents/DROMON-main/source/MeshBase.cpp -o CMakeFiles/DROMON.dir/source/MeshBase.cpp.s

CMakeFiles/DROMON.dir/source/FEBase.cpp.o: CMakeFiles/DROMON.dir/flags.make
CMakeFiles/DROMON.dir/source/FEBase.cpp.o: ../source/FEBase.cpp
CMakeFiles/DROMON.dir/source/FEBase.cpp.o: CMakeFiles/DROMON.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/skasdorf/Documents/DROMON-main/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/DROMON.dir/source/FEBase.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/DROMON.dir/source/FEBase.cpp.o -MF CMakeFiles/DROMON.dir/source/FEBase.cpp.o.d -o CMakeFiles/DROMON.dir/source/FEBase.cpp.o -c /home/skasdorf/Documents/DROMON-main/source/FEBase.cpp

CMakeFiles/DROMON.dir/source/FEBase.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DROMON.dir/source/FEBase.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/skasdorf/Documents/DROMON-main/source/FEBase.cpp > CMakeFiles/DROMON.dir/source/FEBase.cpp.i

CMakeFiles/DROMON.dir/source/FEBase.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DROMON.dir/source/FEBase.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/skasdorf/Documents/DROMON-main/source/FEBase.cpp -o CMakeFiles/DROMON.dir/source/FEBase.cpp.s

CMakeFiles/DROMON.dir/source/FE_HdivMaxOrtho.cpp.o: CMakeFiles/DROMON.dir/flags.make
CMakeFiles/DROMON.dir/source/FE_HdivMaxOrtho.cpp.o: ../source/FE_HdivMaxOrtho.cpp
CMakeFiles/DROMON.dir/source/FE_HdivMaxOrtho.cpp.o: CMakeFiles/DROMON.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/skasdorf/Documents/DROMON-main/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/DROMON.dir/source/FE_HdivMaxOrtho.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/DROMON.dir/source/FE_HdivMaxOrtho.cpp.o -MF CMakeFiles/DROMON.dir/source/FE_HdivMaxOrtho.cpp.o.d -o CMakeFiles/DROMON.dir/source/FE_HdivMaxOrtho.cpp.o -c /home/skasdorf/Documents/DROMON-main/source/FE_HdivMaxOrtho.cpp

CMakeFiles/DROMON.dir/source/FE_HdivMaxOrtho.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DROMON.dir/source/FE_HdivMaxOrtho.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/skasdorf/Documents/DROMON-main/source/FE_HdivMaxOrtho.cpp > CMakeFiles/DROMON.dir/source/FE_HdivMaxOrtho.cpp.i

CMakeFiles/DROMON.dir/source/FE_HdivMaxOrtho.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DROMON.dir/source/FE_HdivMaxOrtho.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/skasdorf/Documents/DROMON-main/source/FE_HdivMaxOrtho.cpp -o CMakeFiles/DROMON.dir/source/FE_HdivMaxOrtho.cpp.s

CMakeFiles/DROMON.dir/source/MultiIndex.cpp.o: CMakeFiles/DROMON.dir/flags.make
CMakeFiles/DROMON.dir/source/MultiIndex.cpp.o: ../source/MultiIndex.cpp
CMakeFiles/DROMON.dir/source/MultiIndex.cpp.o: CMakeFiles/DROMON.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/skasdorf/Documents/DROMON-main/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/DROMON.dir/source/MultiIndex.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/DROMON.dir/source/MultiIndex.cpp.o -MF CMakeFiles/DROMON.dir/source/MultiIndex.cpp.o.d -o CMakeFiles/DROMON.dir/source/MultiIndex.cpp.o -c /home/skasdorf/Documents/DROMON-main/source/MultiIndex.cpp

CMakeFiles/DROMON.dir/source/MultiIndex.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DROMON.dir/source/MultiIndex.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/skasdorf/Documents/DROMON-main/source/MultiIndex.cpp > CMakeFiles/DROMON.dir/source/MultiIndex.cpp.i

CMakeFiles/DROMON.dir/source/MultiIndex.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DROMON.dir/source/MultiIndex.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/skasdorf/Documents/DROMON-main/source/MultiIndex.cpp -o CMakeFiles/DROMON.dir/source/MultiIndex.cpp.s

CMakeFiles/DROMON.dir/source/FECollection.cpp.o: CMakeFiles/DROMON.dir/flags.make
CMakeFiles/DROMON.dir/source/FECollection.cpp.o: ../source/FECollection.cpp
CMakeFiles/DROMON.dir/source/FECollection.cpp.o: CMakeFiles/DROMON.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/skasdorf/Documents/DROMON-main/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/DROMON.dir/source/FECollection.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/DROMON.dir/source/FECollection.cpp.o -MF CMakeFiles/DROMON.dir/source/FECollection.cpp.o.d -o CMakeFiles/DROMON.dir/source/FECollection.cpp.o -c /home/skasdorf/Documents/DROMON-main/source/FECollection.cpp

CMakeFiles/DROMON.dir/source/FECollection.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DROMON.dir/source/FECollection.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/skasdorf/Documents/DROMON-main/source/FECollection.cpp > CMakeFiles/DROMON.dir/source/FECollection.cpp.i

CMakeFiles/DROMON.dir/source/FECollection.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DROMON.dir/source/FECollection.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/skasdorf/Documents/DROMON-main/source/FECollection.cpp -o CMakeFiles/DROMON.dir/source/FECollection.cpp.s

CMakeFiles/DROMON.dir/source/DoFHandler.cpp.o: CMakeFiles/DROMON.dir/flags.make
CMakeFiles/DROMON.dir/source/DoFHandler.cpp.o: ../source/DoFHandler.cpp
CMakeFiles/DROMON.dir/source/DoFHandler.cpp.o: CMakeFiles/DROMON.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/skasdorf/Documents/DROMON-main/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object CMakeFiles/DROMON.dir/source/DoFHandler.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/DROMON.dir/source/DoFHandler.cpp.o -MF CMakeFiles/DROMON.dir/source/DoFHandler.cpp.o.d -o CMakeFiles/DROMON.dir/source/DoFHandler.cpp.o -c /home/skasdorf/Documents/DROMON-main/source/DoFHandler.cpp

CMakeFiles/DROMON.dir/source/DoFHandler.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DROMON.dir/source/DoFHandler.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/skasdorf/Documents/DROMON-main/source/DoFHandler.cpp > CMakeFiles/DROMON.dir/source/DoFHandler.cpp.i

CMakeFiles/DROMON.dir/source/DoFHandler.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DROMON.dir/source/DoFHandler.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/skasdorf/Documents/DROMON-main/source/DoFHandler.cpp -o CMakeFiles/DROMON.dir/source/DoFHandler.cpp.s

CMakeFiles/DROMON.dir/source/DoFGeom.cpp.o: CMakeFiles/DROMON.dir/flags.make
CMakeFiles/DROMON.dir/source/DoFGeom.cpp.o: ../source/DoFGeom.cpp
CMakeFiles/DROMON.dir/source/DoFGeom.cpp.o: CMakeFiles/DROMON.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/skasdorf/Documents/DROMON-main/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object CMakeFiles/DROMON.dir/source/DoFGeom.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/DROMON.dir/source/DoFGeom.cpp.o -MF CMakeFiles/DROMON.dir/source/DoFGeom.cpp.o.d -o CMakeFiles/DROMON.dir/source/DoFGeom.cpp.o -c /home/skasdorf/Documents/DROMON-main/source/DoFGeom.cpp

CMakeFiles/DROMON.dir/source/DoFGeom.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DROMON.dir/source/DoFGeom.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/skasdorf/Documents/DROMON-main/source/DoFGeom.cpp > CMakeFiles/DROMON.dir/source/DoFGeom.cpp.i

CMakeFiles/DROMON.dir/source/DoFGeom.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DROMON.dir/source/DoFGeom.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/skasdorf/Documents/DROMON-main/source/DoFGeom.cpp -o CMakeFiles/DROMON.dir/source/DoFGeom.cpp.s

CMakeFiles/DROMON.dir/source/DoFBase.cpp.o: CMakeFiles/DROMON.dir/flags.make
CMakeFiles/DROMON.dir/source/DoFBase.cpp.o: ../source/DoFBase.cpp
CMakeFiles/DROMON.dir/source/DoFBase.cpp.o: CMakeFiles/DROMON.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/skasdorf/Documents/DROMON-main/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building CXX object CMakeFiles/DROMON.dir/source/DoFBase.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/DROMON.dir/source/DoFBase.cpp.o -MF CMakeFiles/DROMON.dir/source/DoFBase.cpp.o.d -o CMakeFiles/DROMON.dir/source/DoFBase.cpp.o -c /home/skasdorf/Documents/DROMON-main/source/DoFBase.cpp

CMakeFiles/DROMON.dir/source/DoFBase.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DROMON.dir/source/DoFBase.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/skasdorf/Documents/DROMON-main/source/DoFBase.cpp > CMakeFiles/DROMON.dir/source/DoFBase.cpp.i

CMakeFiles/DROMON.dir/source/DoFBase.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DROMON.dir/source/DoFBase.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/skasdorf/Documents/DROMON-main/source/DoFBase.cpp -o CMakeFiles/DROMON.dir/source/DoFBase.cpp.s

CMakeFiles/DROMON.dir/source/DoFGeomBase.cpp.o: CMakeFiles/DROMON.dir/flags.make
CMakeFiles/DROMON.dir/source/DoFGeomBase.cpp.o: ../source/DoFGeomBase.cpp
CMakeFiles/DROMON.dir/source/DoFGeomBase.cpp.o: CMakeFiles/DROMON.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/skasdorf/Documents/DROMON-main/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Building CXX object CMakeFiles/DROMON.dir/source/DoFGeomBase.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/DROMON.dir/source/DoFGeomBase.cpp.o -MF CMakeFiles/DROMON.dir/source/DoFGeomBase.cpp.o.d -o CMakeFiles/DROMON.dir/source/DoFGeomBase.cpp.o -c /home/skasdorf/Documents/DROMON-main/source/DoFGeomBase.cpp

CMakeFiles/DROMON.dir/source/DoFGeomBase.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DROMON.dir/source/DoFGeomBase.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/skasdorf/Documents/DROMON-main/source/DoFGeomBase.cpp > CMakeFiles/DROMON.dir/source/DoFGeomBase.cpp.i

CMakeFiles/DROMON.dir/source/DoFGeomBase.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DROMON.dir/source/DoFGeomBase.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/skasdorf/Documents/DROMON-main/source/DoFGeomBase.cpp -o CMakeFiles/DROMON.dir/source/DoFGeomBase.cpp.s

# Object files for target DROMON
DROMON_OBJECTS = \
"CMakeFiles/DROMON.dir/source/mesh.cpp.o" \
"CMakeFiles/DROMON.dir/source/GeomBase.cpp.o" \
"CMakeFiles/DROMON.dir/source/MeshGenerator.cpp.o" \
"CMakeFiles/DROMON.dir/source/Point.cpp.o" \
"CMakeFiles/DROMON.dir/source/DataOut.cpp.o" \
"CMakeFiles/DROMON.dir/source/MeshBase.cpp.o" \
"CMakeFiles/DROMON.dir/source/FEBase.cpp.o" \
"CMakeFiles/DROMON.dir/source/FE_HdivMaxOrtho.cpp.o" \
"CMakeFiles/DROMON.dir/source/MultiIndex.cpp.o" \
"CMakeFiles/DROMON.dir/source/FECollection.cpp.o" \
"CMakeFiles/DROMON.dir/source/DoFHandler.cpp.o" \
"CMakeFiles/DROMON.dir/source/DoFGeom.cpp.o" \
"CMakeFiles/DROMON.dir/source/DoFBase.cpp.o" \
"CMakeFiles/DROMON.dir/source/DoFGeomBase.cpp.o"

# External object files for target DROMON
DROMON_EXTERNAL_OBJECTS =

libDROMONd.so: CMakeFiles/DROMON.dir/source/mesh.cpp.o
libDROMONd.so: CMakeFiles/DROMON.dir/source/GeomBase.cpp.o
libDROMONd.so: CMakeFiles/DROMON.dir/source/MeshGenerator.cpp.o
libDROMONd.so: CMakeFiles/DROMON.dir/source/Point.cpp.o
libDROMONd.so: CMakeFiles/DROMON.dir/source/DataOut.cpp.o
libDROMONd.so: CMakeFiles/DROMON.dir/source/MeshBase.cpp.o
libDROMONd.so: CMakeFiles/DROMON.dir/source/FEBase.cpp.o
libDROMONd.so: CMakeFiles/DROMON.dir/source/FE_HdivMaxOrtho.cpp.o
libDROMONd.so: CMakeFiles/DROMON.dir/source/MultiIndex.cpp.o
libDROMONd.so: CMakeFiles/DROMON.dir/source/FECollection.cpp.o
libDROMONd.so: CMakeFiles/DROMON.dir/source/DoFHandler.cpp.o
libDROMONd.so: CMakeFiles/DROMON.dir/source/DoFGeom.cpp.o
libDROMONd.so: CMakeFiles/DROMON.dir/source/DoFBase.cpp.o
libDROMONd.so: CMakeFiles/DROMON.dir/source/DoFGeomBase.cpp.o
libDROMONd.so: CMakeFiles/DROMON.dir/build.make
libDROMONd.so: /usr/lib/x86_64-linux-gnu/libmkl_intel_lp64.so
libDROMONd.so: /usr/lib/x86_64-linux-gnu/libmkl_intel_thread.so
libDROMONd.so: /usr/lib/x86_64-linux-gnu/libmkl_core.so
libDROMONd.so: /usr/lib/x86_64-linux-gnu/libiomp5.so
libDROMONd.so: /usr/lib/gcc/x86_64-linux-gnu/11/libgomp.so
libDROMONd.so: /usr/lib/x86_64-linux-gnu/libpthread.a
libDROMONd.so: CMakeFiles/DROMON.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/skasdorf/Documents/DROMON-main/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_15) "Linking CXX shared library libDROMONd.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/DROMON.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/DROMON.dir/build: libDROMONd.so
.PHONY : CMakeFiles/DROMON.dir/build

CMakeFiles/DROMON.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/DROMON.dir/cmake_clean.cmake
.PHONY : CMakeFiles/DROMON.dir/clean

CMakeFiles/DROMON.dir/depend:
	cd /home/skasdorf/Documents/DROMON-main/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/skasdorf/Documents/DROMON-main /home/skasdorf/Documents/DROMON-main /home/skasdorf/Documents/DROMON-main/build /home/skasdorf/Documents/DROMON-main/build /home/skasdorf/Documents/DROMON-main/build/CMakeFiles/DROMON.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/DROMON.dir/depend

