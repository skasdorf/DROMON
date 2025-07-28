# Install script for directory: C:/Users/Crfr/Desktop/Github/DROMON

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "C:/Program Files (x86)/dromon_tests")
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

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "C:/msys64/mingw64/bin/objdump.exe")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY OPTIONAL FILES "C:/Users/Crfr/Desktop/Github/DROMON/DROMON_TESTS/DROMON_BUILD/libDROMON.dll.a")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE SHARED_LIBRARY FILES "C:/Users/Crfr/Desktop/Github/DROMON/DROMON_TESTS/DROMON_BUILD/libDROMON.dll")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/libDROMON.dll" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/libDROMON.dll")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "C:/msys64/mingw64/bin/strip.exe" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/libDROMON.dll")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/DROMON" TYPE FILE FILES
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/AdjointExcitations.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/DIRECTFN_ET_Bounds_Functions.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/DIRECTFN_ET_Singular.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/DIRECTFN_ST_Bounds_Functions.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/DIRECTFN_ST_Singular.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/DIRECTFN_Singular.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/DIRECTFN_VT_Bounds_Functions.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/DIRECTFN_VT_Singular.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/DataOut.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/DoFBase.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/DoFGeom.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/DoFGeomBase.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/DoFHandler.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/DoFMask.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/EFIEIntegrator.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/ErrorEstimation.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/Excitations.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/FEBase.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/FECollection.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/FE_HdivMaxOrtho.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/GalerkinSystem.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/GeomBase.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/IntegratorBase.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/IteratorRanger.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/Kernels.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/LegendrePFast.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/Materials.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/MatrixSolving.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/MatrixSolvingPolicies.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/MeshBase.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/MeshGenerator.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/MultiIndex.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/Point.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/PostProcessing.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/QuadratureCollection.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/Refinement.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/SubMatrix.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/UniformIntegrator.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/config.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/mesh.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DROMON/utility.h"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/cmake/DROMONTargets.cmake")
    file(DIFFERENT _cmake_export_file_changed FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/cmake/DROMONTargets.cmake"
         "C:/Users/Crfr/Desktop/Github/DROMON/DROMON_TESTS/DROMON_BUILD/CMakeFiles/Export/272ceadb8458515b2ae4b5630a6029cc/DROMONTargets.cmake")
    if(_cmake_export_file_changed)
      file(GLOB _cmake_old_config_files "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/cmake/DROMONTargets-*.cmake")
      if(_cmake_old_config_files)
        string(REPLACE ";" ", " _cmake_old_config_files_text "${_cmake_old_config_files}")
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/cmake/DROMONTargets.cmake\" will be replaced.  Removing files [${_cmake_old_config_files_text}].")
        unset(_cmake_old_config_files_text)
        file(REMOVE ${_cmake_old_config_files})
      endif()
      unset(_cmake_old_config_files)
    endif()
    unset(_cmake_export_file_changed)
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake" TYPE FILE FILES "C:/Users/Crfr/Desktop/Github/DROMON/DROMON_TESTS/DROMON_BUILD/CMakeFiles/Export/272ceadb8458515b2ae4b5630a6029cc/DROMONTargets.cmake")
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake" TYPE FILE FILES "C:/Users/Crfr/Desktop/Github/DROMON/DROMON_TESTS/DROMON_BUILD/CMakeFiles/Export/272ceadb8458515b2ae4b5630a6029cc/DROMONTargets-release.cmake")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake" TYPE FILE FILES
    "C:/Users/Crfr/Desktop/Github/DROMON/DROMON_TESTS/DROMON_BUILD/DROMONConfig.cmake"
    "C:/Users/Crfr/Desktop/Github/DROMON/DROMON_TESTS/DROMON_BUILD/DROMONConfigVersion.cmake"
    )
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "C:/Users/Crfr/Desktop/Github/DROMON/DROMON_TESTS/DROMON_BUILD/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
