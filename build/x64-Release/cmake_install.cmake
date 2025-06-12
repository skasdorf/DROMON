# Install script for directory: C:/Users/Crfr/Desktop/Github/DROMON

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "C:/Users/Crfr/Desktop/Github/DROMON/install")
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
  set(CMAKE_OBJDUMP "C:/Strawberry/c/bin/objdump.exe")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY OPTIONAL FILES "C:/Users/Crfr/Desktop/Github/DROMON/build/x64-Release/libDROMON.dll.a")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE SHARED_LIBRARY FILES "C:/Users/Crfr/Desktop/Github/DROMON/build/x64-Release/libDROMON.dll")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/libDROMON.dll" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/libDROMON.dll")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "C:/Strawberry/c/bin/strip.exe" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/libDROMON.dll")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/DROMON" TYPE FILE FILES
    "C:/Users/Crfr/Desktop/Github/DROMON/include/AdjointExcitations.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DoFMask.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/MatrixSolving.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/config.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/EFIEIntegrator.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/MatrixSolvingPolicies.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DataOut.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/ErrorEstimation.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/MeshBase.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DIRECTFN_ET_Bounds_Functions.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/Excitations.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/MeshGenerator.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DIRECTFN_ET_Singular.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/FEBase.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/mesh.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DIRECTFN_Singular.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/FECollection.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/MultiIndex.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DIRECTFN_ST_Bounds_Functions.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/FE_HdivMaxOrtho.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/Point.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DIRECTFN_ST_Singular.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/GalerkinSystem.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/PostProcessing.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DIRECTFN_VT_Bounds_Functions.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/GeomBase.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/QuadratureCollection.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DIRECTFN_VT_Singular.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/IntegratorBase.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/Refinement.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DoFBase.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/IteratorRanger.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/SubMatrix.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DoFGeomBase.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/Kernels.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/UniformIntegrator.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DoFGeom.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/LegendrePFast.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/utility.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/DoFHandler.h"
    "C:/Users/Crfr/Desktop/Github/DROMON/include/Materials.h"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/cmake/DROMONTargets.cmake")
    file(DIFFERENT _cmake_export_file_changed FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/cmake/DROMONTargets.cmake"
         "C:/Users/Crfr/Desktop/Github/DROMON/build/x64-Release/CMakeFiles/Export/272ceadb8458515b2ae4b5630a6029cc/DROMONTargets.cmake")
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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake" TYPE FILE FILES "C:/Users/Crfr/Desktop/Github/DROMON/build/x64-Release/CMakeFiles/Export/272ceadb8458515b2ae4b5630a6029cc/DROMONTargets.cmake")
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake" TYPE FILE FILES "C:/Users/Crfr/Desktop/Github/DROMON/build/x64-Release/CMakeFiles/Export/272ceadb8458515b2ae4b5630a6029cc/DROMONTargets-release.cmake")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake" TYPE FILE FILES
    "C:/Users/Crfr/Desktop/Github/DROMON/build/x64-Release/DROMONConfig.cmake"
    "C:/Users/Crfr/Desktop/Github/DROMON/build/x64-Release/DROMONConfigVersion.cmake"
    )
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "C:/Users/Crfr/Desktop/Github/DROMON/build/x64-Release/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
if(CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_COMPONENT MATCHES "^[a-zA-Z0-9_.+-]+$")
    set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
  else()
    string(MD5 CMAKE_INST_COMP_HASH "${CMAKE_INSTALL_COMPONENT}")
    set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INST_COMP_HASH}.txt")
    unset(CMAKE_INST_COMP_HASH)
  endif()
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "C:/Users/Crfr/Desktop/Github/DROMON/build/x64-Release/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
