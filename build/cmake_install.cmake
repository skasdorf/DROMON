# Install script for directory: /home/skasdorf/Documents/DROMON-main

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
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
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

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libDROMONd.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libDROMONd.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libDROMONd.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/skasdorf/Documents/DROMON-main/build/libDROMONd.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libDROMONd.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libDROMONd.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libDROMONd.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/DROMON" TYPE FILE FILES
    "/home/skasdorf/Documents/DROMON-main/include/AdjointExcitations.h"
    "/home/skasdorf/Documents/DROMON-main/include/DoFMask.h"
    "/home/skasdorf/Documents/DROMON-main/include/MatrixSolving.h"
    "/home/skasdorf/Documents/DROMON-main/include/config.h"
    "/home/skasdorf/Documents/DROMON-main/include/EFIEIntegrator.h"
    "/home/skasdorf/Documents/DROMON-main/include/MatrixSolvingPolicies.h"
    "/home/skasdorf/Documents/DROMON-main/include/DataOut.h"
    "/home/skasdorf/Documents/DROMON-main/include/ErrorEstimation.h"
    "/home/skasdorf/Documents/DROMON-main/include/MeshBase.h"
    "/home/skasdorf/Documents/DROMON-main/include/DIRECTFN_ET_Bounds_Functions.h"
    "/home/skasdorf/Documents/DROMON-main/include/Excitations.h"
    "/home/skasdorf/Documents/DROMON-main/include/MeshGenerator.h"
    "/home/skasdorf/Documents/DROMON-main/include/DIRECTFN_ET_Singular.h"
    "/home/skasdorf/Documents/DROMON-main/include/FEBase.h"
    "/home/skasdorf/Documents/DROMON-main/include/mesh.h"
    "/home/skasdorf/Documents/DROMON-main/include/DIRECTFN_Singular.h"
    "/home/skasdorf/Documents/DROMON-main/include/FECollection.h"
    "/home/skasdorf/Documents/DROMON-main/include/MultiIndex.h"
    "/home/skasdorf/Documents/DROMON-main/include/DIRECTFN_ST_Bounds_Functions.h"
    "/home/skasdorf/Documents/DROMON-main/include/FE_HdivMaxOrtho.h"
    "/home/skasdorf/Documents/DROMON-main/include/Point.h"
    "/home/skasdorf/Documents/DROMON-main/include/DIRECTFN_ST_Singular.h"
    "/home/skasdorf/Documents/DROMON-main/include/GalerkinSystem.h"
    "/home/skasdorf/Documents/DROMON-main/include/PostProcessing.h"
    "/home/skasdorf/Documents/DROMON-main/include/DIRECTFN_VT_Bounds_Functions.h"
    "/home/skasdorf/Documents/DROMON-main/include/GeomBase.h"
    "/home/skasdorf/Documents/DROMON-main/include/QuadratureCollection.h"
    "/home/skasdorf/Documents/DROMON-main/include/DIRECTFN_VT_Singular.h"
    "/home/skasdorf/Documents/DROMON-main/include/IntegratorBase.h"
    "/home/skasdorf/Documents/DROMON-main/include/Refinement.h"
    "/home/skasdorf/Documents/DROMON-main/include/DoFBase.h"
    "/home/skasdorf/Documents/DROMON-main/include/IteratorRanger.h"
    "/home/skasdorf/Documents/DROMON-main/include/SubMatrix.h"
    "/home/skasdorf/Documents/DROMON-main/include/DoFGeomBase.h"
    "/home/skasdorf/Documents/DROMON-main/include/Kernels.h"
    "/home/skasdorf/Documents/DROMON-main/include/UniformIntegrator.h"
    "/home/skasdorf/Documents/DROMON-main/include/DoFGeom.h"
    "/home/skasdorf/Documents/DROMON-main/include/LegendrePFast.h"
    "/home/skasdorf/Documents/DROMON-main/include/utility.h"
    "/home/skasdorf/Documents/DROMON-main/include/DoFHandler.h"
    "/home/skasdorf/Documents/DROMON-main/include/Materials.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/cmake/DROMONTargets.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/cmake/DROMONTargets.cmake"
         "/home/skasdorf/Documents/DROMON-main/build/CMakeFiles/Export/cmake/DROMONTargets.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/cmake/DROMONTargets-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/cmake/DROMONTargets.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake" TYPE FILE FILES "/home/skasdorf/Documents/DROMON-main/build/CMakeFiles/Export/cmake/DROMONTargets.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake" TYPE FILE FILES "/home/skasdorf/Documents/DROMON-main/build/CMakeFiles/Export/cmake/DROMONTargets-debug.cmake")
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/cmake" TYPE FILE FILES
    "/home/skasdorf/Documents/DROMON-main/build/DROMONConfig.cmake"
    "/home/skasdorf/Documents/DROMON-main/build/DROMONConfigVersion.cmake"
    )
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/skasdorf/Documents/DROMON-main/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
