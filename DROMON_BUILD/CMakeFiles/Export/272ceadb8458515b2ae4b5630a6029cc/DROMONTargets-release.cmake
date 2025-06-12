#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "DROMON::DROMON" for configuration "Release"
set_property(TARGET DROMON::DROMON APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(DROMON::DROMON PROPERTIES
  IMPORTED_IMPLIB_RELEASE "${_IMPORT_PREFIX}/lib/libDROMON.dll.a"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/libDROMON.dll"
  )

list(APPEND _cmake_import_check_targets DROMON::DROMON )
list(APPEND _cmake_import_check_files_for_DROMON::DROMON "${_IMPORT_PREFIX}/lib/libDROMON.dll.a" "${_IMPORT_PREFIX}/bin/libDROMON.dll" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
