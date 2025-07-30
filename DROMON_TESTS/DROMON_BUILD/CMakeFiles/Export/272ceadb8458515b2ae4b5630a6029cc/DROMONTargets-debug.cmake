#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "DROMON::DROMON" for configuration "Debug"
set_property(TARGET DROMON::DROMON APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(DROMON::DROMON PROPERTIES
  IMPORTED_IMPLIB_DEBUG "${_IMPORT_PREFIX}/lib/libDROMONd.dll.a"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/bin/libDROMONd.dll"
  )

list(APPEND _cmake_import_check_targets DROMON::DROMON )
list(APPEND _cmake_import_check_files_for_DROMON::DROMON "${_IMPORT_PREFIX}/lib/libDROMONd.dll.a" "${_IMPORT_PREFIX}/bin/libDROMONd.dll" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
