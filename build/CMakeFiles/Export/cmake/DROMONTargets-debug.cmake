#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "DROMON::DROMON" for configuration "Debug"
set_property(TARGET DROMON::DROMON APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(DROMON::DROMON PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libDROMONd.so"
  IMPORTED_SONAME_DEBUG "libDROMONd.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS DROMON::DROMON )
list(APPEND _IMPORT_CHECK_FILES_FOR_DROMON::DROMON "${_IMPORT_PREFIX}/lib/libDROMONd.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
