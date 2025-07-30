cmake_minimum_required(VERSION 3.30)
# Set CMP0167 to allow config-based Boost
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()

project(dromon_tests)
set(CMAKE_CXX_STANDARD 17)

# Use CMake modules from your local project
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR} "C:/Users/Crfr/Desktop/Github/DROMON/cmake")

# Find dependencies
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(OpenMP REQUIRED)
find_package(boost_filesystem CONFIG REQUIRED)
find_package(boost_system CONFIG REQUIRED)
find_package(gsl REQUIRED)
add_subdirectory("C:/Users/Crfr/Desktop/Github/DROMON" DROMON_BUILD)


# Compiler flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -DMKL_Complex16=\"std::complex<double>\"")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG -DMKL_Complex16=\"std::complex<double>\"")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-DNDEBUG -DMKL_Complex16=\"std::complex<double>\"")

# Collect all test source files
file(GLOB TEST_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

foreach(TEST_FILE ${TEST_SOURCES})
    # Extract the filename without extension
    get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
    
    # Create executable for each test
    add_executable(${TEST_NAME} ${TEST_FILE})
    
    # Link libraries
    target_link_libraries(${TEST_NAME}
        DROMON
        Boost::filesystem
        Boost::system
        OpenMP::OpenMP_CXX
        LAPACK::LAPACK
        BLAS::BLAS
    )

    # Include directories
    target_include_directories(${TEST_NAME} PRIVATE
        "C:/Program Files (x86)/Intel/oneAPI/mkl/latest/include"
    )

    # Optionally register with CTest
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach()
