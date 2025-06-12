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

# Test target
add_executable(hops test_HOPS.cpp)

# Link libraries
target_link_libraries(hops
    DROMON
    Boost::filesystem
    Boost::system
    OpenMP::OpenMP_CXX
    LAPACK::LAPACK
    BLAS::BLAS
)

target_include_directories(hops PRIVATE "C:/Program Files (x86)/Intel/oneAPI/mkl/latest/include")
