cmake_minimum_required(VERSION 3.30)

# Set CMP0167 to allow config-based Boost
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()

project(DROMON
    DESCRIPTION "DROMON: BEM/MoM Library"
    VERSION 0.0.1
)
set(namespace "DROMON")

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Check if CMAKE_INSTALL_PREFIX was set
if (DEFINED CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	message(FATAL_ERROR "CMAKE_INSTALL_PREFIX must be set! You must set this by appending -DCMAKE_INSTALL_PREFIX=<location of install folder>\n")
endif()

#
#Find MKL BLAS and LAPACK
#
#set(BLA_VENDOR "Intel10_64lp")
find_package(BLAS REQUIRED)
#set(BLA_VENDOR "Intel10_64lp")
find_package(LAPACK REQUIRED)
find_package(OpenMP REQUIRED)
find_package(boost_filesystem CONFIG REQUIRED)
find_package(boost_system CONFIG REQUIRED)
find_package(gsl REQUIRED)

add_library(${PROJECT_NAME} SHARED)

file(GLOB SOURCE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp")

target_sources(${PROJECT_NAME}
    PRIVATE
        "${SOURCE_FILES}"
)

target_include_directories(${PROJECT_NAME}
	PRIVATE
	    ${CMAKE_CURRENT_SOURCE_DIR}/source
	PUBLIC
	    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
	    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
	    )
file(GLOB_RECURSE public_headers "${CMAKE_CURRENT_SOURCE_DIR}/include/DROMON/*.h")
    
target_link_libraries(${PROJECT_NAME}
    LAPACK::LAPACK
    BLAS::BLAS
    OpenMP::OpenMP_CXX
    Boost::filesystem
    Boost::system
    GSL::gsl
)
    
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -DMKL_Complex16=\"std::complex<double>\"")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG -DMKL_Complex16=\"std::complex<double>\"")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-DNDEBUG -DMKL_Complex16=\"std::complex<double>\"")

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(Installing)
