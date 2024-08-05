cmake_minimum_required(VERSION 3.16)
project(dromon_tests)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}
$HOME/Documents/DROMON-main/cmake)
set(CMAKE_CXX_STANDARD 17)

# Set up the various  libraries  we require:
#
#Find MKL BLAS and LAPACK
#
#set(BLA_VENDOR "Intel10_64lp")
find_package(BLAS REQUIRED)
#set(BLA_VENDOR "Intel10_64lp")
find_package(LAPACK REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED)
find_package(GSL REQUIRED)
find_package(DROMON REQUIRED)

# For comparing output files
find_program(NUMDIFF numdiff)

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -DMKL_Complex16=\"std::complex<double>\"")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG -DMKL_Complex16=\"std::complex<double>\"")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-DNDEBUG -DMKL_Complex16=\"std::complex<double>\"")


#HOPS
add_executable(hops test_HOPS.cpp)
target_link_libraries(hops DROMON::DROMON)



	
