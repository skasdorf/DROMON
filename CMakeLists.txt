cmake_minimum_required(VERSION 3.16)

set(namespace "DROMON")
project("DROMON"
	DESCRIPTION "DROMON: BEM/MoM Library"
	VERSION 0.0.1
	)
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
find_package(Boost REQUIRED)
find_package(GSL REQUIRED)

add_library(${PROJECT_NAME} SHARED)

set(SOURCE_FILES 
    source/mesh.cpp
    source/GeomBase.cpp
    source/MeshGenerator.cpp
    source/Point.cpp
    source/DataOut.cpp
    source/MeshBase.cpp 
    source/FEBase.cpp 
    source/FE_HdivMaxOrtho.cpp 
    source/MultiIndex.cpp 
    source/FECollection.cpp 
    source/DoFHandler.cpp 
    source/DoFGeom.cpp 
    source/DoFBase.cpp 
    source/DoFGeomBase.cpp 
)

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
set(public_headers
    include/AdjointExcitations.h            
    include/DoFMask.h          
    include/MatrixSolving.h
    include/config.h                        
    include/EFIEIntegrator.h   
    include/MatrixSolvingPolicies.h
    include/DataOut.h                       
    include/ErrorEstimation.h  
    include/MeshBase.h
    include/DIRECTFN_ET_Bounds_Functions.h  
    include/Excitations.h      
    include/MeshGenerator.h
    include/DIRECTFN_ET_Singular.h          
    include/FEBase.h           
    include/mesh.h 
    include/DIRECTFN_Singular.h             
    include/FECollection.h     
    include/MultiIndex.h 
    include/DIRECTFN_ST_Bounds_Functions.h  
    include/FE_HdivMaxOrtho.h  
    include/Point.h 
    include/DIRECTFN_ST_Singular.h          
    include/GalerkinSystem.h   
    include/PostProcessing.h 
    include/DIRECTFN_VT_Bounds_Functions.h  
    include/GeomBase.h         
    include/QuadratureCollection.h 
    include/DIRECTFN_VT_Singular.h          
    include/IntegratorBase.h   
    include/Refinement.h 
    include/DoFBase.h                       
    include/IteratorRanger.h   
    include/SubMatrix.h 
    include/DoFGeomBase.h                   
    include/Kernels.h          
    include/UniformIntegrator.h 
    include/DoFGeom.h                       
    include/LegendrePFast.h    
    include/utility.h 
    include/DoFHandler.h                    
    include/Materials.h
    )
    
    
target_link_libraries(${PROJECT_NAME} LAPACK::LAPACK BLAS::BLAS OpenMP::OpenMP_CXX gsl)
    
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -DMKL_Complex16=\"std::complex<double>\"")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG -DMKL_Complex16=\"std::complex<double>\"")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-DNDEBUG -DMKL_Complex16=\"std::complex<double>\"")

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(Installing)
