################################################################################
# General Information
################################################################################

cmake_minimum_required(VERSION 3.1)
project(WrinkleInterpolationCli)

include(CheckLanguage)
check_language(CUDA)

set(CLI_CUDA_SOURCES "")
set(CLI_CUDA_LIBS "")

if(CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
        find_package(CUDAToolkit REQUIRED)
        list(APPEND CLI_CUDA_SOURCES ../../src/ComplexLoop/ComplexLoopZuenkoCUDA.cu)
        list(APPEND CLI_CUDA_LIBS CUDA::cudart)
        add_compile_definitions(CLI_HAS_CUDA_UPSAMPLER=1)
        message(STATUS "WrinkleInterpolationCli: CUDA upsampling enabled")
else()
        message(WARNING "WrinkleInterpolationCli: CUDA compiler not found; building CPU-only upsampling")
endif()

# Add your project files
# Create an executable
file(GLOB SRCFILES  
    ../../dep/SecStencils/*.h
    ../../include/*.h
    ../../include/*.hpp
    ../../include/MeshLib/*.h
    ../../include/Visualization/*.h
    ../../include/Optimization/*.h
    ../../include/IntrinsicFormula/*.h
    ../../include/ComplexLoop/*.h
    *.h

    ../../dep/SecStencils/*.cpp
    ../../src/*.cpp
    ../../src/MeshLib/*.cpp
    ../../src/Visualization/*.cpp
    ../../src/Optimization/*.cpp
    ../../src/IntrinsicFormula/*.cpp
    ../../src/ComplexLoop/*.cpp
    *.cpp)
    
add_executable(
        ${PROJECT_NAME}_bin
        ${SRCFILES}
        )

set_target_properties(${PROJECT_NAME}_bin PROPERTIES CXX_STANDARD 17)
set_target_properties(${PROJECT_NAME}_bin PROPERTIES CXX_STANDARD_REQUIRED ON)

# Link settings
target_include_directories(${PROJECT_NAME}_bin PUBLIC ${spectra_SOURCE_DIR}/include/Spectra)
target_include_directories(${PROJECT_NAME}_bin PUBLIC ${discregrid_SOURCE_DIR}/discregrid/include)
set(CLI_OPTIONAL_METIS_LIBS "")
if(METIS_FOUND)
        list(APPEND CLI_OPTIONAL_METIS_LIBS ${METIS_LIBRARY})
endif()

target_link_libraries( ${PROJECT_NAME}_bin 
        polyscope 
        geometry-central
        TBB::tbb 
        tbb_static 
        tbbmalloc_static
        igl::core 
        igl_restricted::triangle 
        igl_copyleft::comiso 
        ${SUITESPARSE_LIBRARIES} ${CLI_OPTIONAL_METIS_LIBS}
        CLI11::CLI11
        LBFGSpp
        Discregrid
        ${CLI_CUDA_LIBS})