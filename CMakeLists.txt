# =============================================================================
# 1. Project Declaration
# =============================================================================
cmake_minimum_required(VERSION 3.16) # Increased for modern FetchContent and policies
project(
    openargus
    VERSION 1.0
    DESCRIPTION "Cloth simulation framework"
    LANGUAGES CXX
)

# =============================================================================
# 2. Standard and Output Configuration
# =============================================================================
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) # Enforce standard compliance
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Use modern variables to control where libraries and executables are built
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# Enable nice folder structures in IDEs like Visual Studio
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# =============================================================================
# 3. Project Build Options
# =============================================================================
message(STATUS "Configuring OpenArgus...")
option(ARGUS_USE_FLOATS "Use floats instead of doubles" OFF)
option(ARGUS_USE_BOGUS "Use so-bogus for sparse matrix operations" OFF)
option(ARGUS_ENABLE_MULTITHREADING "Enable multithreading with OpenMP" OFF)
option(ARGUS_ENABLE_VISUALIZATION "Enable Polyscope visualization" OFF)
option(ARGUS_BUILD_DOCUMENTATION "Build Doxygen API documentation" OFF)
option(ARGUS_CHECKPOINT "Save checkpoints for visualization" OFF)
option(ARGUS_PRINT_DEBINFO "Print debug information" OFF)
option(ARGUS_QUIET "Silence messages inside simulation loop" OFF)
option(ARGUS_ENABLE_SELF_COLLISION "Compile FCL for self-collisions" ON)
option(ARGUS_BUILD_PYBINDINGS "Build python bindings" OFF)

# =============================================================================
# 4. Core Configuration Target
# =============================================================================
# An INTERFACE library is used to manage all common project settings.
# Targets in `src/` and `examples/` will link to this to inherit everything
# below, like compile definitions, linked libraries, and compiler options.
add_library(argus_config INTERFACE)

# --- Add Compile Definitions based on Options ---
target_compile_definitions(argus_config INTERFACE "PROJECT_ROOT_DIR=\"${CMAKE_SOURCE_DIR}/\"")
if(ARGUS_ENABLE_VISUALIZATION)
    target_compile_definitions(argus_config INTERFACE VIZ)
endif()
if(ARGUS_USE_BOGUS)
    target_compile_definitions(argus_config INTERFACE BOGUS)
endif()
if(ARGUS_QUIET)
    target_compile_definitions(argus_config INTERFACE ARGUS_QUIET)
endif()
if(ARGUS_PRINT_DEBINFO)
    target_compile_definitions(argus_config INTERFACE ARGUS_DEBUG)
endif()
if(ARGUS_CHECKPOINT)
    target_compile_definitions(argus_config INTERFACE ARGUS_CHECKPOINT)
endif()
if(ARGUS_USE_FLOATS)
    target_compile_definitions(argus_config INTERFACE ARGUS_USE_FLOATS)
endif()
if(ARGUS_ENABLE_SELF_COLLISION)
    target_compile_definitions(argus_config INTERFACE ARGUS_FCL)
endif()

# --- Multithreading Configuration (OpenMP) ---
if(ARGUS_ENABLE_MULTITHREADING)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        message(STATUS "OpenMP found, enabling multithreading.")
        target_link_libraries(argus_config INTERFACE OpenMP::OpenMP_CXX)
        target_compile_definitions(argus_config INTERFACE USE_OMP)
        if(ARGUS_USE_BOGUS)
            set(OPENMP ON CACHE BOOL "Enable OpenMP in so-bogus" FORCE)
        endif()
    else()
        message(WARNING "OpenMP not found, multithreading disabled.")
        target_compile_definitions(argus_config INTERFACE EIGEN_DONT_PARALLELIZE)
    endif()
else()
    message(STATUS "Multithreading disabled by option.")
    target_compile_definitions(argus_config INTERFACE EIGEN_DONT_PARALLELIZE)
    if(ARGUS_USE_BOGUS)
        set(OPENMP OFF CACHE BOOL "Disable OpenMP in so-bogus" FORCE)
    endif()
endif()

target_compile_definitions(argus_config INTERFACE EIGEN_DONT_PARALLELIZE)
# --- Platform and Build Type Specific Options ---
# For MSVC compiler
if(MSVC)
    target_compile_options(argus_config INTERFACE /MP /bigobj)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)

    # Disable PDB generation on Windows and instead stuff the debug info into the object file
    # This is required for sccache speedups.
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        string(REPLACE "/Zi" "/Z7" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
        string(REPLACE "/Zi" "/Z7" CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}")
    elseif (CMAKE_BUILD_TYPE STREQUAL "Release")
        string(REPLACE "/Zi" "/Z7" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
        string(REPLACE "/Zi" "/Z7" CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")
    elseif (CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
        string(REPLACE "/Zi" "/Z7" CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
        string(REPLACE "/Zi" "/Z7" CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}")
    endif ()
endif()

# For Debug builds (Address Sanitizer)
# A generator expression applies these options only for the Debug configuration
target_compile_options(argus_config INTERFACE $<$<CONFIG:Debug>:-fsanitize=address>)
target_link_options(argus_config INTERFACE $<$<CONFIG:Debug>:-fsanitize=address>)

# =============================================================================
# 5. Dependencies
# =============================================================================
include(FetchContent)
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_LIST_DIR}/cmake")

# --- Eigen system ---
find_package (Eigen3 3.4 QUIET)
if(EIGEN3_FOUND)
    message(STATUS "--    success. find_package() succeeded with Eigen at ${EIGEN3_INCLUDE_DIR}")
# The find_package() failed, manually download the package
else()
    message(STATUS "--    Could not find Eigen anywhere. Downloading it (or using cached download).")

    if(${CMAKE_VERSION} VERSION_LESS 3.11)
        include(FetchContentLocal)
    else()
        include(FetchContent)
    endif()
    FetchContent_Declare(
        eigen
        GIT_REPOSITORY  https://gitlab.com/libeigen/eigen.git
        GIT_TAG         3.4.0
        SOURCE_DIR      "${CMAKE_CURRENT_BINARY_DIR}/eigen-src"
        BINARY_DIR      "${CMAKE_CURRENT_BINARY_DIR}/eigen-build"
    )
    FetchContent_GetProperties(eigen)
    if(NOT eigen_POPULATED)
        FetchContent_Populate(eigen)
        include(EigenChecker)
        eigen3checker(${eigen_SOURCE_DIR} 3.4)
    endif()
    message(STATUS "--    success. Downloaded Eigen in ${eigen_SOURCE_DIR}")
endif()

add_library (eigen INTERFACE)
add_library (Eigen3::Eigen ALIAS eigen)
# added SYSTEM here to supress warnings
target_include_directories(eigen INTERFACE
    $<BUILD_INTERFACE:${EIGEN3_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:${INCLUDE_INSTALL_DIR}>
)


# --- FCL (system) ---
# find_package(FCL REQUIRED)

# --- SuiteSparse (system) ---
message(STATUS "Finding SuiteSparse...")
find_package(SuiteSparse REQUIRED)
add_library(SuiteSparse::SuiteSparse INTERFACE IMPORTED)
target_include_directories(SuiteSparse::SuiteSparse INTERFACE ${SUITESPARSE_INCLUDE_DIRS})
target_link_libraries(SuiteSparse::SuiteSparse INTERFACE ${SUITESPARSE_LIBRARIES})

# --- nlohmann/json (System or Cached) ---
# Try to find the system-installed nlohmann/json package
find_package(nlohmann_json CONFIG)
find_package(fcl REQUIRED)
# If the package is not found, fall back to FetchContent
if (NOT nlohmann_json_FOUND)
    message(STATUS "nlohmann/json not found, falling back to FetchContent.")
    FetchContent_Declare(
        json
        URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz
        # Optionally, set this to true to use the cached version if available
        # if the external content has been downloaded once.
        # This still requires internet access for the first download.
        # Use SYSTEM or NO_CACHE to prevent caching in the build directory.
        # Here, we keep the default (caching is usually on) for offline use after 1st build.
    )
    FetchContent_MakeAvailable(json)
    
    # Define a target alias for consistency, if FetchContent was used
    # The actual target name from FetchContent might be `nlohmann_json::nlohmann_json` or similar, 
    # but the package approach is often cleaner. Assuming `nlohmann_json` is the target name.
    # Note: For simple single-header libs like nlohmann/json, often it's just an include path.
else()
    message(STATUS "Found system nlohmann/json.")
endif()
# --- libigl (System or Cached) ---
# Try to find the system-installed libigl package
find_package(libigl CONFIG)

# If the package is not found, fall back to FetchContent
if (NOT libigl_FOUND)
    message(STATUS "libigl not found, falling back to FetchContent.")
    
    FetchContent_Declare(
        libigl
        GIT_REPOSITORY https://github.com/libigl/libigl.git
        GIT_TAG main
        # Again, this will use the cached version (in the build tree) on subsequent runs 
        # *after* the initial successful fetch, allowing for offline compilation.
    )

    # Propagate options to the fetched content before making it available
    option(LIBIGL_WITH_GLFW "Build with GLFW support" OFF)
    option(LIBIGL_WITH_IMGUI "Build with ImGui support" OFF)
    option(LIBIGL_WITH_PNG "Build with PNG support" OFF)
    option(LIBIGL_WITH_OPENGL "Build with OpenGL support" OFF)
    
    FetchContent_MakeAvailable(libigl)
else()
    message(STATUS "Found system libigl.")
    
    # If using the system version, you don't need to set the options above 
    # as they should have been set during the system build/installation.
endif()
# =============================================================================
# 6. Documentation (Doxygen)
# =============================================================================
if(ARGUS_BUILD_DOCUMENTATION)
    find_package(Doxygen REQUIRED)
    message(STATUS "Doxygen found, 'docs' target will be available.")

    set(DOXYGEN_IN "${PROJECT_SOURCE_DIR}/docs/Doxyfile.in")
    set(DOXYGEN_OUT "${PROJECT_BINARY_DIR}/Doxyfile")
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

    add_custom_target(docs ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen..."
        VERBATIM
    )
endif()

# =============================================================================
# 7. Add Source Code Subdirectories
# =============================================================================
# NOTE: Your targets inside these directories should link to argus_config, e.g.:
#
#   add_library(my_lib ...)
#   target_link_libraries(my_lib PUBLIC
#       argus_config
#       SuiteSparse::SuiteSparse
#       igl::core
#   )
#

# =============================================================================
# Bundled Boost (self-contained; no system Boost installation required)
# =============================================================================
set(BOOST_BUNDLE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/vendor/boost")
target_include_directories(argus_config INTERFACE "${BOOST_BUNDLE_DIR}/include")
set(CMAKE_BUILD_RPATH "${BOOST_BUNDLE_DIR}/lib")
set(CMAKE_INSTALL_RPATH "$ORIGIN/../../vendor/boost/lib")

add_subdirectory(vendor)
add_subdirectory(src)
add_subdirectory(examples)
