add_executable(wrinkle_dynamics wrinkle_dynamics.cpp)

find_package(glm REQUIRED)

if (MSVC)
    set_target_properties(
        wrinkle_dynamics PROPERTIES
        VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/")
endif ()

target_link_libraries(wrinkle_dynamics
        argus
        polyscope
        )
