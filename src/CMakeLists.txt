
find_package(glm REQUIRED) 
set(ARGUS_SOURCES 
        wrinkle_reconstruction/wrinkle_helpers.cpp
        wrinkle_reconstruction/phase_reconstruction.cpp
        get_forces.cpp
        wrinkles_solver.cpp
        thin_shell_stretching.cpp
        thin_shell_bending.cpp
        self_collisions.cpp
        self_collisions_helper.cpp
        vbd.cpp
        wrinkles_vbd.cpp
        sheet_state_io.cpp
        graph_coloring.cpp
        lbfgs_wrapper.cpp
        lbfgsb_wrapper.cpp
        symplectic_euler.cpp
        symplectic_euler_wrinkles.cpp
        symplectic_euler_static_wrinkles.cpp
        material.cpp
        sheet.cpp
        world.cpp
        util.cpp
        gravity.cpp
        rest_shape.cpp
        geometry_functions.cpp
        tfw_stretching.cpp
        dphi_penalty.cpp
        tfw_bending.cpp
        tfw_gravity.cpp
        obstacles.cpp
        obstacle_helpers.cpp
        )

add_library(argus SHARED ${ARGUS_SOURCES} ${ARGUS_HEADERS} )
target_include_directories(argus SYSTEM PUBLIC
        # ${PROJECT_SOURCE_DIR}/include)# ${DISCREGRID_HEADERS})
        ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/vendor/Discregrid/discregrid/include)

target_link_libraries(argus 
        argus_config
        tmd
        lbfgspp
        eigen
        halfedge
        igl::core      
        SuiteSparse::SuiteSparse
        SFFFormulations
        MeshLib
        nlohmann_json::nlohmann_json
        TinyAD
        fcl
        Discregrid
        # OpenVDB::openvdb
        # openvdb
        ${OPTIONAL_DEPENDENCIES}
        ${BOOST_BUNDLE_DIR}/lib/libboost_filesystem.so.1.90.0
        # ${SUITESPARSE_LIBRARIES}
        # CoMISo
        # umfpack
        dlib
        )

if(ARGUS_USE_BOGUS)
    target_link_libraries(
            argus
            bogus
        )
endif()

