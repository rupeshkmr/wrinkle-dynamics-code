set(ARGUS_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/argus.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/argus_collision_forces.hpp
    
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/nanoflann.hpp

    ${CMAKE_CURRENT_SOURCE_DIR}/argus/wrinkle_reconstruction/wrinkle_helpers.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/wrinkle_reconstruction/phase_reconstruction.hpp

    ${CMAKE_CURRENT_SOURCE_DIR}/argus/get_forces.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/obstacle_helpers.hpp

    ${CMAKE_CURRENT_SOURCE_DIR}/argus/scene_animator.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/obstacles.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/geometry_functions.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/symplectic_euler.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/symplectic_euler_wrinkles.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/symplectic_euler_static_wrinkles.hpp

    ${CMAKE_CURRENT_SOURCE_DIR}/argus/material.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/sheet.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/world.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/util.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/constraints.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/common.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/dynamics_object.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/time_integrator.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/sheet_forces.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/gravity.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/rest_shape.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/dphi_penalty.hpp
    # ${CMAKE_CURRENT_SOURCE_DIR}/argus/wrinkle_shell.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/tfw_bending.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/thin_shell_bending.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/thin_shell_stretching.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/graph_coloring.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/optimizers.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/lbfgs_wrapper.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/lbfgsb_wrapper.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/static_solver.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/wrinkles_vbd.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/vbd.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/self_collisions.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/self_collisions_helper.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/argus/wrinkles_solver.hpp
    CACHE INTERNAL "")

