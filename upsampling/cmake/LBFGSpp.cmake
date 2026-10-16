include(FetchContent)

FetchContent_Declare(
    LBFGSpp
    GIT_REPOSITORY "https://github.com/yixuan/LBFGSpp.git"
    GIT_TAG "c524a407fb85b74807f53de5a3ca2ddbcc164e54"
)

FetchContent_GetProperties(LBFGSpp)
if(NOT lbfgspp_POPULATED)
    FetchContent_Populate(LBFGSpp)
endif()

# Ensure Eigen is available
find_package(Eigen3 3.3 REQUIRED)

# Create a simple INTERFACE target manually (header-only)
add_library(LBFGSpp INTERFACE)
target_include_directories(LBFGSpp INTERFACE ${lbfgspp_SOURCE_DIR}/include/)
target_link_libraries(LBFGSpp INTERFACE Eigen3::Eigen)
