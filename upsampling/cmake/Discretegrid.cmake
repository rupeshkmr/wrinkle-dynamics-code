include(FetchContent)

# Discregrid options
set(BUILD_CMD_EXECUTABLE OFF CACHE BOOL "Don't build Discregrid command line tools")

FetchContent_Declare(
    Discregrid
    GIT_REPOSITORY "https://github.com/InteractiveComputerGraphics/Discregrid.git"
    GIT_TAG "ddf20dc0480874bf02e0bdc6ded76c1f101b17fb"
)

# Make Discregrid available - this will build the library
FetchContent_MakeAvailable(Discregrid)
