cmake_minimum_required(VERSION 3.5...3.28)
include(FetchContent)
project(thinshells)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

FetchContent_Declare(
    libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG v2.4.0
)
FetchContent_MakeAvailable(libigl)
FetchContent_MakeAvailable(eigen)
igl_include(glfw)
igl_include(imgui)

# Add autodiff to the project
add_subdirectory(autodiff)  # Replace with the actual path where autodiff is located


# Add your project files
FILE(GLOB SRCFILES src/*.cpp)
add_executable(${PROJECT_NAME} ${SRCFILES}
        src/utils.h
        src/constants.h)

# Link libigl and autodiff libraries to your project
target_link_libraries(${PROJECT_NAME} 
    PRIVATE igl::core igl::imgui igl::glfw autodiff::autodiff
)