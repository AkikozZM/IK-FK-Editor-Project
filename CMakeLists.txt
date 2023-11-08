cmake_minimum_required(VERSION 3.10)

# set the project name
project(CIS562HW)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})

# Set up libraries

add_library(curve STATIC
    ./src/animation/aSplineVec3.h
    ./src/animation/aSplineVec3.cpp
    ./src/animation/aVector.h
    ./src/animation/aVector.cpp
    ./src/animation/aRotation.h
    ./src/animation/aRotation.cpp
    ./src/animation/aSplineQuat.h
    ./src/animation/aSplineQuat.cpp
)

target_include_directories(curve PUBLIC
    ./3rdparty/eigen-3.2.8
    ./src/animation
)

add_library(FKIK STATIC
    ./src/animation/aActor.h
    ./src/animation/aActor.cpp
    ./src/animation/aBVHController.h
    ./src/animation/aBVHController.cpp
    ./src/animation/aIKController.h
    ./src/animation/aIKController.cpp
    ./src/animation/aJoint.h
    ./src/animation/aJoint.cpp
    ./src/animation/aSkeleton.h
    ./src/animation/aSkeleton.cpp
    ./src/animation/aTarget.h
    ./src/animation/aTarget.cpp
    ./src/animation/aTransform.h
    ./src/animation/aTransform.cpp
)

target_include_directories(FKIK PUBLIC
    ./3rdparty/eigen-3.2.8
    ./src/animation
)

add_library(behavior STATIC
    ./src/animation/aBehaviorController.h
    ./src/animation/aBehaviorController.cpp
    ./src/animation/aBehaviors.h
    ./src/animation/aBehaviors.cpp
)

target_include_directories(behavior PUBLIC
    ./3rdparty/eigen-3.2.8
    ./src/animation
)

target_link_libraries(FKIK PUBLIC behavior)

# Set up executables/viewers
# Find OpenGL
find_package(OpenGL REQUIRED)

# Add GLFW
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
add_subdirectory(./3rdparty/glfw-3.3.2)

# Add GLAD
add_subdirectory(./3rdparty/glad)

# Add ImGUI
add_subdirectory(./3rdparty/imgui-1.77)

# Add tinyobj
add_library(tinyobj STATIC
    ./3rdparty/tinyobj/tiny_obj_loader.h
    ./3rdparty/tinyobj/tiny_obj_loader.cc
)

# Add OpenFBX
add_subdirectory(./3rdparty/OpenFBX)

add_executable(FKViewer
    ./src/viewer/FKMain.cpp
    ./src/viewer/FKViewer.h
    ./src/viewer/FKViewer.cpp
    ./src/viewer/camera.cpp
    ./src/viewer/camera.h
    ./src/viewer/objmodel.h
    ./src/viewer/objmodel.cpp
    ./src/viewer/shader.h
    ./src/viewer/viewer.cpp
    ./src/viewer/viewer.h
    ./src/viewer/FBXModel.h
    ./src/viewer/FBXModel.cpp
    ./src/viewer/drawable.h
    ./src/viewer/drawable.cpp
    ./src/viewer/utils.cpp
    ./src/viewer/utils.h
)

target_include_directories(FKViewer PUBLIC 
    ./src/viewer
    ./src/animation
    ./3rdparty/imgui-1.77
    ./3rdparty/glm
    ./3rdparty/tinyobj
)

target_link_libraries(FKViewer PUBLIC curve FKIK behavior glad glfw imgui tinyobj OpenFBX)

# Set up Unity plugins
add_library(CurvePlugin SHARED
    ./src/plugin/Plugin.h
    ./src/plugin/CurvePlugin.cpp
)

add_library(FKIKPlugin SHARED
    ./src/plugin/Plugin.h
    ./src/plugin/FKIKPlugin.cpp
)

target_link_libraries(CurvePlugin PUBLIC curve)
target_link_libraries(FKIKPlugin PUBLIC curve FKIK behavior)


set(PLUGIN_OUTPUT_DIR "${PROJECT_SOURCE_DIR}/UnityPlugin/Assets/Plugins")
message(${PLUGIN_OUTPUT_DIR})

set_target_properties(CurvePlugin PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY_DEBUG "${PLUGIN_OUTPUT_DIR}"
    LIBRARY_OUTPUT_DIRECTORY_DEBUG "${PLUGIN_OUTPUT_DIR}"
    ARCHIVE_OUTPUT_DIRECTORY_DEBUG "${PLUGIN_OUTPUT_DIR}"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${PLUGIN_OUTPUT_DIR}"
    LIBRARY_OUTPUT_DIRECTORY_RELEASE "${PLUGIN_OUTPUT_DIR}"
    ARCHIVE_OUTPUT_DIRECTORY_RELEASE "${PLUGIN_OUTPUT_DIR}"
    PREFIX "__"
)
set_target_properties(FKIKPlugin PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY_DEBUG "${PLUGIN_OUTPUT_DIR}"
    LIBRARY_OUTPUT_DIRECTORY_DEBUG "${PLUGIN_OUTPUT_DIR}"
    ARCHIVE_OUTPUT_DIRECTORY_DEBUG "${PLUGIN_OUTPUT_DIR}"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${PLUGIN_OUTPUT_DIR}"
    LIBRARY_OUTPUT_DIRECTORY_RELEASE "${PLUGIN_OUTPUT_DIR}"
    ARCHIVE_OUTPUT_DIRECTORY_RELEASE "${PLUGIN_OUTPUT_DIR}"
    PREFIX "__"
)

