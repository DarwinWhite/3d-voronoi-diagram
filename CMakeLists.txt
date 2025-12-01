cmake_minimum_required(VERSION 3.10)
project(Voronoi3D VERSION 1.0.0 LANGUAGES C CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler-specific options
if(MSVC)
    add_compile_options(/W4 /WX)
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
else()
    add_compile_options(-Wall -Wextra -Wpedantic -Werror)
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(-g -fsanitize=address -fsanitize=undefined)
        add_link_options(-fsanitize=address -fsanitize=undefined)
    else()
        add_compile_options(-O3 -DNDEBUG)
    endif()
endif()

# Find packages
find_package(OpenGL REQUIRED)
find_package(glfw3 REQUIRED)
find_package(glm REQUIRED)

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/external)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/external/glad/include)

# GLAD (OpenGL loader)
set(GLAD_SOURCES
    external/glad/src/glad.c
)

add_library(glad STATIC ${GLAD_SOURCES})
target_include_directories(glad PUBLIC external/glad/include)
set_target_properties(glad PROPERTIES C_STANDARD 99)

# Disable warnings for GLAD compilation
if(NOT MSVC)
    target_compile_options(glad PRIVATE -Wno-pedantic -Wno-unused-parameter)
endif()

# Source files
file(GLOB_RECURSE GEOMETRY_SOURCES "src/geometry/*.cpp")
file(GLOB_RECURSE ALGORITHM_SOURCES "src/algorithms/*.cpp")
file(GLOB_RECURSE VISUALIZATION_SOURCES "src/visualization/*.cpp")
file(GLOB_RECURSE UTIL_SOURCES "src/utils/*.cpp")

# Create library for shared code
set(LIBRARY_SOURCES
    ${GEOMETRY_SOURCES}
    ${ALGORITHM_SOURCES}
    ${VISUALIZATION_SOURCES}
    ${UTIL_SOURCES}
)

# Header files
file(GLOB_RECURSE HEADERS "include/*.h" "include/*.hpp")

# Create static library for shared functionality
add_library(voronoi3d_core STATIC ${LIBRARY_SOURCES} ${HEADERS})
target_include_directories(voronoi3d_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/external
    ${CMAKE_CURRENT_SOURCE_DIR}/external/glad/include
)

# Create main visualization executable
add_executable(voronoi3d src/main.cpp)
target_link_libraries(voronoi3d voronoi3d_core)

# Create benchmark executable
add_executable(voronoi3d_benchmark src/benchmark_main.cpp)
target_link_libraries(voronoi3d_benchmark voronoi3d_core)

# Link libraries for shared library
target_link_libraries(voronoi3d_core 
    OpenGL::GL
    glfw
    glad
    ${CMAKE_DL_LIBS}
)

# Link libraries for visualization executable
target_link_libraries(voronoi3d 
    voronoi3d_core
)

# Link libraries for benchmark executable (no OpenGL needed)
target_link_libraries(voronoi3d_benchmark
    voronoi3d_core
)

# GLM is header-only, but we need to ensure it's found
if(NOT glm_FOUND)
    message(STATUS "GLM not found via find_package, attempting manual detection")
    find_path(GLM_INCLUDE_DIR glm/glm.hpp
        HINTS
        /usr/include
        /usr/local/include
        /opt/local/include
        ${GLM_ROOT}
    )
    if(GLM_INCLUDE_DIR)
        target_include_directories(voronoi3d PRIVATE ${GLM_INCLUDE_DIR})
        message(STATUS "Found GLM at ${GLM_INCLUDE_DIR}")
    else()
        message(FATAL_ERROR "GLM not found. Please install libglm-dev or set GLM_ROOT")
    endif()
endif()

# Platform-specific configurations
if(WIN32)
    target_compile_definitions(voronoi3d PRIVATE GLFW_DLL)
endif()

if(APPLE)
    target_link_libraries(voronoi3d "-framework Cocoa" "-framework IOKit")
endif()

# Copy shaders to build directory
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/shaders 
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Testing (optional)
option(BUILD_TESTS "Build test suite" ON)
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Documentation (optional)
option(BUILD_DOCS "Build documentation with Doxygen" OFF)
if(BUILD_DOCS)
    find_package(Doxygen)
    if(DOXYGEN_FOUND)
        set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in)
        set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
        configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
        
        add_custom_target(doc_doxygen ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM
        )
    endif()
endif()

# Install rules
install(TARGETS voronoi3d voronoi3d_benchmark DESTINATION bin)
install(DIRECTORY shaders/ DESTINATION share/voronoi3d/shaders)

# CPack configuration
set(CPACK_PACKAGE_NAME "Voronoi3D")
set(CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "3D Voronoi Diagram Visualization")
set(CPACK_PACKAGE_VENDOR "CSCE 620 Project")
set(CPACK_PACKAGE_CONTACT "darwin_white@example.edu")

if(WIN32)
    set(CPACK_GENERATOR "ZIP;NSIS")
elseif(APPLE)
    set(CPACK_GENERATOR "TGZ;DragNDrop")
else()
    set(CPACK_GENERATOR "TGZ;DEB;RPM")
endif()

include(CPack)

# Print configuration summary
message(STATUS "")
message(STATUS "Configuration Summary:")
message(STATUS "  Project: ${PROJECT_NAME} ${PROJECT_VERSION}")
message(STATUS "  Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  OpenGL: ${OPENGL_VERSION_STRING}")
message(STATUS "  GLFW: Found")
if(glm_FOUND)
    message(STATUS "  GLM: Found")
endif()
message(STATUS "  Tests: ${BUILD_TESTS}")
message(STATUS "  Documentation: ${BUILD_DOCS}")
message(STATUS "")