cmake_minimum_required(VERSION 3.10)

# Find required packages for testing
find_package(Threads REQUIRED)

# Test source files
set(TEST_SOURCES
    test_geometric_predicates.cpp
    test_delaunay_triangulation.cpp
    test_voronoi_diagram.cpp
    test_point3d.cpp
    test_main.cpp
)

# Source files from main project (needed for linking)
file(GLOB_RECURSE GEOMETRY_SOURCES "${CMAKE_SOURCE_DIR}/src/geometry/*.cpp")
file(GLOB_RECURSE ALGORITHM_SOURCES "${CMAKE_SOURCE_DIR}/src/algorithms/*.cpp")

# Create test executable
add_executable(tests ${TEST_SOURCES} ${GEOMETRY_SOURCES} ${ALGORITHM_SOURCES})

# Include directories
target_include_directories(tests PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/external
    ${CMAKE_SOURCE_DIR}/external/glad/include
)

# Link libraries
target_link_libraries(tests
    Threads::Threads
    OpenGL::GL
    glfw
    glad
)

# Disable specific warnings for tests
if(NOT MSVC)
    target_compile_options(tests PRIVATE -Wno-unused-variable -Wno-unused-but-set-variable)
endif()

# Platform-specific configurations
if(WIN32)
    target_compile_definitions(tests PRIVATE GLFW_DLL)
endif()

if(APPLE)
    target_link_libraries(tests "-framework Cocoa" "-framework IOKit")
endif()

# Copy test data to build directory if needed
# file(COPY test_data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Enable testing
enable_testing()

# Add test
add_test(NAME unit_tests COMMAND tests)

# Set test properties
set_tests_properties(unit_tests PROPERTIES
    TIMEOUT 300  # 5 minutes timeout
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)