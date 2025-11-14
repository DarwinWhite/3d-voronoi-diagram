# 3D Voronoi Diagram Construction and Visualization

A comprehensive implementation of 3D Voronoi diagram algorithms with interactive OpenGL visualization.

## Project Overview

This project implements a robust 3D Voronoi diagram construction algorithm from scratch, focusing on both theoretical correctness and practical robustness challenges in three dimensions. The implementation includes:

- **Delaunay Triangulation**: Incremental Bowyer-Watson algorithm for 3D tetrahedralization
- **Voronoi Diagram Extraction**: Dual graph construction from Delaunay triangulation
- **Interactive Visualization**: Real-time 3D rendering with OpenGL and GLFW
- **Robust Geometry**: Careful handling of numerical precision and edge cases

## Features

### Core Algorithm
- 3D Delaunay triangulation using Bowyer-Watson incremental insertion
- Voronoi diagram extraction as the dual of Delaunay triangulation
- Robust geometric predicates with numerical stability
- Support for various point distributions (random, grid, spherical)

### Visualization
- Interactive 3D rendering with OpenGL
- Multiple viewing modes (points, wireframes, solid cells)
- Real-time construction animation
- Camera controls (orbit, pan, zoom)
- Dynamic point insertion and removal

### Testing & Evaluation
- Comprehensive unit tests for geometric operations
- Performance benchmarking tools
- Correctness validation on known configurations
- Support for stress testing with large datasets

## Project Structure

```
csce620-project/
├── src/                    # Source files
│   ├── geometry/          # Geometric classes and algorithms
│   ├── visualization/     # OpenGL rendering and UI
│   ├── algorithms/        # Core Voronoi/Delaunay algorithms
│   └── utils/             # Utility functions and helpers
├── include/               # Header files
├── shaders/              # GLSL shader files
├── tests/                # Unit tests and benchmarks
├── docs/                 # Documentation and reports
├── external/             # Third-party dependencies
├── CMakeLists.txt        # CMake build configuration
└── README.md            # This file
```

## Dependencies

### Required Libraries
- **OpenGL 3.3+**: Core graphics API
- **GLFW 3.3+**: Window management and input handling
- **GLM 0.9.9+**: Mathematics library for graphics
- **GLAD**: OpenGL function loader (included)

### Build Requirements
- **CMake 3.10+**: Build system
- **GCC 7+ or Clang 6+**: C++17 compatible compiler
- **Linux/Windows/macOS**: Cross-platform support

## Building the Project

### Prerequisites (Ubuntu/Debian)
```bash
sudo apt update
sudo apt install build-essential cmake
sudo apt install libglfw3-dev libgl1-mesa-dev libglu1-mesa-dev
sudo apt install libglm-dev
```

### Prerequisites (Other Systems)
- **Windows**: Install Visual Studio with CMake support
- **macOS**: Install Xcode and CMake via Homebrew

### Compilation
```bash
# Clone and navigate to project directory
cd csce620-project

# Create build directory
mkdir build && cd build

# Configure with CMake
cmake ..

# Build the project
make -j$(nproc)

# Run the application
./voronoi3d
```

### CMake Options
```bash
# Debug build
cmake -DCMAKE_BUILD_TYPE=Debug ..

# Release build (optimized)
cmake -DCMAKE_BUILD_TYPE=Release ..

# Enable testing
cmake -DBUILD_TESTS=ON ..
```

## Usage

### Basic Operation
```bash
# Run with default settings
./voronoi3d

# Generate random points
./voronoi3d --points 100 --distribution random

# Load points from file
./voronoi3d --input points.txt

# Enable animation mode
./voronoi3d --animate --speed 1.0
```

### Interactive Controls
- **Mouse**: 
  - Left click + drag: Rotate camera
  - Right click + drag: Pan view
  - Scroll wheel: Zoom in/out
- **Keyboard**:
  - `Space`: Toggle animation
  - `R`: Reset camera view
  - `1-3`: Switch rendering modes
  - `P`: Add random point
  - `C`: Clear all points
  - `ESC`: Exit application

### Configuration Files
Create `config.json` to customize default settings:
```json
{
  "window": {
    "width": 1200,
    "height": 800,
    "title": "3D Voronoi Visualization"
  },
  "rendering": {
    "background_color": [0.1, 0.1, 0.1],
    "point_size": 5.0,
    "line_width": 1.0
  },
  "algorithm": {
    "epsilon": 1e-10,
    "max_points": 10000
  }
}
```

## Algorithm Implementation

### Delaunay Triangulation
The implementation uses the Bowyer-Watson algorithm:
1. Start with a large bounding tetrahedron
2. Insert points incrementally
3. For each point, find containing tetrahedron
4. Remove conflicting tetrahedra
5. Re-triangulate the cavity

### Voronoi Diagram Construction
Voronoi cells are extracted as the dual of the Delaunay triangulation:
1. Each Delaunay vertex corresponds to a Voronoi cell
2. Each Delaunay tetrahedron corresponds to a Voronoi vertex
3. Voronoi edges connect circumcenters of adjacent tetrahedra
4. Handle infinite cells at the convex hull boundary

### Geometric Predicates
Robust implementation of critical geometric tests:
- **Orientation Test**: Determines relative position of points
- **In-Sphere Test**: Checks if a point lies inside a tetrahedron's circumsphere
- **Exact Arithmetic**: Uses adaptive precision for numerical stability

## Testing

### Unit Tests
```bash
# Build and run all tests
cd build
make test

# Run specific test categories
./tests/geometry_tests
./tests/algorithm_tests
./tests/visualization_tests
```

### Performance Benchmarks
```bash
# Run performance tests
./tests/benchmark_suite

# Generate performance report
./tests/benchmark_suite --report performance_report.json
```

### Correctness Validation
```bash
# Test against known configurations
./tests/correctness_tests --validate

# Visual debugging mode
./voronoi3d --debug --step-by-step
```

## Documentation

- **[Algorithm Details](docs/algorithm.md)**: In-depth explanation of the implementation
- **[API Reference](docs/api.md)**: Class and function documentation
- **[Performance Analysis](docs/performance.md)**: Benchmarking results and optimization notes
- **[Implementation Challenges](docs/challenges.md)**: Numerical issues and solutions

## Contributing

### Development Guidelines
1. Follow the existing code style and naming conventions
2. Write comprehensive unit tests for new features
3. Document all public APIs and algorithms
4. Ensure cross-platform compatibility
5. Run performance tests before submitting changes

### Code Style
- Use C++17 features where appropriate
- Follow RAII principles for resource management
- Prefer const-correctness and immutability
- Use meaningful variable and function names
- Keep functions focused and modular

## References

- Sugihara, K., & Iri, M. (1992). Computing the 3D Voronoi Diagram Robustly: An Easy Explanation. *Computer-Aided Design*, 24(8), 437–442.
- Barber, C. B., Dobkin, D. P., & Huhdanpää, H. (1996). The Quickhull Algorithm for Convex Hulls. *ACM Transactions on Mathematical Software*, 22(4), 469–483.
- Aurenhammer, F. (1991). Voronoi Diagrams: A Survey of a Fundamental Geometric Data Structure. *ACM Computing Surveys*, 23(3), 345–405.

## License

This project is developed for educational purposes as part of CSCE 620 - Computational Geometry.

## Author

Darwin White - CSCE 620 Project Implementation

---

*Last updated: November 2025*