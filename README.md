# 3D Voronoi Diagram Implementation

## Project Overview

This project implements a complete 3D Voronoi diagram system with the following achievements:

- **Core Algorithms**: Bowyer-Watson incremental Delaunay triangulation
- **Dual Extraction**: Voronoi diagram construction from Delaunay dual
- **Robust Predicates**: Numerically stable geometric computations
- **Interactive Visualization**: Real-time 3D rendering with OpenGL
- **Comprehensive Testing**: Unit test suite covering all major components
- **Modern C++**: C++17 implementation with CMake build system

## Key Features

### Algorithmic Components
1. **3D Delaunay Triangulation**:
   - Incremental Bowyer-Watson algorithm
   - Robust circumsphere testing
   - Cavity triangulation and tetrahedron management
   - Handles degenerate cases and numerical precision

2. **Voronoi Diagram Construction**:
   - Dual graph extraction from Delaunay triangulation
   - Voronoi cell construction from circumcenters
   - Site-to-cell mapping and neighbor relationships
   - Spatial queries and geometric extent calculation

3. **Geometric Predicates**:
   - 3D orientation testing with determinant computation
   - In-sphere predicate for Delaunay property validation
   - Circumcenter calculation with numerical stability
   - Adaptive precision for floating-point comparisons

### Visualization Features
1. **Interactive 3D Rendering**:
   - OpenGL 3.3+ based visualization system
   - Real-time point cloud display
   - Camera system with orbit/pan/zoom controls
   - Immediate mode rendering for development

2. **User Interface**:
   - Mouse-based camera manipulation
   - Keyboard controls for point generation
   - Console output of algorithm statistics
   - Real-time computation and display

## Installation & Building

### Prerequisites
```bash
# Required dependencies
sudo apt update
sudo apt install -y build-essential cmake libgl1-mesa-dev libglfw3-dev libglew-dev
```

### Build Process
```bash
# Clone repository
git clone <repository-url>
cd csce620-project

# Create build directory
mkdir build && cd build

# Configure with CMake
cmake ..

# Build main application
make

# Build with tests enabled (optional)
cmake -DBUILD_TESTS=ON ..
make

# Run application
./voronoi3d

# Run test suite (if built with tests)
./tests/tests
```

### Build Outputs
- `voronoi3d` - Main interactive visualization application
- `tests/tests` - Comprehensive unit test suite (if enabled)

### System Requirements
- **OS**: Linux (tested on Ubuntu/Debian)
- **Compiler**: GCC 7+ or Clang 5+ (C++17 support required)
- **Graphics**: OpenGL 3.3+ compatible graphics driver
- **Memory**: Minimum 1GB RAM for typical point sets

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
**Current Implementation:**
```bash
# Run the 3D Voronoi visualization
./voronoi3d

# Run comprehensive test suite
./tests/tests

# Enable tests during build
cmake -DBUILD_TESTS=ON ..
make && make test
```

**Runtime Behavior:**
- Generates random 3D point sets (8 points by default)
- Computes Delaunay triangulation using Bowyer-Watson algorithm
- Extracts Voronoi diagram as dual graph
- Displays points with interactive 3D camera
- Reports algorithm statistics to console
- Supports real-time point regeneration with Space key

### Interactive Controls
- **Mouse**: 
  - Left click + drag: Orbit camera around target
  - Right click + drag: Pan camera view
  - Mouse wheel: Zoom in/out
- **Keyboard**:
  - `Space`: Generate new random point set and recompute Voronoi diagram
  - `R`: Reset camera to default position
  - `ESC`: Exit application

**Application Features:**
- Real-time 3D point visualization
- Automatic Voronoi diagram computation on point generation
- Interactive camera system with orbit/pan/zoom controls
- Console output of algorithm statistics and computation status

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

### Delaunay Triangulation (Bowyer-Watson Algorithm)
**Fully Implemented in VoronoiStructures.cpp:**

1. **Incremental Construction**: 
   - `insertPointBowyerWatson()` - Core insertion algorithm
   - `findConflictingTetrahedra()` - Robust circumsphere testing
   - `retriangulateRemovedRegion()` - Cavity retriangulation

2. **Geometric Robustness**:
   - Proper orientation handling with `hasPositiveOrientation()`
   - Circumsphere intersection using exact predicates
   - Degenerate case handling and error recovery

3. **Data Structure Management**:
   - Efficient tetrahedron storage and neighbor tracking
   - Active tetrahedron set maintenance
   - Incremental point insertion with `addPoint()`

### Voronoi Diagram Extraction
**Dual Graph Construction:**

1. **Cell Construction**:
   - `extractVoronoiDiagram()` - Main dual extraction method
   - `constructCell()` - Per-site cell building from incident tetrahedra
   - Circumcenter collection as Voronoi vertices

2. **Topological Consistency**:
   - Site-to-cell mapping with proper indexing
   - Adjacent tetrahedra discovery for each site
   - Face and edge construction from dual relationships

3. **Query Interface**:
   - `findCell()` - Nearest cell location
   - `getCellsInRegion()` - Spatial range queries
   - `getBoundingBox()` - Geometric extent calculation

### Geometric Predicates
**Numerically Robust Computations:**

1. **Orientation Testing**:
   - `orient3D()` - 3D point orientation with determinant calculation
   - Coplanarity detection and degenerate case handling

2. **In-Sphere Testing**:
   - `inSphere()` - Critical Delaunay property validation
   - Circumcenter calculation with `circumcenter()`
   - Adaptive precision for numerical stability

3. **Determinant Computation**:
   - Optimized 3x3 and 4x4 matrix determinants
   - Epsilon-based floating-point comparisons

## Testing

### Unit Tests
```bash
## Testing

### Comprehensive Test Suite

```bash
# Build with tests enabled
cmake -DBUILD_TESTS=ON ..
make

# Run all tests
./tests/tests

# Expected output:
# === 3D Voronoi Diagram Test Suite ===
# [PASS] Point3D Tests
# [PASS] Geometric Predicates Tests  
# [PASS] Delaunay Triangulation Tests
# [PASS] Voronoi Diagram Tests
```

**Test Coverage:**

1. **Point3D Tests** (`test_point3d.cpp`):
   - Arithmetic operations (+, -, *, /)
   - Comparison operators (==, !=, <)
   - Distance and magnitude calculations
   - Vector operations (dot, cross, normalize)

2. **Geometric Predicates Tests** (`test_geometric_predicates.cpp`):
   - 3x3 matrix determinant computation
   - 3D orientation testing
   - In-sphere predicate validation
   - Circumcenter calculation accuracy
   - Coplanarity detection

3. **Delaunay Triangulation Tests** (`test_delaunay_triangulation.cpp`):
   - Tetrahedron basic functionality
   - Circumsphere intersection testing
   - Simple and incremental triangulation
   - Random point set handling
   - Delaunay property validation

4. **Voronoi Diagram Tests** (`test_voronoi_diagram.cpp`):
   - VoronoiCell construction and properties
   - Diagram computation from point sets
   - Spatial query operations
   - Delaunay-Voronoi duality verification
   - Incremental site addition

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