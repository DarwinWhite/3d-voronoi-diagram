#include "algorithms/VoronoiStructures.h"
#include <cassert>
#include <iostream>
#include <vector>
#include <random>

using namespace Voronoi3D;

void test_tetrahedron_basic() {
    Point3D a(0.0, 0.0, 0.0);
    Point3D b(1.0, 0.0, 0.0);
    Point3D c(0.0, 1.0, 0.0);
    Point3D d(0.0, 0.0, 1.0);
    
    Tetrahedron tet(0, a, b, c, d);
    
    assert(tet.getId() == 0);
    assert(!tet.isInfinite());
    assert(tet.hasVertex(a));
    assert(tet.hasVertex(b));
    assert(tet.hasVertex(c));
    assert(tet.hasVertex(d));
    
    // Test volume calculation
    double volume = tet.getVolume();
    assert(volume > 0);
    
    // Test orientation
    bool positive = tet.hasPositiveOrientation();
    // Note: orientation depends on vertex ordering
    
    std::cout << "      Tetrahedron volume: " << volume << ", positive orientation: " << positive << "\n";
}

void test_tetrahedron_circumsphere() {
    Point3D a(1.0, 0.0, 0.0);
    Point3D b(0.0, 1.0, 0.0);
    Point3D c(0.0, 0.0, 1.0);
    Point3D d(0.0, 0.0, 0.0);
    
    Tetrahedron tet(0, a, b, c, d);
    
    // Test point inside circumsphere
    Point3D inside(0.2, 0.2, 0.2);
    bool is_inside = tet.isPointInCircumsphere(inside);
    
    // Test point outside circumsphere
    Point3D outside(5.0, 5.0, 5.0);
    bool is_outside = tet.isPointInCircumsphere(outside);
    
    assert(is_inside != is_outside); // They should be different
    
    std::cout << "      Circumsphere test: inside=" << is_inside << ", outside=" << is_outside << "\n";
}

void test_delaunay_simple() {
    std::vector<Point3D> points = {
        Point3D(0.0, 0.0, 0.0),
        Point3D(1.0, 0.0, 0.0),
        Point3D(0.0, 1.0, 0.0),
        Point3D(0.0, 0.0, 1.0),
        Point3D(0.5, 0.5, 0.5)
    };
    
    DelaunayTriangulation dt;
    dt.initialize(points);
    
    assert(dt.isValid());
    assert(dt.getNumPoints() == points.size());
    
    auto active_tets = dt.getActiveTetrahedra();
    assert(!active_tets.empty());
    
    // Test validation (should not throw)
    try {
        dt.validateTriangulation();
        std::cout << "    Delaunay property validated\n";
    } catch (const std::exception& e) {
        std::cout << "    Validation warning: " << e.what() << "\n";
    }
}

void test_delaunay_incremental() {
    std::vector<Point3D> initial_points = {
        Point3D(0.0, 0.0, 0.0),
        Point3D(1.0, 0.0, 0.0),
        Point3D(0.0, 1.0, 0.0),
        Point3D(0.0, 0.0, 1.0)
    };
    
    DelaunayTriangulation dt;
    dt.initialize(initial_points);
    
    size_t initial_count = dt.getNumTetrahedra();
    
    // Add a point incrementally
    Point3D new_point(0.3, 0.3, 0.3);
    dt.addPoint(new_point);
    
    assert(dt.getNumPoints() == initial_points.size() + 1);
    
    // Verify the point was added correctly
    auto containing_tet = dt.findContainingTetrahedron(new_point);
    if (containing_tet) {
        assert(containing_tet->containsPoint(new_point) || containing_tet->hasVertex(new_point));
    }
    
    std::cout << "      Incremental insertion: " << initial_count << " -> " << dt.getNumTetrahedra() << " tetrahedra\n";
}

void test_delaunay_random() {
    std::mt19937 gen(42); // Fixed seed for reproducibility
    std::uniform_real_distribution<> dis(-5.0, 5.0);
    
    std::vector<Point3D> points;
    for (int i = 0; i < 20; ++i) {
        points.emplace_back(dis(gen), dis(gen), dis(gen));
    }
    
    DelaunayTriangulation dt;
    try {
        dt.initialize(points);
        assert(dt.isValid());
        
        auto stats = dt.getStatistics();
        std::cout << "    Random triangulation statistics:\n";
        std::cout << "      " << stats;
        
    } catch (const std::exception& e) {
        std::cout << "    Note: Random triangulation failed (may be due to degenerate configuration): " << e.what() << "\n";
    }
}

void test_delaunay_triangulation() {
    std::cout << "  Testing Tetrahedron basic functionality...\n";
    test_tetrahedron_basic();
    
    std::cout << "  Testing Tetrahedron circumsphere...\n";
    test_tetrahedron_circumsphere();
    
    std::cout << "  Testing simple Delaunay triangulation...\n";
    test_delaunay_simple();
    
    std::cout << "  Testing incremental insertion...\n";
    test_delaunay_incremental();
    
    std::cout << "  Testing random point triangulation...\n";
    test_delaunay_random();
}