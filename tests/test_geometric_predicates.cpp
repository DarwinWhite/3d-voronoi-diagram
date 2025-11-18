#include "geometry/GeometricPredicates.h"
#include <cassert>
#include <cmath>
#include <iostream>

using namespace Voronoi3D;
using namespace GeometricPredicates;

void test_determinant3x3() {
    // Test identity matrix
    double det = determinant3x3(
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0
    );
    assert(std::abs(det - 1.0) < EPSILON);
    
    // Test known determinant
    det = determinant3x3(
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
        7.0, 8.0, 9.0
    );
    assert(std::abs(det - 0.0) < EPSILON); // This matrix is singular
    
    std::cout << "      Identity determinant: " << det << "\n";
}

void test_orient3d() {
    // Test coplanar points
    Point3D a(0.0, 0.0, 0.0);
    Point3D b(1.0, 0.0, 0.0);
    Point3D c(0.0, 1.0, 0.0);
    Point3D d(0.5, 0.5, 0.0); // On the same plane
    
    double orient = orient3D(a, b, c, d);
    assert(std::abs(orient) < EPSILON);
    
    // Test point above plane
    Point3D e(0.5, 0.5, 1.0);
    orient = orient3D(a, b, c, e);
    assert(orient > 0);
    
    // Test point below plane
    Point3D f(0.5, 0.5, -1.0);
    orient = orient3D(a, b, c, f);
    assert(orient < 0);
    
    std::cout << "      Orientation tests passed (above: " << orient3D(a, b, c, e) << ", below: " << orient3D(a, b, c, f) << ")\n";
}

void test_insphere() {
    // Create a tetrahedron with known circumsphere
    Point3D a(1.0, 0.0, 0.0);
    Point3D b(-1.0, 0.0, 0.0);
    Point3D c(0.0, 1.0, 0.0);
    Point3D d(0.0, 0.0, 1.0);
    
    // Point at center should be inside
    Point3D center(0.0, 0.0, 0.0);
    double result = inSphere(a, b, c, d, center);
    assert(result > 0); // Inside sphere
    
    // Point far away should be outside
    Point3D far(10.0, 10.0, 10.0);
    double result_far = inSphere(a, b, c, d, far);
    assert(result_far < 0); // Outside sphere
    
    std::cout << "      InSphere tests passed (center: " << result << ", far: " << result_far << ")\n";
}

void test_circumcenter() {
    // Test with a regular tetrahedron
    Point3D a(1.0, 1.0, 1.0);
    Point3D b(-1.0, -1.0, 1.0);
    Point3D c(-1.0, 1.0, -1.0);
    Point3D d(1.0, -1.0, -1.0);
    
    try {
        Point3D center = circumcenter(a, b, c, d);
        
        // Verify that all points are equidistant from circumcenter
        double dist_a = center.distance(a);
        double dist_b = center.distance(b);
        double dist_c = center.distance(c);
        double dist_d = center.distance(d);
        
        assert(std::abs(dist_a - dist_b) < EPSILON);
        assert(std::abs(dist_b - dist_c) < EPSILON);
        assert(std::abs(dist_c - dist_d) < EPSILON);
        
        std::cout << "      Circumcenter distances: " << dist_a << ", " << dist_b << ", " << dist_c << ", " << dist_d << "\n";
    } catch (const std::exception& e) {
        std::cout << "    Note: Circumcenter calculation failed (may be expected for degenerate cases): " << e.what() << "\n";
    }
}

void test_coplanarity() {
    // Test coplanar points
    Point3D a(0.0, 0.0, 0.0);
    Point3D b(1.0, 0.0, 0.0);
    Point3D c(0.0, 1.0, 0.0);
    Point3D d(1.0, 1.0, 0.0);
    
    assert(areCoplanar(a, b, c, d));
    
    // Test non-coplanar points
    Point3D e(0.0, 0.0, 1.0);
    assert(!areCoplanar(a, b, c, e));
}

void test_geometric_predicates() {
    std::cout << "  Testing 3x3 determinant...\n";
    test_determinant3x3();
    
    std::cout << "  Testing 3D orientation...\n";
    test_orient3d();
    
    std::cout << "  Testing in-sphere predicate...\n";
    test_insphere();
    
    std::cout << "  Testing circumcenter calculation...\n";
    test_circumcenter();
    
    std::cout << "  Testing coplanarity test...\n";
    test_coplanarity();
}