#include "geometry/Point3D.h"
#include <cassert>
#include <cmath>
#include <iostream>

using namespace Voronoi3D;

void test_point3d_construction() {
    Point3D p1;
    assert(p1.x == 0.0 && p1.y == 0.0 && p1.z == 0.0);
    
    Point3D p2(1.0, 2.0, 3.0);
    assert(p2.x == 1.0 && p2.y == 2.0 && p2.z == 3.0);
    
    Point3D p3(p2);
    assert(p3 == p2);
}

void test_point3d_arithmetic() {
    Point3D p1(1.0, 2.0, 3.0);
    Point3D p2(4.0, 5.0, 6.0);
    
    Point3D sum = p1 + p2;
    assert(sum.x == 5.0 && sum.y == 7.0 && sum.z == 9.0);
    
    Point3D diff = p2 - p1;
    assert(diff.x == 3.0 && diff.y == 3.0 && diff.z == 3.0);
    
    Point3D scaled = p1 * 2.0;
    assert(scaled.x == 2.0 && scaled.y == 4.0 && scaled.z == 6.0);
    
    Point3D divided = p1 / 2.0;
    assert(std::abs(divided.x - 0.5) < EPSILON &&
           std::abs(divided.y - 1.0) < EPSILON &&
           std::abs(divided.z - 1.5) < EPSILON);
}

void test_point3d_comparison() {
    Point3D p1(1.0, 2.0, 3.0);
    Point3D p2(1.0, 2.0, 3.0);
    Point3D p3(1.1, 2.0, 3.0);
    
    assert(p1 == p2);
    assert(p1 != p3);
    assert(p1 < p3);
}

void test_point3d_distance() {
    Point3D p1(0.0, 0.0, 0.0);
    Point3D p2(3.0, 4.0, 0.0);
    
    double dist = p1.distance(p2);
    assert(std::abs(dist - 5.0) < EPSILON);
    
    double dist_sq = p1.distanceSquared(p2);
    assert(std::abs(dist_sq - 25.0) < EPSILON);
}

void test_point3d_magnitude() {
    Point3D p(3.0, 4.0, 0.0);
    
    double mag = p.magnitude();
    assert(std::abs(mag - 5.0) < EPSILON);
    
    double mag_sq = p.magnitudeSquared();
    assert(std::abs(mag_sq - 25.0) < EPSILON);
}

void test_point3d_normalization() {
    Point3D p(3.0, 4.0, 0.0);
    Point3D normalized = p.normalized();
    
    assert(std::abs(normalized.magnitude() - 1.0) < EPSILON);
    assert(std::abs(normalized.x - 0.6) < EPSILON);
    assert(std::abs(normalized.y - 0.8) < EPSILON);
    assert(std::abs(normalized.z - 0.0) < EPSILON);
}

void test_point3d_dot_cross() {
    Point3D p1(1.0, 0.0, 0.0);
    Point3D p2(0.0, 1.0, 0.0);
    
    double dot = p1.dot(p2);
    assert(std::abs(dot - 0.0) < EPSILON);
    
    Point3D cross = p1.cross(p2);
    assert(std::abs(cross.x - 0.0) < EPSILON);
    assert(std::abs(cross.y - 0.0) < EPSILON);
    assert(std::abs(cross.z - 1.0) < EPSILON);
}

void test_point3d() {
    std::cout << "  Testing Point3D construction...\n";
    test_point3d_construction();
    
    std::cout << "  Testing Point3D arithmetic...\n";
    test_point3d_arithmetic();
    
    std::cout << "  Testing Point3D comparison...\n";
    test_point3d_comparison();
    
    std::cout << "  Testing Point3D distance...\n";
    test_point3d_distance();
    
    std::cout << "  Testing Point3D magnitude...\n";
    test_point3d_magnitude();
    
    std::cout << "  Testing Point3D normalization...\n";
    test_point3d_normalization();
    
    std::cout << "  Testing Point3D dot/cross products...\n";
    test_point3d_dot_cross();
}