#include "geometry/GeometricPredicates.h"
#include <cmath>

namespace Voronoi3D {

namespace GeometricPredicates {

Point3D circumcenter(const Point3D& a, const Point3D& b, const Point3D& c, const Point3D& d) {
    // Calculate circumcenter using the formula for the center of a sphere passing through 4 points
    // This involves solving a system of linear equations
    
    // Translate points so that 'a' is at the origin
    Point3D b_rel = b - a;
    Point3D c_rel = c - a;
    Point3D d_rel = d - a;
    
    // Calculate squared magnitudes
    double b_mag2 = b_rel.magnitudeSquared();
    double c_mag2 = c_rel.magnitudeSquared();
    double d_mag2 = d_rel.magnitudeSquared();
    
    // Calculate the determinant for the denominator
    double denom = 2.0 * determinant3x3(
        b_rel.x, b_rel.y, b_rel.z,
        c_rel.x, c_rel.y, c_rel.z,
        d_rel.x, d_rel.y, d_rel.z
    );
    
    if (std::abs(denom) < EPSILON) {
        // Points are coplanar or nearly coplanar
        throw std::runtime_error("Cannot compute circumcenter: points are coplanar");
    }
    
    // Calculate circumcenter coordinates using Cramer's rule
    double cx = determinant3x3(
        b_mag2, b_rel.y, b_rel.z,
        c_mag2, c_rel.y, c_rel.z,
        d_mag2, d_rel.y, d_rel.z
    ) / denom;
    
    double cy = determinant3x3(
        b_rel.x, b_mag2, b_rel.z,
        c_rel.x, c_mag2, c_rel.z,
        d_rel.x, d_mag2, d_rel.z
    ) / denom;
    
    double cz = determinant3x3(
        b_rel.x, b_rel.y, b_mag2,
        c_rel.x, c_rel.y, c_mag2,
        d_rel.x, d_rel.y, d_mag2
    ) / denom;
    
    // Translate back to original coordinate system
    return Point3D(cx, cy, cz) + a;
}

double circumradiusSquared(const Point3D& a, const Point3D& b, const Point3D& c, const Point3D& d) {
    Point3D center = circumcenter(a, b, c, d);
    return center.distanceSquared(a);
}

} // namespace GeometricPredicates

} // namespace Voronoi3D