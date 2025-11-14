#pragma once

#include "Point3D.h"
#include <array>

namespace Voronoi3D {

/**
 * @brief Robust geometric predicates for 3D Voronoi/Delaunay computations
 * 
 * These functions implement critical geometric tests needed for Delaunay triangulation
 * and Voronoi diagram construction, with careful attention to numerical stability.
 */
namespace GeometricPredicates {

/**
 * @brief Determinant of a 3x3 matrix
 * Used internally by orientation and in-sphere tests
 */
inline double determinant3x3(
    double a11, double a12, double a13,
    double a21, double a22, double a23,
    double a31, double a32, double a33
) {
    return a11 * (a22 * a33 - a23 * a32) -
           a12 * (a21 * a33 - a23 * a31) +
           a13 * (a21 * a32 - a22 * a31);
}

/**
 * @brief Determinant of a 4x4 matrix
 * Used for in-sphere test
 */
inline double determinant4x4(
    double a11, double a12, double a13, double a14,
    double a21, double a22, double a23, double a24,
    double a31, double a32, double a33, double a34,
    double a41, double a42, double a43, double a44
) {
    return a11 * determinant3x3(a22, a23, a24, a32, a33, a34, a42, a43, a44) -
           a12 * determinant3x3(a21, a23, a24, a31, a33, a34, a41, a43, a44) +
           a13 * determinant3x3(a21, a22, a24, a31, a32, a34, a41, a42, a44) -
           a14 * determinant3x3(a21, a22, a23, a31, a32, a33, a41, a42, a43);
}

/**
 * @brief 3D orientation test
 * 
 * Determines the orientation of four points in 3D space.
 * Returns positive if d is above the plane defined by a, b, c (when viewed from above),
 * negative if below, and zero if coplanar.
 * 
 * @param a First point
 * @param b Second point
 * @param c Third point
 * @param d Fourth point to test
 * @return Orientation indicator (positive/negative/zero)
 */
inline double orient3D(const Point3D& a, const Point3D& b, const Point3D& c, const Point3D& d) {
    // Calculate the signed volume of the tetrahedron formed by the four points
    // This is equivalent to the determinant:
    // |ax  ay  az  1|
    // |bx  by  bz  1|
    // |cx  cy  cz  1|
    // |dx  dy  dz  1|
    
    double adx = a.x - d.x;
    double ady = a.y - d.y;
    double adz = a.z - d.z;
    double bdx = b.x - d.x;
    double bdy = b.y - d.y;
    double bdz = b.z - d.z;
    double cdx = c.x - d.x;
    double cdy = c.y - d.y;
    double cdz = c.z - d.z;

    return determinant3x3(
        adx, ady, adz,
        bdx, bdy, bdz,
        cdx, cdy, cdz
    );
}

/**
 * @brief In-sphere test (critical for Delaunay triangulation)
 * 
 * Tests whether point e lies inside the sphere defined by points a, b, c, d.
 * Returns positive if e is inside the sphere, negative if outside, zero if on sphere.
 * 
 * @param a First sphere point
 * @param b Second sphere point
 * @param c Third sphere point
 * @param d Fourth sphere point
 * @param e Point to test
 * @return Sphere test result (positive/negative/zero)
 */
inline double inSphere(const Point3D& a, const Point3D& b, const Point3D& c, 
                      const Point3D& d, const Point3D& e) {
    // Calculate the determinant:
    // |ax  ay  az  ax²+ay²+az²  1|
    // |bx  by  bz  bx²+by²+bz²  1|
    // |cx  cy  cz  cx²+cy²+cz²  1|
    // |dx  dy  dz  dx²+dy²+dz²  1|
    // |ex  ey  ez  ex²+ey²+ez²  1|
    
    double aex = a.x - e.x;
    double aey = a.y - e.y;
    double aez = a.z - e.z;
    double bex = b.x - e.x;
    double bey = b.y - e.y;
    double bez = b.z - e.z;
    double cex = c.x - e.x;
    double cey = c.y - e.y;
    double cez = c.z - e.z;
    double dex = d.x - e.x;
    double dey = d.y - e.y;
    double dez = d.z - e.z;

    double aexbey = aex * bey;
    double bexaey = bex * aey;
    double ab = aexbey - bexaey;
    double bexcey = bex * cey;
    double cexbey = cex * bey;
    double bc = bexcey - cexbey;
    double cexdey = cex * dey;
    double dexcey = dex * cey;
    double cd = cexdey - dexcey;
    double dexaey = dex * aey;
    double aexdey = aex * dey;
    double da = dexaey - aexdey;
    double aexcey = aex * cey;
    double cexaey = cex * aey;
    double ac = aexcey - cexaey;
    double bexdey = bex * dey;
    double dexbey = dex * bey;
    double bd = bexdey - dexbey;

    double abc = aez * bc - bez * ac + cez * ab;
    double bcd = bez * cd - cez * bd + dez * bc;
    double cda = cez * da + dez * ac + aez * cd;
    double dab = dez * ab + aez * bd + bez * da;

    double alift = aex * aex + aey * aey + aez * aez;
    double blift = bex * bex + bey * bey + bez * bez;
    double clift = cex * cex + cey * cey + cez * cez;
    double dlift = dex * dex + dey * dey + dez * dez;

    return dlift * abc - clift * dab + blift * cda - alift * bcd;
}

/**
 * @brief Check if four points are coplanar
 * @param a First point
 * @param b Second point
 * @param c Third point
 * @param d Fourth point
 * @return True if points are coplanar
 */
inline bool areCoplanar(const Point3D& a, const Point3D& b, const Point3D& c, const Point3D& d) {
    return std::abs(orient3D(a, b, c, d)) < EPSILON;
}

/**
 * @brief Check if four points are cocircular (lie on same sphere)
 * @param a First point
 * @param b Second point
 * @param c Third point
 * @param d Fourth point
 * @return True if points are cocircular
 */
inline bool areCocircular(const Point3D& /*a*/, const Point3D& /*b*/, const Point3D& /*c*/, const Point3D& /*d*/) {
    // For this we need a fifth point to test against the sphere of a,b,c,d
    // This is mainly used for degeneracy detection
    return false; // Placeholder - implement if needed for specific cases
}

/**
 * @brief Calculate circumcenter of four points (center of circumsphere)
 * @param a First point
 * @param b Second point
 * @param c Third point
 * @param d Fourth point
 * @return Circumcenter point
 */
Point3D circumcenter(const Point3D& a, const Point3D& b, const Point3D& c, const Point3D& d);

/**
 * @brief Calculate circumradius squared of four points
 * @param a First point
 * @param b Second point
 * @param c Third point
 * @param d Fourth point
 * @return Circumradius squared
 */
double circumradiusSquared(const Point3D& a, const Point3D& b, const Point3D& c, const Point3D& d);

/**
 * @brief Calculate the signed volume of a tetrahedron
 * @param a First vertex
 * @param b Second vertex
 * @param c Third vertex
 * @param d Fourth vertex
 * @return Signed volume (1/6 of orient3D result)
 */
inline double tetrahedronVolume(const Point3D& a, const Point3D& b, const Point3D& c, const Point3D& d) {
    return orient3D(a, b, c, d) / 6.0;
}

/**
 * @brief Check if a tetrahedron has positive orientation
 * @param a First vertex
 * @param b Second vertex
 * @param c Third vertex
 * @param d Fourth vertex
 * @return True if positively oriented
 */
inline bool isPositiveOrientation(const Point3D& a, const Point3D& b, const Point3D& c, const Point3D& d) {
    return orient3D(a, b, c, d) > EPSILON;
}

} // namespace GeometricPredicates

} // namespace Voronoi3D