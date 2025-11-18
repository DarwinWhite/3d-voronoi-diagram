#pragma once

#include <cmath>
#include <iostream>
#include <limits>
#include <vector>
#include <string>

namespace Voronoi3D {

/**
 * @brief Epsilon value for floating-point comparisons
 */
constexpr double EPSILON = 1e-10;

/**
 * @brief 3D point class with robust floating-point operations
 */
class Point3D {
public:
    double x, y, z;

    /**
     * @brief Default constructor - creates point at origin
     */
    Point3D() : x(0.0), y(0.0), z(0.0) {}

    /**
     * @brief Constructor with coordinates
     * @param x X coordinate
     * @param y Y coordinate  
     * @param z Z coordinate
     */
    Point3D(double x, double y, double z) : x(x), y(y), z(z) {}

    /**
     * @brief Copy constructor
     */
    Point3D(const Point3D& other) = default;

    /**
     * @brief Assignment operator
     */
    Point3D& operator=(const Point3D& other) = default;

    // Arithmetic operators
    Point3D operator+(const Point3D& other) const {
        return Point3D(x + other.x, y + other.y, z + other.z);
    }

    Point3D operator-(const Point3D& other) const {
        return Point3D(x - other.x, y - other.y, z - other.z);
    }

    Point3D operator*(double scalar) const {
        return Point3D(x * scalar, y * scalar, z * scalar);
    }

    Point3D operator/(double scalar) const {
        if (std::abs(scalar) < EPSILON) {
            throw std::runtime_error("Division by zero in Point3D");
        }
        return Point3D(x / scalar, y / scalar, z / scalar);
    }

    Point3D& operator+=(const Point3D& other) {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }

    Point3D& operator-=(const Point3D& other) {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }

    Point3D& operator*=(double scalar) {
        x *= scalar;
        y *= scalar;
        z *= scalar;
        return *this;
    }

    Point3D& operator/=(double scalar) {
        if (std::abs(scalar) < EPSILON) {
            throw std::runtime_error("Division by zero in Point3D");
        }
        x /= scalar;
        y /= scalar;
        z /= scalar;
        return *this;
    }

    // Comparison operators with epsilon tolerance
    bool operator==(const Point3D& other) const {
        return std::abs(x - other.x) < EPSILON &&
               std::abs(y - other.y) < EPSILON &&
               std::abs(z - other.z) < EPSILON;
    }

    bool operator!=(const Point3D& other) const {
        return !(*this == other);
    }

    bool operator<(const Point3D& other) const {
        if (std::abs(x - other.x) > EPSILON) return x < other.x;
        if (std::abs(y - other.y) > EPSILON) return y < other.y;
        return z < other.z;
    }

    /**
     * @brief Calculate squared distance to another point
     * @param other The other point
     * @return Squared Euclidean distance
     */
    double distanceSquared(const Point3D& other) const {
        double dx = x - other.x;
        double dy = y - other.y;
        double dz = z - other.z;
        return dx * dx + dy * dy + dz * dz;
    }

    /**
     * @brief Calculate distance to another point
     * @param other The other point
     * @return Euclidean distance
     */
    double distance(const Point3D& other) const {
        return std::sqrt(distanceSquared(other));
    }

    /**
     * @brief Calculate squared magnitude of point as vector from origin
     * @return Squared magnitude
     */
    double magnitudeSquared() const {
        return x * x + y * y + z * z;
    }

    /**
     * @brief Calculate magnitude of point as vector from origin
     * @return Magnitude
     */
    double magnitude() const {
        return std::sqrt(magnitudeSquared());
    }

    /**
     * @brief Normalize the point as a vector (in-place)
     * @return Reference to this point
     */
    Point3D& normalize() {
        double mag = magnitude();
        if (mag > EPSILON) {
            x /= mag;
            y /= mag;
            z /= mag;
        }
        return *this;
    }

    /**
     * @brief Get normalized version of this point as vector
     * @return Normalized point
     */
    Point3D normalized() const {
        Point3D result(*this);
        return result.normalize();
    }

    /**
     * @brief Dot product with another point (treating as vectors)
     * @param other The other point/vector
     * @return Dot product result
     */
    double dot(const Point3D& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    /**
     * @brief Cross product with another point (treating as vectors)
     * @param other The other point/vector
     * @return Cross product result
     */
    Point3D cross(const Point3D& other) const {
        return Point3D(
            y * other.z - z * other.y,
            z * other.x - x * other.z,
            x * other.y - y * other.x
        );
    }

    /**
     * @brief Check if point is at origin (within epsilon tolerance)
     * @return True if at origin
     */
    bool isZero() const {
        return std::abs(x) < EPSILON && std::abs(y) < EPSILON && std::abs(z) < EPSILON;
    }

    /**
     * @brief Check if point coordinates are finite
     * @return True if all coordinates are finite
     */
    bool isFinite() const {
        return std::isfinite(x) && std::isfinite(y) && std::isfinite(z);
    }

    /**
     * @brief String representation of point
     * @return String representation
     */
    std::string toString() const {
        return "(" + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ")";
    }

    // Static utility functions
    /**
     * @brief Create point at origin
     * @return Origin point
     */
    static Point3D origin() {
        return Point3D(0.0, 0.0, 0.0);
    }

    /**
     * @brief Calculate centroid of multiple points
     * @param points Vector of points
     * @return Centroid point
     */
    static Point3D centroid(const std::vector<Point3D>& points) {
        if (points.empty()) {
            return Point3D::origin();
        }
        
        Point3D center(0.0, 0.0, 0.0);
        for (const auto& point : points) {
            center += point;
        }
        return center / static_cast<double>(points.size());
    }
};

// Non-member operators
inline Point3D operator*(double scalar, const Point3D& point) {
    return point * scalar;
}

inline std::ostream& operator<<(std::ostream& os, const Point3D& point) {
    os << point.toString();
    return os;
}

// Type aliases for convenience
using Vector3D = Point3D;  // Same mathematical operations apply

} // namespace Voronoi3D