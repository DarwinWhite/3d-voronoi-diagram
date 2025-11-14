#include "algorithms/VoronoiStructures.h"
#include "geometry/GeometricPredicates.h"
#include <algorithm>
#include <sstream>
#include <stdexcept>

namespace Voronoi3D {

// =============================================================================
// Tetrahedron Implementation
// =============================================================================

Tetrahedron::Tetrahedron(ID id, const Point3D& a, const Point3D& b, const Point3D& c, const Point3D& d)
    : id_(id), vertices_({a, b, c, d}), is_infinite_(false), circumcenter_computed_(false) {
    
    // Ensure positive orientation
    if (!GeometricPredicates::isPositiveOrientation(a, b, c, d)) {
        // Swap vertices to fix orientation
        std::swap(vertices_[2], vertices_[3]);
    }
}

Tetrahedron::Tetrahedron(ID id, const Point3D& a, const Point3D& b, const Point3D& c)
    : id_(id), vertices_({a, b, c, Point3D()}), is_infinite_(true), circumcenter_computed_(false) {
    // Fourth vertex is undefined for infinite tetrahedron
}

void Tetrahedron::setNeighbor(int face_index, std::shared_ptr<Tetrahedron> neighbor) {
    if (face_index < 0 || face_index >= 4) {
        throw std::out_of_range("Face index out of range");
    }
    neighbors_[face_index] = neighbor;
}

std::shared_ptr<Tetrahedron> Tetrahedron::getNeighbor(int face_index) const {
    if (face_index < 0 || face_index >= 4) {
        throw std::out_of_range("Face index out of range");
    }
    return neighbors_[face_index].lock();
}

int Tetrahedron::findNeighborIndex(std::shared_ptr<Tetrahedron> neighbor) const {
    for (int i = 0; i < 4; ++i) {
        if (neighbors_[i].lock() == neighbor) {
            return i;
        }
    }
    return -1;
}

const Point3D& Tetrahedron::getCircumcenter() {
    if (!circumcenter_computed_) {
        computeCircumcenter();
    }
    return circumcenter_;
}

double Tetrahedron::getCircumradiusSquared() {
    if (!circumcenter_computed_) {
        computeCircumcenter();
    }
    return circumradius_squared_;
}

void Tetrahedron::computeCircumcenter() {
    if (is_infinite_) {
        throw std::runtime_error("Cannot compute circumcenter for infinite tetrahedron");
    }
    
    try {
        circumcenter_ = GeometricPredicates::circumcenter(
            vertices_[0], vertices_[1], vertices_[2], vertices_[3]
        );
        circumradius_squared_ = circumcenter_.distanceSquared(vertices_[0]);
        circumcenter_computed_ = true;
    } catch (const std::exception& e) {
        throw std::runtime_error("Failed to compute circumcenter: " + std::string(e.what()));
    }
}

double Tetrahedron::getVolume() const {
    if (is_infinite_) {
        return std::numeric_limits<double>::infinity();
    }
    return std::abs(GeometricPredicates::tetrahedronVolume(
        vertices_[0], vertices_[1], vertices_[2], vertices_[3]
    ));
}

bool Tetrahedron::hasPositiveOrientation() const {
    if (is_infinite_) {
        return true; // Convention
    }
    return GeometricPredicates::isPositiveOrientation(
        vertices_[0], vertices_[1], vertices_[2], vertices_[3]
    );
}

bool Tetrahedron::containsPoint(const Point3D& point) const {
    if (is_infinite_) {
        return false; // Simplified for now
    }
    
    // Check if point is on the same side of each face as the opposite vertex
    for (int i = 0; i < 4; ++i) {
        std::array<Point3D, 3> face = getFace(i);
        Point3D opposite = vertices_[i];
        
        double orient_opposite = GeometricPredicates::orient3D(face[0], face[1], face[2], opposite);
        double orient_point = GeometricPredicates::orient3D(face[0], face[1], face[2], point);
        
        if (orient_opposite * orient_point < 0) {
            return false;
        }
    }
    return true;
}

bool Tetrahedron::isPointInCircumsphere(const Point3D& point) const {
    if (is_infinite_) {
        return false; // Simplified
    }
    
    return GeometricPredicates::inSphere(
        vertices_[0], vertices_[1], vertices_[2], vertices_[3], point
    ) > 0;
}

bool Tetrahedron::hasVertex(const Point3D& point) const {
    for (const auto& vertex : vertices_) {
        if (vertex == point) {
            return true;
        }
    }
    return false;
}

int Tetrahedron::getVertexIndex(const Point3D& point) const {
    for (int i = 0; i < 4; ++i) {
        if (vertices_[i] == point) {
            return i;
        }
    }
    return -1;
}

std::array<Point3D, 3> Tetrahedron::getFace(int face_index) const {
    if (face_index < 0 || face_index >= 4) {
        throw std::out_of_range("Face index out of range");
    }
    
    // Face opposite to vertex face_index
    std::array<Point3D, 3> face;
    int j = 0;
    for (int i = 0; i < 4; ++i) {
        if (i != face_index) {
            face[j++] = vertices_[i];
        }
    }
    return face;
}

Point3D Tetrahedron::getFaceNormal(int face_index) const {
    std::array<Point3D, 3> face = getFace(face_index);
    Vector3D edge1 = face[1] - face[0];
    Vector3D edge2 = face[2] - face[0];
    return edge1.cross(edge2).normalized();
}

Point3D Tetrahedron::getFaceCenter(int face_index) const {
    std::array<Point3D, 3> face = getFace(face_index);
    return (face[0] + face[1] + face[2]) / 3.0;
}

bool Tetrahedron::isValid() const {
    // Check for degenerate cases
    if (is_infinite_) {
        return true; // Simplified validation for infinite tetrahedra
    }
    
    // Check volume
    double volume = getVolume();
    return volume > EPSILON;
}

std::string Tetrahedron::toString() const {
    std::ostringstream oss;
    oss << "Tetrahedron " << id_ << ": ";
    if (is_infinite_) {
        oss << "(infinite) ";
    }
    oss << vertices_[0] << ", " << vertices_[1] << ", " << vertices_[2];
    if (!is_infinite_) {
        oss << ", " << vertices_[3];
    }
    return oss.str();
}

// =============================================================================
// VoronoiCell Implementation
// =============================================================================

VoronoiCell::VoronoiCell(ID id, const Point3D& site)
    : id_(id), site_(site), is_infinite_(false), volume_(0.0), volume_computed_(false) {
}

void VoronoiCell::addVertex(const Point3D& vertex) {
    vertices_.push_back(vertex);
    volume_computed_ = false;
}

void VoronoiCell::addFace(const std::vector<int>& vertex_indices) {
    faces_.push_back(vertex_indices);
    volume_computed_ = false;
}

void VoronoiCell::addNeighbor(ID cell_id) {
    if (std::find(neighboring_cells_.begin(), neighboring_cells_.end(), cell_id) == neighboring_cells_.end()) {
        neighboring_cells_.push_back(cell_id);
    }
}

double VoronoiCell::getVolume() {
    if (!volume_computed_) {
        computeVolume();
    }
    return volume_;
}

Point3D VoronoiCell::getCentroid() const {
    if (vertices_.empty()) {
        return site_;
    }
    return Point3D::centroid(vertices_);
}

bool VoronoiCell::containsPoint(const Point3D& /*point*/) const {
    // Simple distance test - point belongs to this cell if it's closer to this site
    // than to any other site (this is a simplified implementation)
    return true; // TODO: Implement proper containment test
}

double VoronoiCell::distanceToSite(const Point3D& point) const {
    return site_.distance(point);
}

void VoronoiCell::computeVolume() {
    if (is_infinite_ || faces_.empty()) {
        volume_ = std::numeric_limits<double>::infinity();
        volume_computed_ = true;
        return;
    }
    
    // Compute volume using divergence theorem
    volume_ = 0.0;
    for (const auto& face : faces_) {
        if (face.size() >= 3) {
            // Triangulate face and sum contributions
            Point3D v0 = vertices_[face[0]];
            for (size_t i = 1; i < face.size() - 1; ++i) {
                Point3D v1 = vertices_[face[i]];
                Point3D v2 = vertices_[face[i + 1]];
                
                // Volume contribution from this triangle
                Vector3D to_v0 = v0 - site_;
                Vector3D edge1 = v1 - v0;
                Vector3D edge2 = v2 - v0;
                Vector3D normal = edge1.cross(edge2);
                
                volume_ += to_v0.dot(normal) / 6.0;
            }
        }
    }
    
    volume_ = std::abs(volume_);
    volume_computed_ = true;
}

bool VoronoiCell::isValid() const {
    return !vertices_.empty() && !faces_.empty();
}

void VoronoiCell::clear() {
    vertices_.clear();
    faces_.clear();
    neighboring_cells_.clear();
    volume_computed_ = false;
    volume_ = 0.0;
}

std::string VoronoiCell::toString() const {
    std::ostringstream oss;
    oss << "VoronoiCell " << id_ << " at " << site_ 
        << " with " << vertices_.size() << " vertices and " 
        << faces_.size() << " faces";
    return oss.str();
}

// =============================================================================
// DelaunayTriangulation Implementation (Skeleton)
// =============================================================================

DelaunayTriangulation::DelaunayTriangulation()
    : next_tetrahedron_id_(0), is_initialized_(false) {
}

void DelaunayTriangulation::initialize(const std::vector<Point3D>& points) {
    clear();
    points_ = points;
    
    if (points.size() < 4) {
        throw std::invalid_argument("Need at least 4 points for 3D triangulation");
    }
    
    createBoundingTetrahedron();
    
    // TODO: Implement incremental insertion algorithm in Phase 3
    is_initialized_ = true;
}

void DelaunayTriangulation::addPoint(const Point3D& point) {
    if (!is_initialized_) {
        throw std::runtime_error("Triangulation not initialized");
    }
    
    points_.push_back(point);
    // TODO: Implement Bowyer-Watson insertion in Phase 3
}

void DelaunayTriangulation::clear() {
    points_.clear();
    tetrahedra_.clear();
    active_tetrahedra_.clear();
    next_tetrahedron_id_ = 0;
    is_initialized_ = false;
}

std::vector<DelaunayTriangulation::TetrahedronPtr> DelaunayTriangulation::getActiveTetrahedra() const {
    std::vector<TetrahedronPtr> active;
    for (const auto& tet : tetrahedra_) {
        if (tet && active_tetrahedra_.count(tet->getId()) > 0) {
            active.push_back(tet);
        }
    }
    return active;
}

void DelaunayTriangulation::createBoundingTetrahedron() {
    // Create a large tetrahedron that contains all points
    // This is a placeholder implementation
    double max_coord = 0.0;
    for (const auto& point : points_) {
        max_coord = std::max(max_coord, std::abs(point.x));
        max_coord = std::max(max_coord, std::abs(point.y));
        max_coord = std::max(max_coord, std::abs(point.z));
    }
    
    double scale = max_coord * 10.0;
    bounding_vertices_ = {
        Point3D(-scale, -scale, -scale),
        Point3D(scale, -scale, -scale),
        Point3D(0, scale, -scale),
        Point3D(0, 0, scale)
    };
    
    auto bounding_tet = createTetrahedron(
        bounding_vertices_[0], bounding_vertices_[1], 
        bounding_vertices_[2], bounding_vertices_[3]
    );
    active_tetrahedra_.insert(bounding_tet->getId());
}

DelaunayTriangulation::TetrahedronPtr DelaunayTriangulation::createTetrahedron(
    const Point3D& a, const Point3D& b, const Point3D& c, const Point3D& d) {
    
    auto tet = std::make_shared<Tetrahedron>(getNextTetrahedronId(), a, b, c, d);
    tetrahedra_.push_back(tet);
    return tet;
}

bool DelaunayTriangulation::isValid() const {
    // TODO: Implement validation checks
    return is_initialized_;
}

std::string DelaunayTriangulation::getStatistics() const {
    std::ostringstream oss;
    oss << "DelaunayTriangulation Statistics:\n";
    oss << "  Points: " << points_.size() << "\n";
    oss << "  Tetrahedra: " << active_tetrahedra_.size() << "\n";
    return oss.str();
}

// =============================================================================
// VoronoiDiagram Implementation (Skeleton)
// =============================================================================

VoronoiDiagram::VoronoiDiagram()
    : delaunay_(std::make_unique<DelaunayTriangulation>()), is_computed_(false) {
}

void VoronoiDiagram::compute(const std::vector<Point3D>& sites) {
    clear();
    sites_ = sites;
    
    // Compute Delaunay triangulation
    delaunay_->initialize(sites);
    
    // Extract Voronoi diagram as dual
    extractVoronoiDiagram();
    
    is_computed_ = true;
}

void VoronoiDiagram::addSite(const Point3D& site) {
    sites_.push_back(site);
    if (is_computed_) {
        // TODO: Incremental update
        compute(sites_);
    }
}

void VoronoiDiagram::clear() {
    sites_.clear();
    cells_.clear();
    delaunay_->clear();
    is_computed_ = false;
}

void VoronoiDiagram::extractVoronoiDiagram() {
    // TODO: Implement dual extraction in Phase 3
    // For now, create empty cells for each site
    cells_.clear();
    for (size_t i = 0; i < sites_.size(); ++i) {
        auto cell = std::make_shared<VoronoiCell>(i, sites_[i]);
        cells_.push_back(cell);
    }
}

std::string VoronoiDiagram::getStatistics() const {
    std::ostringstream oss;
    oss << "VoronoiDiagram Statistics:\n";
    oss << "  Sites: " << sites_.size() << "\n";
    oss << "  Cells: " << cells_.size() << "\n";
    oss << delaunay_->getStatistics();
    return oss.str();
}

} // namespace Voronoi3D