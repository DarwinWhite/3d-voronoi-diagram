#include "algorithms/VoronoiStructures.h"
#include "geometry/GeometricPredicates.h"
#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <set>
#include <array>
#include <limits>

namespace Voronoi3D {

// =============================================================================
// Tetrahedron Implementation
// =============================================================================

Tetrahedron::Tetrahedron(ID id, const Point3D& a, const Point3D& b, const Point3D& c, const Point3D& d)
    : id_(id), vertices_({a, b, c, d}), is_infinite_(false), circumcenter_computed_(false) {
    // Initialize all neighbors to null
    for (auto& neighbor : neighbors_) {
        neighbor.reset();
    }
}

Tetrahedron::Tetrahedron(ID id, const Point3D& a, const Point3D& b, const Point3D& c)
    : id_(id), vertices_({a, b, c, Point3D()}), is_infinite_(true), circumcenter_computed_(false) {
    for (auto& neighbor : neighbors_) {
        neighbor.reset();
    }
}

void Tetrahedron::setNeighbor(int face_index, std::shared_ptr<Tetrahedron> neighbor) {
    if (face_index >= 0 && face_index < 4) {
        neighbors_[face_index] = neighbor;
    }
}

std::shared_ptr<Tetrahedron> Tetrahedron::getNeighbor(int face_index) const {
    if (face_index >= 0 && face_index < 4) {
        return neighbors_[face_index].lock();
    }
    return nullptr;
}

int Tetrahedron::findNeighborIndex(std::shared_ptr<Tetrahedron> neighbor) const {
    for (int i = 0; i < 4; ++i) {
        auto n = neighbors_[i].lock();
        if (n && n->getId() == neighbor->getId()) {
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
        throw std::runtime_error("Cannot compute circumcenter of infinite tetrahedron");
    }
    
    try {
        circumcenter_ = GeometricPredicates::circumcenter(
            vertices_[0], vertices_[1], vertices_[2], vertices_[3]
        );
        circumradius_squared_ = circumcenter_.distanceSquared(vertices_[0]);
        circumcenter_computed_ = true;
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("Failed to compute circumcenter: ") + e.what());
    }
}

double Tetrahedron::getVolume() const {
    if (is_infinite_) {
        return std::numeric_limits<double>::infinity();
    }
    
    // Volume = |det(b-a, c-a, d-a)| / 6
    Point3D v1 = vertices_[1] - vertices_[0];
    Point3D v2 = vertices_[2] - vertices_[0];
    Point3D v3 = vertices_[3] - vertices_[0];
    
    return std::abs(v1.dot(v2.cross(v3))) / 6.0;
}

bool Tetrahedron::hasPositiveOrientation() const {
    if (is_infinite_) {
        return false;
    }
    
    return GeometricPredicates::orient3D(
        vertices_[0], vertices_[1], vertices_[2], vertices_[3]
    ) > 0;
}

bool Tetrahedron::containsPoint(const Point3D& point) const {
    if (is_infinite_) {
        return false;
    }
    
    // Check if point is on same side of each face as the opposite vertex
    for (int i = 0; i < 4; ++i) {
        std::array<Point3D, 3> face = getFace(i);
        double orient_opposite = GeometricPredicates::orient3D(face[0], face[1], face[2], vertices_[i]);
        double orient_point = GeometricPredicates::orient3D(face[0], face[1], face[2], point);
        
        if (orient_opposite * orient_point < 0) {
            return false;
        }
    }
    
    return true;
}

bool Tetrahedron::isPointInCircumsphere(const Point3D& point) const {
    if (is_infinite_) {
        return false;
    }
    
    return GeometricPredicates::inSphere(
        vertices_[0], vertices_[1], vertices_[2], vertices_[3], point
    ) > 0;
}

bool Tetrahedron::hasVertex(const Point3D& point) const {
    for (int i = 0; i < (is_infinite_ ? 3 : 4); ++i) {
        if (vertices_[i] == point) {
            return true;
        }
    }
    return false;
}

int Tetrahedron::getVertexIndex(const Point3D& point) const {
    for (int i = 0; i < (is_infinite_ ? 3 : 4); ++i) {
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
    
    // Return the face opposite to vertex face_index
    std::array<Point3D, 3> face;
    int idx = 0;
    for (int i = 0; i < 4; ++i) {
        if (i != face_index) {
            face[idx++] = vertices_[i];
        }
    }
    return face;
}

Point3D Tetrahedron::getFaceNormal(int face_index) const {
    std::array<Point3D, 3> face = getFace(face_index);
    Point3D v1 = face[1] - face[0];
    Point3D v2 = face[2] - face[0];
    return v1.cross(v2).normalized();
}

Point3D Tetrahedron::getFaceCenter(int face_index) const {
    std::array<Point3D, 3> face = getFace(face_index);
    return (face[0] + face[1] + face[2]) / 3.0;
}

bool Tetrahedron::isValid() const {
    if (is_infinite_) {
        return true; // Infinite tetrahedra are always valid for our purposes
    }
    
    // Check that all vertices are different
    for (int i = 0; i < 4; ++i) {
        for (int j = i + 1; j < 4; ++j) {
            if (vertices_[i] == vertices_[j]) {
                return false;
            }
        }
    }
    
    // Check positive volume
    return getVolume() > EPSILON;
}

std::string Tetrahedron::toString() const {
    std::ostringstream oss;
    oss << "Tetrahedron " << id_;
    if (is_infinite_) {
        oss << " (infinite)";
    }
    oss << " with vertices: ";
    for (int i = 0; i < (is_infinite_ ? 3 : 4); ++i) {
        oss << vertices_[i];
        if (i < (is_infinite_ ? 2 : 3)) oss << ", ";
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
}

void VoronoiCell::addNeighbor(ID cell_id) {
    if (std::find(neighboring_cells_.begin(), neighboring_cells_.end(), cell_id) 
        == neighboring_cells_.end()) {
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
    
    Point3D centroid(0, 0, 0);
    for (const Point3D& vertex : vertices_) {
        centroid += vertex;
    }
    return centroid / static_cast<double>(vertices_.size());
}

bool VoronoiCell::containsPoint(const Point3D& /* point */) const {
    // TODO: Implement point-in-polyhedron test
    // This requires proper face ordering and normal computation
    return false;
}

double VoronoiCell::distanceToSite(const Point3D& point) const {
    return point.distance(site_);
}

void VoronoiCell::computeVolume() {
    if (is_infinite_ || faces_.empty() || vertices_.empty()) {
        volume_ = std::numeric_limits<double>::infinity();
        volume_computed_ = true;
        return;
    }
    
    // Compute volume using divergence theorem / surface integration
    // For now, use a simple approximation
    volume_ = 1.0; // Placeholder
    
    // TODO: Implement proper polyhedral volume calculation
    // This involves triangulating faces and summing tetrahedral volumes
    
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
// DelaunayTriangulation Implementation
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
    
    // Insert points incrementally using Bowyer-Watson algorithm
    for (const Point3D& point : points) {
        insertPointBowyerWatson(point);
    }
    
    // Remove bounding tetrahedron
    removeBoundingTetrahedron();
    
    is_initialized_ = true;
}

void DelaunayTriangulation::addPoint(const Point3D& point) {
    if (!is_initialized_) {
        throw std::runtime_error("Triangulation not initialized");
    }
    
    points_.push_back(point);
    insertPointBowyerWatson(point);
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

void DelaunayTriangulation::removeTetrahedron(TetrahedronID id) {
    active_tetrahedra_.erase(id);
}

void DelaunayTriangulation::removeBoundingTetrahedron() {
    // Remove tetrahedra that contain any bounding vertices
    std::vector<TetrahedronID> to_remove;
    
    for (const auto& tet : tetrahedra_) {
        if (!tet || active_tetrahedra_.find(tet->getId()) == active_tetrahedra_.end()) {
            continue;
        }
        
        bool contains_bounding = false;
        for (int i = 0; i < 4; ++i) {
            const Point3D& vertex = tet->getVertex(i);
            for (const Point3D& bounding : bounding_vertices_) {
                if (vertex == bounding) {
                    contains_bounding = true;
                    break;
                }
            }
            if (contains_bounding) break;
        }
        
        if (contains_bounding) {
            to_remove.push_back(tet->getId());
        }
    }
    
    for (TetrahedronID id : to_remove) {
        removeTetrahedron(id);
    }
}

// Core Bowyer-Watson algorithm implementation
void DelaunayTriangulation::insertPointBowyerWatson(const Point3D& point) {
    // Step 1: Find all tetrahedra in conflict with the new point
    std::vector<TetrahedronPtr> conflicting = findConflictingTetrahedra(point);
    
    if (conflicting.empty()) {
        return; // Point might be duplicate or outside
    }
    
    // Step 2: Collect boundary faces of the conflicting region
    std::vector<std::array<Point3D, 3>> boundary_faces;
    std::set<std::array<Point3D, 3>> face_set;
    
    for (const auto& tet : conflicting) {
        for (int face_idx = 0; face_idx < 4; ++face_idx) {
            auto neighbor = tet->getNeighbor(face_idx);
            
            // Check if this face is on the boundary (neighbor not in conflict)
            bool is_boundary = true;
            if (neighbor && active_tetrahedra_.count(neighbor->getId())) {
                for (const auto& conflict_tet : conflicting) {
                    if (neighbor->getId() == conflict_tet->getId()) {
                        is_boundary = false;
                        break;
                    }
                }
            }
            
            if (is_boundary) {
                std::array<Point3D, 3> face = tet->getFace(face_idx);
                // Sort face vertices for canonical representation
                std::sort(face.begin(), face.end(), [](const Point3D& a, const Point3D& b) {
                    if (std::abs(a.x - b.x) > EPSILON) return a.x < b.x;
                    if (std::abs(a.y - b.y) > EPSILON) return a.y < b.y;
                    return a.z < b.z;
                });
                
                if (face_set.find(face) == face_set.end()) {
                    face_set.insert(face);
                    boundary_faces.push_back(face);
                }
            }
        }
    }
    
    // Step 3: Remove conflicting tetrahedra
    removeConflictingTetrahedra(conflicting);
    
    // Step 4: Retriangulate the cavity
    retriangulateRemovedRegion(point, boundary_faces);
}

std::vector<DelaunayTriangulation::TetrahedronPtr> DelaunayTriangulation::findConflictingTetrahedra(const Point3D& point) const {
    std::vector<TetrahedronPtr> conflicting;
    
    for (const auto& tet : tetrahedra_) {
        if (!tet || active_tetrahedra_.find(tet->getId()) == active_tetrahedra_.end()) {
            continue;
        }
        
        if (tet->isPointInCircumsphere(point)) {
            conflicting.push_back(tet);
        }
    }
    
    return conflicting;
}

void DelaunayTriangulation::removeConflictingTetrahedra(const std::vector<TetrahedronPtr>& conflicting) {
    for (const auto& tet : conflicting) {
        if (tet && active_tetrahedra_.find(tet->getId()) != active_tetrahedra_.end()) {
            // Clear neighbor relationships
            for (int i = 0; i < 4; ++i) {
                auto neighbor = tet->getNeighbor(i);
                if (neighbor) {
                    int neighbor_face = neighbor->findNeighborIndex(tet);
                    if (neighbor_face >= 0) {
                        neighbor->setNeighbor(neighbor_face, nullptr);
                    }
                }
            }
            
            removeTetrahedron(tet->getId());
        }
    }
}

void DelaunayTriangulation::retriangulateRemovedRegion(const Point3D& point, 
                                                      const std::vector<std::array<Point3D, 3>>& boundary_faces) {
    // Create new tetrahedra by connecting the point to each boundary face
    for (const auto& face : boundary_faces) {
        try {
            TetrahedronPtr new_tet = createTetrahedron(point, face[0], face[1], face[2]);
            
            // Ensure positive orientation
            if (!new_tet->hasPositiveOrientation()) {
                new_tet = createTetrahedron(point, face[0], face[2], face[1]);
            }
            
            if (new_tet->hasPositiveOrientation()) {
                active_tetrahedra_.insert(new_tet->getId());
            }
        } catch (const std::exception&) {
            // Skip degenerate tetrahedra
            continue;
        }
    }
    
    // Update neighbor relationships
    updateNeighborships();
}

void DelaunayTriangulation::updateNeighborships() {
    // Clear all neighbor relationships
    for (auto& tet : tetrahedra_) {
        if (tet && active_tetrahedra_.count(tet->getId())) {
            for (int i = 0; i < 4; ++i) {
                tet->setNeighbor(i, nullptr);
            }
        }
    }
    
    // Rebuild neighbor relationships
    auto active_tets = getActiveTetrahedra();
    for (size_t i = 0; i < active_tets.size(); ++i) {
        for (size_t j = i + 1; j < active_tets.size(); ++j) {
            auto tet1 = active_tets[i];
            auto tet2 = active_tets[j];
            
            // Check if they share a face
            for (int face1 = 0; face1 < 4; ++face1) {
                for (int face2 = 0; face2 < 4; ++face2) {
                    auto f1 = tet1->getFace(face1);
                    auto f2 = tet2->getFace(face2);
                    
                    // Sort faces for comparison
                    std::sort(f1.begin(), f1.end(), [](const Point3D& a, const Point3D& b) {
                        if (std::abs(a.x - b.x) > EPSILON) return a.x < b.x;
                        if (std::abs(a.y - b.y) > EPSILON) return a.y < b.y;
                        return a.z < b.z;
                    });
                    std::sort(f2.begin(), f2.end(), [](const Point3D& a, const Point3D& b) {
                        if (std::abs(a.x - b.x) > EPSILON) return a.x < b.x;
                        if (std::abs(a.y - b.y) > EPSILON) return a.y < b.y;
                        return a.z < b.z;
                    });
                    
                    if (f1 == f2) {
                        linkTetrahedra(tet1, face1, tet2, face2);
                    }
                }
            }
        }
    }
}

void DelaunayTriangulation::linkTetrahedra(TetrahedronPtr tet1, int face1, 
                                          TetrahedronPtr tet2, int face2) {
    tet1->setNeighbor(face1, tet2);
    tet2->setNeighbor(face2, tet1);
}

// Geometric query methods
DelaunayTriangulation::TetrahedronPtr DelaunayTriangulation::findContainingTetrahedron(const Point3D& point) const {
    for (const auto& tet : tetrahedra_) {
        if (tet && active_tetrahedra_.count(tet->getId()) && tet->containsPoint(point)) {
            return tet;
        }
    }
    return nullptr;
}

std::vector<DelaunayTriangulation::TetrahedronPtr> DelaunayTriangulation::getConflictingTetrahedra(const Point3D& point) const {
    return findConflictingTetrahedra(point);
}

std::vector<DelaunayTriangulation::TetrahedronPtr> DelaunayTriangulation::getAdjacentTetrahedra(const Point3D& point) const {
    std::vector<TetrahedronPtr> adjacent;
    for (const auto& tet : tetrahedra_) {
        if (tet && active_tetrahedra_.count(tet->getId()) && tet->hasVertex(point)) {
            adjacent.push_back(tet);
        }
    }
    return adjacent;
}

void DelaunayTriangulation::validateTriangulation() const {
    // Check Delaunay property: no point should be inside any circumsphere
    for (const auto& tet : tetrahedra_) {
        if (!tet || !active_tetrahedra_.count(tet->getId())) continue;
        
        for (const Point3D& point : points_) {
            if (!tet->hasVertex(point) && tet->isPointInCircumsphere(point)) {
                throw std::runtime_error("Delaunay property violation detected");
            }
        }
    }
}

bool DelaunayTriangulation::isValid() const {
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
// VoronoiDiagram Implementation
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
    cells_.clear();
    
    if (!delaunay_->isValid()) {
        throw std::runtime_error("Cannot extract Voronoi diagram from invalid triangulation");
    }
    
    // Create cells for each site
    for (size_t i = 0; i < sites_.size(); ++i) {
        auto cell = std::make_shared<VoronoiCell>(i, sites_[i]);
        cells_.push_back(cell);
    }
    
    // For each site, find all adjacent tetrahedra and construct the Voronoi cell
    for (size_t i = 0; i < sites_.size(); ++i) {
        const Point3D& site = sites_[i];
        auto cell = cells_[i];
        
        // Find all tetrahedra that have this site as a vertex
        auto adjacent_tetrahedra = delaunay_->getAdjacentTetrahedra(site);
        
        if (!adjacent_tetrahedra.empty()) {
            constructCell(cell, adjacent_tetrahedra);
        }
    }
}

VoronoiDiagram::CellPtr VoronoiDiagram::createCell(const Point3D& site) {
    static CellID next_id = 0;
    return std::make_shared<VoronoiCell>(next_id++, site);
}

void VoronoiDiagram::constructCell(CellPtr cell, const std::vector<DelaunayTriangulation::TetrahedronPtr>& incident_tetrahedra) {
    std::vector<Point3D> voronoi_vertices;
    
    // Collect circumcenters of all incident tetrahedra
    for (const auto& tet : incident_tetrahedra) {
        if (tet && !tet->isInfinite()) {
            try {
                Point3D circumcenter = tet->getCircumcenter();
                voronoi_vertices.push_back(circumcenter);
            } catch (const std::exception&) {
                // Skip degenerate tetrahedra
                continue;
            }
        }
    }
    
    // Add vertices to the cell
    for (const Point3D& vertex : voronoi_vertices) {
        cell->addVertex(vertex);
    }
    
    // Construct faces by connecting appropriate vertices
    // This is a simplified implementation - a full implementation would
    // need to properly order vertices and construct faces
    if (voronoi_vertices.size() >= 4) {
        // Create a simple face as an example
        std::vector<int> face_indices;
        for (size_t i = 0; i < std::min(voronoi_vertices.size(), size_t(4)); ++i) {
            face_indices.push_back(static_cast<int>(i));
        }
        cell->addFace(face_indices);
    }
}

// Query methods implementation
VoronoiDiagram::CellPtr VoronoiDiagram::findCell(const Point3D& query_point) const {
    double min_distance = std::numeric_limits<double>::max();
    CellPtr closest_cell = nullptr;
    
    for (const auto& cell : cells_) {
        double distance = cell->distanceToSite(query_point);
        if (distance < min_distance) {
            min_distance = distance;
            closest_cell = cell;
        }
    }
    
    return closest_cell;
}

std::vector<VoronoiDiagram::CellPtr> VoronoiDiagram::getCellsInRegion(const Point3D& center, double radius) const {
    std::vector<CellPtr> cells_in_region;
    
    for (const auto& cell : cells_) {
        if (cell->distanceToSite(center) <= radius) {
            cells_in_region.push_back(cell);
        }
    }
    
    return cells_in_region;
}

std::vector<VoronoiDiagram::CellPtr> VoronoiDiagram::getNeighborCells(CellID cell_id) const {
    std::vector<CellPtr> neighbors;
    
    if (cell_id < cells_.size()) {
        const auto& cell = cells_[cell_id];
        for (CellID neighbor_id : cell->getNeighboringCells()) {
            if (neighbor_id < cells_.size()) {
                neighbors.push_back(cells_[neighbor_id]);
            }
        }
    }
    
    return neighbors;
}

double VoronoiDiagram::getTotalVolume() const {
    double total = 0.0;
    for (const auto& cell : cells_) {
        if (!cell->isInfinite()) {
            total += cell->getVolume();
        }
    }
    return total;
}

Point3D VoronoiDiagram::getCentroid() const {
    if (sites_.empty()) {
        return Point3D(0, 0, 0);
    }
    
    Point3D centroid(0, 0, 0);
    for (const Point3D& site : sites_) {
        centroid += site;
    }
    centroid /= static_cast<double>(sites_.size());
    
    return centroid;
}

std::pair<Point3D, Point3D> VoronoiDiagram::getBoundingBox() const {
    if (sites_.empty()) {
        return {Point3D(0, 0, 0), Point3D(0, 0, 0)};
    }
    
    Point3D min_pt = sites_[0];
    Point3D max_pt = sites_[0];
    
    for (const Point3D& site : sites_) {
        min_pt.x = std::min(min_pt.x, site.x);
        min_pt.y = std::min(min_pt.y, site.y);
        min_pt.z = std::min(min_pt.z, site.z);
        max_pt.x = std::max(max_pt.x, site.x);
        max_pt.y = std::max(max_pt.y, site.y);
        max_pt.z = std::max(max_pt.z, site.z);
    }
    
    return {min_pt, max_pt};
}

bool VoronoiDiagram::isValid() const {
    return is_computed_ && delaunay_->isValid() && !cells_.empty();
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