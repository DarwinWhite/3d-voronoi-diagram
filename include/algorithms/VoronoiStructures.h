#pragma once

#include "geometry/Point3D.h"
#include <vector>
#include <memory>
#include <unordered_set>

namespace Voronoi3D {

// Forward declarations
class DelaunayTriangulation;

/**
 * @brief Represents a tetrahedron in 3D Delaunay triangulation
 */
class Tetrahedron {
public:
    using ID = size_t;
    static constexpr ID INVALID_ID = std::numeric_limits<ID>::max();

private:
    ID id_;
    std::array<Point3D, 4> vertices_;
    std::array<std::weak_ptr<Tetrahedron>, 4> neighbors_;
    Point3D circumcenter_;
    double circumradius_squared_;
    bool is_infinite_;
    bool circumcenter_computed_;

public:
    /**
     * @brief Constructor for finite tetrahedron
     * @param id Unique identifier
     * @param a First vertex
     * @param b Second vertex
     * @param c Third vertex
     * @param d Fourth vertex
     */
    Tetrahedron(ID id, const Point3D& a, const Point3D& b, const Point3D& c, const Point3D& d);

    /**
     * @brief Constructor for infinite tetrahedron (on convex hull)
     * @param id Unique identifier
     * @param a First vertex
     * @param b Second vertex
     * @param c Third vertex
     */
    Tetrahedron(ID id, const Point3D& a, const Point3D& b, const Point3D& c);

    // Accessors
    ID getId() const { return id_; }
    const std::array<Point3D, 4>& getVertices() const { return vertices_; }
    const Point3D& getVertex(int index) const { return vertices_[index]; }
    bool isInfinite() const { return is_infinite_; }

    // Neighbor management
    void setNeighbor(int face_index, std::shared_ptr<Tetrahedron> neighbor);
    std::shared_ptr<Tetrahedron> getNeighbor(int face_index) const;
    int findNeighborIndex(std::shared_ptr<Tetrahedron> neighbor) const;

    // Geometric properties
    const Point3D& getCircumcenter();
    double getCircumradiusSquared();
    double getVolume() const;
    bool hasPositiveOrientation() const;

    // Geometric tests
    bool containsPoint(const Point3D& point) const;
    bool isPointInCircumsphere(const Point3D& point) const;
    bool hasVertex(const Point3D& point) const;
    int getVertexIndex(const Point3D& point) const;

    // Face operations
    std::array<Point3D, 3> getFace(int face_index) const;
    Point3D getFaceNormal(int face_index) const;
    Point3D getFaceCenter(int face_index) const;

    // Utility functions
    bool isValid() const;
    std::string toString() const;

private:
    void computeCircumcenter();
};

/**
 * @brief Represents a Voronoi cell (polyhedron)
 */
class VoronoiCell {
public:
    using ID = size_t;

private:
    ID id_;
    Point3D site_;
    std::vector<Point3D> vertices_;
    std::vector<std::vector<int>> faces_;
    std::vector<ID> neighboring_cells_;
    bool is_infinite_;
    double volume_;
    bool volume_computed_;

public:
    /**
     * @brief Constructor
     * @param id Unique identifier
     * @param site The generating point (Voronoi site)
     */
    VoronoiCell(ID id, const Point3D& site);

    // Accessors
    ID getId() const { return id_; }
    const Point3D& getSite() const { return site_; }
    const std::vector<Point3D>& getVertices() const { return vertices_; }
    const std::vector<std::vector<int>>& getFaces() const { return faces_; }
    bool isInfinite() const { return is_infinite_; }

    // Geometric properties
    double getVolume();
    Point3D getCentroid() const;
    std::vector<Point3D> getFaceNormals() const;

    // Cell construction (to be called by VoronoiDiagram)
    void addVertex(const Point3D& vertex);
    void addFace(const std::vector<int>& vertex_indices);
    void setInfinite(bool infinite) { is_infinite_ = infinite; }

    // Neighbor management
    void addNeighbor(ID cell_id);
    const std::vector<ID>& getNeighbors() const { return neighboring_cells_; }

    // Geometric tests
    bool containsPoint(const Point3D& point) const;
    double distanceToSite(const Point3D& point) const;

    // Utility functions
    bool isValid() const;
    void clear();
    std::string toString() const;

private:
    void computeVolume();
};

/**
 * @brief Main class for 3D Delaunay triangulation
 */
class DelaunayTriangulation {
public:
    using TetrahedronPtr = std::shared_ptr<Tetrahedron>;
    using TetrahedronID = Tetrahedron::ID;

private:
    std::vector<Point3D> points_;
    std::vector<TetrahedronPtr> tetrahedra_;
    std::unordered_set<TetrahedronID> active_tetrahedra_;
    TetrahedronID next_tetrahedron_id_;
    
    // Bounding tetrahedron for initialization
    std::array<Point3D, 4> bounding_vertices_;
    bool is_initialized_;

public:
    /**
     * @brief Constructor
     */
    DelaunayTriangulation();

    /**
     * @brief Destructor
     */
    ~DelaunayTriangulation() = default;

    // Main interface
    void initialize(const std::vector<Point3D>& points);
    void addPoint(const Point3D& point);
    void clear();

    // Accessors
    const std::vector<Point3D>& getPoints() const { return points_; }
    const std::vector<TetrahedronPtr>& getTetrahedra() const { return tetrahedra_; }
    std::vector<TetrahedronPtr> getActiveTetrahedra() const;
    size_t getNumPoints() const { return points_.size(); }
    size_t getNumTetrahedra() const { return active_tetrahedra_.size(); }

    // Geometric queries
    TetrahedronPtr findContainingTetrahedron(const Point3D& point) const;
    std::vector<TetrahedronPtr> getConflictingTetrahedra(const Point3D& point) const;
    std::vector<TetrahedronPtr> getAdjacentTetrahedra(const Point3D& point) const;

    // Validation and debugging
    bool isValid() const;
    void validateTriangulation() const;
    std::string getStatistics() const;

private:
    // Initialization helpers
    void createBoundingTetrahedron();
    void removeBoundingTetrahedron();

    // Bowyer-Watson algorithm implementation (to be implemented in Phase 3)
    void insertPointBowyerWatson(const Point3D& point);
    std::vector<TetrahedronPtr> findConflictingTetrahedra(const Point3D& point) const;
    void removeConflictingTetrahedra(const std::vector<TetrahedronPtr>& conflicting);
    void retriangulateRemovedRegion(const Point3D& point, 
                                   const std::vector<std::array<Point3D, 3>>& boundary_faces);

    // Tetrahedron management
    TetrahedronPtr createTetrahedron(const Point3D& a, const Point3D& b, 
                                   const Point3D& c, const Point3D& d);
    void removeTetrahedron(TetrahedronID id);
    TetrahedronID getNextTetrahedronId() { return next_tetrahedron_id_++; }

    // Neighbor management
    void updateNeighborships();
    void linkTetrahedra(TetrahedronPtr tet1, int face1, TetrahedronPtr tet2, int face2);
};

/**
 * @brief Main class for 3D Voronoi diagram
 */
class VoronoiDiagram {
public:
    using CellPtr = std::shared_ptr<VoronoiCell>;
    using CellID = VoronoiCell::ID;

private:
    std::vector<Point3D> sites_;
    std::vector<CellPtr> cells_;
    std::unique_ptr<DelaunayTriangulation> delaunay_;
    bool is_computed_;

public:
    /**
     * @brief Constructor
     */
    VoronoiDiagram();

    /**
     * @brief Destructor
     */
    ~VoronoiDiagram() = default;

    // Main interface
    void compute(const std::vector<Point3D>& sites);
    void addSite(const Point3D& site);
    void clear();

    // Accessors
    const std::vector<Point3D>& getSites() const { return sites_; }
    const std::vector<CellPtr>& getCells() const { return cells_; }
    const DelaunayTriangulation& getDelaunayTriangulation() const { return *delaunay_; }
    size_t getNumSites() const { return sites_.size(); }
    size_t getNumCells() const { return cells_.size(); }

    // Queries
    CellPtr findCell(const Point3D& query_point) const;
    std::vector<CellPtr> getCellsInRegion(const Point3D& center, double radius) const;
    std::vector<CellPtr> getNeighborCells(CellID cell_id) const;

    // Geometric properties
    double getTotalVolume() const;
    Point3D getCentroid() const;
    std::pair<Point3D, Point3D> getBoundingBox() const;

    // Validation and debugging
    bool isValid() const;
    std::string getStatistics() const;

private:
    // Dual extraction from Delaunay triangulation (to be implemented in Phase 3)
    void extractVoronoiDiagram();
    void computeVoronoiVertices();
    void computeVoronoiEdges();
    void computeVoronoiFaces();
    void handleInfiniteCells();

    // Cell construction helpers
    CellPtr createCell(const Point3D& site);
    void constructCell(CellPtr cell, const std::vector<DelaunayTriangulation::TetrahedronPtr>& incident_tetrahedra);
};

} // namespace Voronoi3D