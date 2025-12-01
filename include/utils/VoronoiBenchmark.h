#pragma once

#include "geometry/Point3D.h"
#include "algorithms/VoronoiStructures.h"
#include <vector>
#include <string>
#include <chrono>
#include <fstream>
#include <random>

namespace Voronoi3D {

/**
 * @brief Point generation strategies for benchmarking
 */
enum class PointDistribution {
    UNIFORM_CUBE,       // Uniform random in unit cube
    UNIFORM_SPHERE,     // Uniform random on unit sphere surface
    GAUSSIAN_CLUSTER,   // Gaussian distribution around center
    GRID_WITH_NOISE,    // Regular grid with random noise
    PATHOLOGICAL       // Points designed to trigger worst-case behavior
};

/**
 * @brief Performance metrics for benchmark results
 */
struct BenchmarkResult {
    size_t point_count;
    PointDistribution distribution;
    double construction_time_ms;
    double voronoi_extraction_time_ms;
    double total_time_ms;
    size_t tetrahedra_count;
    size_t voronoi_cells_count;
    size_t peak_memory_mb;
    bool delaunay_property_verified;
    double avg_tetrahedra_per_point;
    double avg_circumradius_ratio;
    
    // Quality metrics
    double min_cell_volume;
    double max_cell_volume;
    double avg_cell_volume;
    double cell_volume_stddev;
    int min_neighbors_per_cell;
    int max_neighbors_per_cell;
    double avg_neighbors_per_cell;
};

/**
 * @brief Comprehensive benchmarking system for Voronoi diagrams
 */
class VoronoiBenchmark {
private:
    std::mt19937 rng_;
    
public:
    /**
     * @brief Constructor
     * @param seed Random seed for reproducible results
     */
    explicit VoronoiBenchmark(uint32_t seed = std::chrono::steady_clock::now().time_since_epoch().count());
    
    /**
     * @brief Generate points according to specified distribution
     * @param count Number of points to generate
     * @param distribution Point distribution strategy
     * @param distribution_param Optional parameter for distribution (e.g., sigma for Gaussian)
     * @return Vector of generated points
     */
    std::vector<Point3D> generatePoints(size_t count, 
                                       PointDistribution distribution,
                                       double distribution_param = 1.0);
    
    /**
     * @brief Run comprehensive benchmark on a point set
     * @param points Input points
     * @param distribution Distribution type (for result recording)
     * @param verify_delaunay Whether to verify Delaunay property (expensive for large sets)
     * @return Benchmark results
     */
    BenchmarkResult runBenchmark(const std::vector<Point3D>& points,
                                 PointDistribution distribution,
                                 bool verify_delaunay = true);
    
    /**
     * @brief Run scaling benchmark across different point counts
     * @param min_points Minimum number of points
     * @param max_points Maximum number of points
     * @param steps Number of steps in scaling test
     * @param distribution Point distribution to use
     * @param runs_per_size Number of runs to average per point count
     * @return Vector of benchmark results
     */
    std::vector<BenchmarkResult> runScalingBenchmark(size_t min_points,
                                                    size_t max_points,
                                                    size_t steps,
                                                    PointDistribution distribution,
                                                    size_t runs_per_size = 3);
    
    /**
     * @brief Run comparison benchmark across different distributions
     * @param point_count Number of points to use
     * @param runs_per_distribution Number of runs per distribution
     * @return Vector of benchmark results for each distribution
     */
    std::vector<BenchmarkResult> runDistributionBenchmark(size_t point_count,
                                                         size_t runs_per_distribution = 5);
    
    /**
     * @brief Export benchmark results to CSV file
     * @param results Vector of benchmark results
     * @param filename Output filename
     */
    void exportResultsToCSV(const std::vector<BenchmarkResult>& results,
                           const std::string& filename);
    
    /**
     * @brief Print benchmark results to console
     * @param results Vector of benchmark results
     */
    void printResults(const std::vector<BenchmarkResult>& results);
    
    /**
     * @brief Convert distribution enum to string
     */
    static std::string distributionToString(PointDistribution dist);
    
private:
    /**
     * @brief Generate uniform random points in unit cube
     */
    std::vector<Point3D> generateUniformCube(size_t count);
    
    /**
     * @brief Generate uniform random points on unit sphere surface
     */
    std::vector<Point3D> generateUniformSphere(size_t count);
    
    /**
     * @brief Generate Gaussian clustered points
     */
    std::vector<Point3D> generateGaussianCluster(size_t count, double sigma);
    
    /**
     * @brief Generate grid points with random noise
     */
    std::vector<Point3D> generateGridWithNoise(size_t count, double noise_factor);
    
    /**
     * @brief Generate pathological point set (e.g., points on convex hull)
     */
    std::vector<Point3D> generatePathological(size_t count);
    
    /**
     * @brief Verify Delaunay property for triangulation
     */
    bool verifyDelaunayProperty(const DelaunayTriangulation& triangulation);
    
    /**
     * @brief Calculate quality metrics for triangulation
     */
    void calculateQualityMetrics(const DelaunayTriangulation& triangulation,
                               const VoronoiDiagram& voronoi,
                               BenchmarkResult& result);
    
    /**
     * @brief Estimate memory usage
     */
    size_t estimateMemoryUsage(const DelaunayTriangulation& triangulation,
                             const VoronoiDiagram& voronoi);
};

} // namespace Voronoi3D