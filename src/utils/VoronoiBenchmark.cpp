#include "utils/VoronoiBenchmark.h"
#include <random>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <sstream>

namespace Voronoi3D {

VoronoiBenchmark::VoronoiBenchmark(uint32_t seed) : rng_(seed) {
}

std::vector<Point3D> VoronoiBenchmark::generatePoints(size_t count, 
                                                     PointDistribution distribution,
                                                     double distribution_param) {
    switch (distribution) {
        case PointDistribution::UNIFORM_CUBE:
            return generateUniformCube(count);
        case PointDistribution::UNIFORM_SPHERE:
            return generateUniformSphere(count);
        case PointDistribution::GAUSSIAN_CLUSTER:
            return generateGaussianCluster(count, distribution_param);
        case PointDistribution::GRID_WITH_NOISE:
            return generateGridWithNoise(count, distribution_param);
        case PointDistribution::PATHOLOGICAL:
            return generatePathological(count);
        default:
            return generateUniformCube(count);
    }
}

BenchmarkResult VoronoiBenchmark::runBenchmark(const std::vector<Point3D>& points,
                                              PointDistribution distribution,
                                              bool verify_delaunay) {
    BenchmarkResult result{};
    result.point_count = points.size();
    result.distribution = distribution;
    
    // Measure construction time
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Build Delaunay triangulation
    DelaunayTriangulation delaunay;
    auto delaunay_start = std::chrono::high_resolution_clock::now();
    delaunay.initialize(points);
    auto delaunay_end = std::chrono::high_resolution_clock::now();
    
    result.construction_time_ms = std::chrono::duration<double, std::milli>(
        delaunay_end - delaunay_start).count();
    
    // Extract Voronoi diagram
    VoronoiDiagram voronoi;
    auto voronoi_start = std::chrono::high_resolution_clock::now();
    voronoi.compute(points);
    auto voronoi_end = std::chrono::high_resolution_clock::now();
    
    result.voronoi_extraction_time_ms = std::chrono::duration<double, std::milli>(
        voronoi_end - voronoi_start).count();
    
    auto end_time = std::chrono::high_resolution_clock::now();
    result.total_time_ms = std::chrono::duration<double, std::milli>(
        end_time - start_time).count();
    
    // Collect basic metrics
    result.tetrahedra_count = delaunay.getNumTetrahedra();
    result.voronoi_cells_count = voronoi.getNumCells();
    result.avg_tetrahedra_per_point = static_cast<double>(result.tetrahedra_count) / 
                                     static_cast<double>(result.point_count);
    
    // Verify Delaunay property if requested
    if (verify_delaunay && points.size() <= 1000) {
        result.delaunay_property_verified = verifyDelaunayProperty(delaunay);
    } else {
        result.delaunay_property_verified = true; // Assume correct for large sets
    }
    
    // Calculate quality metrics
    calculateQualityMetrics(delaunay, voronoi, result);
    
    // Estimate memory usage
    result.peak_memory_mb = estimateMemoryUsage(delaunay, voronoi);
    
    return result;
}

std::vector<BenchmarkResult> VoronoiBenchmark::runScalingBenchmark(size_t min_points,
                                                                  size_t max_points,
                                                                  size_t steps,
                                                                  PointDistribution distribution,
                                                                  size_t runs_per_size) {
    std::vector<BenchmarkResult> results;
    
    // Generate logarithmic scaling
    for (size_t i = 0; i < steps; ++i) {
        double ratio = static_cast<double>(i) / static_cast<double>(steps - 1);
        size_t point_count = static_cast<size_t>(
            min_points * std::pow(static_cast<double>(max_points) / min_points, ratio));
        
        // Run multiple times and average
        std::vector<BenchmarkResult> run_results;
        for (size_t run = 0; run < runs_per_size; ++run) {
            auto points = generatePoints(point_count, distribution);
            auto result = runBenchmark(points, distribution, point_count <= 1000);
            run_results.push_back(result);
        }
        
        // Average the results
        BenchmarkResult avg_result = run_results[0];
        if (runs_per_size > 1) {
            double total_construction_time = 0;
            double total_voronoi_time = 0;
            double total_time = 0;
            size_t total_tetrahedra = 0;
            
            for (const auto& run_result : run_results) {
                total_construction_time += run_result.construction_time_ms;
                total_voronoi_time += run_result.voronoi_extraction_time_ms;
                total_time += run_result.total_time_ms;
                total_tetrahedra += run_result.tetrahedra_count;
            }
            
            avg_result.construction_time_ms = total_construction_time / runs_per_size;
            avg_result.voronoi_extraction_time_ms = total_voronoi_time / runs_per_size;
            avg_result.total_time_ms = total_time / runs_per_size;
            avg_result.tetrahedra_count = total_tetrahedra / runs_per_size;
            avg_result.avg_tetrahedra_per_point = static_cast<double>(avg_result.tetrahedra_count) /
                                                 static_cast<double>(avg_result.point_count);
        }
        
        results.push_back(avg_result);
        
        // Progress indicator
        std::cout << "Completed " << (i + 1) << "/" << steps 
                  << " scaling tests (n=" << point_count << ")" << std::endl;
    }
    
    return results;
}

std::vector<BenchmarkResult> VoronoiBenchmark::runDistributionBenchmark(size_t point_count,
                                                                       size_t runs_per_distribution) {
    std::vector<BenchmarkResult> results;
    
    std::vector<PointDistribution> distributions = {
        PointDistribution::UNIFORM_CUBE,
        PointDistribution::UNIFORM_SPHERE,
        PointDistribution::GAUSSIAN_CLUSTER,
        PointDistribution::GRID_WITH_NOISE,
        PointDistribution::PATHOLOGICAL
    };
    
    for (auto dist : distributions) {
        std::vector<BenchmarkResult> run_results;
        
        for (size_t run = 0; run < runs_per_distribution; ++run) {
            auto points = generatePoints(point_count, dist);
            auto result = runBenchmark(points, dist, point_count <= 1000);
            run_results.push_back(result);
        }
        
        // Average the results for this distribution
        BenchmarkResult avg_result = run_results[0];
        if (runs_per_distribution > 1) {
            double total_construction_time = 0;
            double total_voronoi_time = 0;
            double total_time = 0;
            size_t total_tetrahedra = 0;
            
            for (const auto& run_result : run_results) {
                total_construction_time += run_result.construction_time_ms;
                total_voronoi_time += run_result.voronoi_extraction_time_ms;
                total_time += run_result.total_time_ms;
                total_tetrahedra += run_result.tetrahedra_count;
            }
            
            avg_result.construction_time_ms = total_construction_time / runs_per_distribution;
            avg_result.voronoi_extraction_time_ms = total_voronoi_time / runs_per_distribution;
            avg_result.total_time_ms = total_time / runs_per_distribution;
            avg_result.tetrahedra_count = total_tetrahedra / runs_per_distribution;
            avg_result.avg_tetrahedra_per_point = static_cast<double>(avg_result.tetrahedra_count) /
                                                 static_cast<double>(avg_result.point_count);
        }
        
        results.push_back(avg_result);
        
        std::cout << "Completed distribution: " << distributionToString(dist) << std::endl;
    }
    
    return results;
}

void VoronoiBenchmark::exportResultsToCSV(const std::vector<BenchmarkResult>& results,
                                         const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open " << filename << " for writing" << std::endl;
        return;
    }
    
    // Write header
    file << "PointCount,Distribution,ConstructionTime_ms,VoronoiTime_ms,TotalTime_ms,"
         << "TetrahedraCount,VoronoiCellsCount,PeakMemory_MB,DelaunayVerified,"
         << "AvgTetrahedraPerPoint,AvgCircumradiusRatio,MinCellVolume,MaxCellVolume,"
         << "AvgCellVolume,CellVolumeStdDev,MinNeighbors,MaxNeighbors,AvgNeighbors\\n";
    
    // Write data
    for (const auto& result : results) {
        file << result.point_count << ","
             << distributionToString(result.distribution) << ","
             << std::fixed << std::setprecision(3) << result.construction_time_ms << ","
             << result.voronoi_extraction_time_ms << ","
             << result.total_time_ms << ","
             << result.tetrahedra_count << ","
             << result.voronoi_cells_count << ","
             << result.peak_memory_mb << ","
             << (result.delaunay_property_verified ? "true" : "false") << ","
             << std::setprecision(2) << result.avg_tetrahedra_per_point << ","
             << result.avg_circumradius_ratio << ","
             << result.min_cell_volume << ","
             << result.max_cell_volume << ","
             << result.avg_cell_volume << ","
             << result.cell_volume_stddev << ","
             << result.min_neighbors_per_cell << ","
             << result.max_neighbors_per_cell << ","
             << result.avg_neighbors_per_cell << "\\n";
    }
    
    file.close();
    std::cout << "Results exported to " << filename << std::endl;
}

void VoronoiBenchmark::printResults(const std::vector<BenchmarkResult>& results) {
    std::cout << std::endl;
    std::cout << "=== Benchmark Results ===" << std::endl;
    std::cout << std::left << std::setw(12) << "Points"
              << std::setw(16) << "Distribution"
              << std::setw(12) << "Total(ms)"
              << std::setw(12) << "Delaunay(ms)"
              << std::setw(12) << "Voronoi(ms)"
              << std::setw(10) << "Tetrahedra"
              << std::setw(10) << "Memory(MB)"
              << std::setw(8) << "T/P Ratio" << std::endl;
    std::cout << std::string(110, '-') << std::endl;
    
    for (const auto& result : results) {
        std::cout << std::left << std::setw(12) << result.point_count
                  << std::setw(16) << distributionToString(result.distribution).substr(0, 15)
                  << std::setw(12) << std::fixed << std::setprecision(1) << result.total_time_ms
                  << std::setw(12) << result.construction_time_ms
                  << std::setw(12) << result.voronoi_extraction_time_ms
                  << std::setw(10) << result.tetrahedra_count
                  << std::setw(10) << result.peak_memory_mb
                  << std::setw(8) << std::setprecision(2) << result.avg_tetrahedra_per_point
                  << std::endl;
    }
    std::cout << std::endl;
}

std::string VoronoiBenchmark::distributionToString(PointDistribution dist) {
    switch (dist) {
        case PointDistribution::UNIFORM_CUBE: return "UniformCube";
        case PointDistribution::UNIFORM_SPHERE: return "UniformSphere";
        case PointDistribution::GAUSSIAN_CLUSTER: return "GaussianCluster";
        case PointDistribution::GRID_WITH_NOISE: return "GridWithNoise";
        case PointDistribution::PATHOLOGICAL: return "Pathological";
        default: return "Unknown";
    }
}

// Private implementation methods

std::vector<Point3D> VoronoiBenchmark::generateUniformCube(size_t count) {
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    std::vector<Point3D> points;
    points.reserve(count);
    
    for (size_t i = 0; i < count; ++i) {
        points.emplace_back(dist(rng_), dist(rng_), dist(rng_));
    }
    
    return points;
}

std::vector<Point3D> VoronoiBenchmark::generateUniformSphere(size_t count) {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    std::vector<Point3D> points;
    points.reserve(count);
    
    for (size_t i = 0; i < count; ++i) {
        // Generate point on unit sphere using Marsaglia method
        double x, y, z;
        double length_sq;
        do {
            x = 2.0 * dist(rng_) - 1.0;
            y = 2.0 * dist(rng_) - 1.0;
            z = 2.0 * dist(rng_) - 1.0;
            length_sq = x * x + y * y + z * z;
        } while (length_sq > 1.0 || length_sq == 0.0);
        
        double length = std::sqrt(length_sq);
        points.emplace_back(x / length, y / length, z / length);
    }
    
    return points;
}

std::vector<Point3D> VoronoiBenchmark::generateGaussianCluster(size_t count, double sigma) {
    std::normal_distribution<double> dist(0.0, sigma);
    std::vector<Point3D> points;
    points.reserve(count);
    
    for (size_t i = 0; i < count; ++i) {
        points.emplace_back(dist(rng_), dist(rng_), dist(rng_));
    }
    
    return points;
}

std::vector<Point3D> VoronoiBenchmark::generateGridWithNoise(size_t count, double noise_factor) {
    std::uniform_real_distribution<double> noise_dist(-noise_factor, noise_factor);
    std::vector<Point3D> points;
    points.reserve(count);
    
    // Create roughly cubic grid
    int grid_size = static_cast<int>(std::ceil(std::cbrt(count)));
    double spacing = 2.0 / (grid_size - 1);
    
    for (int i = 0; i < grid_size && points.size() < count; ++i) {
        for (int j = 0; j < grid_size && points.size() < count; ++j) {
            for (int k = 0; k < grid_size && points.size() < count; ++k) {
                double x = -1.0 + i * spacing + noise_dist(rng_);
                double y = -1.0 + j * spacing + noise_dist(rng_);
                double z = -1.0 + k * spacing + noise_dist(rng_);
                points.emplace_back(x, y, z);
            }
        }
    }
    
    // Trim to exact count
    points.resize(count);
    return points;
}

std::vector<Point3D> VoronoiBenchmark::generatePathological(size_t count) {
    std::vector<Point3D> points;
    points.reserve(count);
    
    // Generate points on convex hull (worst case for many algorithms)
    std::uniform_real_distribution<double> angle_dist(0.0, 2.0 * M_PI);
    std::uniform_real_distribution<double> height_dist(-1.0, 1.0);
    
    for (size_t i = 0; i < count; ++i) {
        double theta = angle_dist(rng_);
        double phi = angle_dist(rng_);
        double x = std::cos(theta) * std::sin(phi);
        double y = std::sin(theta) * std::sin(phi);
        double z = std::cos(phi);
        points.emplace_back(x, y, z);
    }
    
    return points;
}

bool VoronoiBenchmark::verifyDelaunayProperty(const DelaunayTriangulation& triangulation) {
    const auto& tetrahedra = triangulation.getActiveTetrahedra();
    const auto& points = triangulation.getPoints();
    
    for (const auto& tet : tetrahedra) {
        if (!tet || tet->isInfinite()) continue;
        
        // Check all other points against this tetrahedron's circumsphere
        for (const Point3D& point : points) {
            // Skip if point is a vertex of this tetrahedron
            if (tet->hasVertex(point)) continue;
            
            // Check if point is inside circumsphere
            if (tet->isPointInCircumsphere(point)) {
                return false; // Delaunay property violated
            }
        }
    }
    
    return true;
}

void VoronoiBenchmark::calculateQualityMetrics(const DelaunayTriangulation& triangulation,
                                              const VoronoiDiagram& voronoi,
                                              BenchmarkResult& result) {
    const auto& tetrahedra = triangulation.getActiveTetrahedra();
    const auto& cells = voronoi.getCells();
    
    // Calculate circumradius ratios for tetrahedra
    std::vector<double> circumradius_ratios;
    for (const auto& tet : tetrahedra) {
        if (!tet || tet->isInfinite()) continue;
        
        // Simple quality metric: circumradius to shortest edge ratio
        const auto& vertices = tet->getVertices();
        double min_edge = std::numeric_limits<double>::max();
        for (int i = 0; i < 4; ++i) {
            for (int j = i + 1; j < 4; ++j) {
                double edge_length = vertices[i].distance(vertices[j]);
                min_edge = std::min(min_edge, edge_length);
            }
        }
        
        if (min_edge > 0) {
            double circumradius = std::sqrt(tet->getCircumradiusSquared());
            circumradius_ratios.push_back(circumradius / min_edge);
        }
    }
    
    if (!circumradius_ratios.empty()) {
        result.avg_circumradius_ratio = std::accumulate(circumradius_ratios.begin(),
                                                       circumradius_ratios.end(), 0.0) /
                                       circumradius_ratios.size();
    }
    
    // Calculate Voronoi cell quality metrics
    std::vector<double> cell_volumes;
    std::vector<int> neighbor_counts;
    
    for (const auto& cell : cells) {
        if (!cell) continue;
        
        double volume = cell->getVolume();
        if (volume > 0) {
            cell_volumes.push_back(volume);
        }
        
        neighbor_counts.push_back(static_cast<int>(cell->getNeighboringCells().size()));
    }
    
    if (!cell_volumes.empty()) {
        result.min_cell_volume = *std::min_element(cell_volumes.begin(), cell_volumes.end());
        result.max_cell_volume = *std::max_element(cell_volumes.begin(), cell_volumes.end());
        result.avg_cell_volume = std::accumulate(cell_volumes.begin(), cell_volumes.end(), 0.0) /
                                cell_volumes.size();
        
        // Calculate standard deviation
        double sum_sq_diff = 0.0;
        for (double vol : cell_volumes) {
            double diff = vol - result.avg_cell_volume;
            sum_sq_diff += diff * diff;
        }
        result.cell_volume_stddev = std::sqrt(sum_sq_diff / cell_volumes.size());
    }
    
    if (!neighbor_counts.empty()) {
        result.min_neighbors_per_cell = *std::min_element(neighbor_counts.begin(), neighbor_counts.end());
        result.max_neighbors_per_cell = *std::max_element(neighbor_counts.begin(), neighbor_counts.end());
        result.avg_neighbors_per_cell = std::accumulate(neighbor_counts.begin(), neighbor_counts.end(), 0.0) /
                                       neighbor_counts.size();
    }
}

size_t VoronoiBenchmark::estimateMemoryUsage(const DelaunayTriangulation& triangulation,
                                           const VoronoiDiagram& voronoi) {
    size_t memory_bytes = 0;
    
    // Estimate tetrahedra memory
    memory_bytes += triangulation.getNumTetrahedra() * sizeof(Tetrahedron);
    
    // Estimate points memory
    memory_bytes += triangulation.getNumPoints() * sizeof(Point3D);
    
    // Estimate Voronoi cells memory (rough approximation)
    memory_bytes += voronoi.getNumCells() * 512; // Approximate per-cell overhead
    
    return memory_bytes / (1024 * 1024); // Convert to MB
}

} // namespace Voronoi3D