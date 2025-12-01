#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>
#include <chrono>

#include "utils/VoronoiBenchmark.h"
#include "geometry/Point3D.h"

using namespace Voronoi3D;

void printUsage(const char* program_name) {
    std::cout << "Usage: " << program_name << " [OPTIONS]" << std::endl;
    std::cout << std::endl;
    std::cout << "OPTIONS:" << std::endl;
    std::cout << "  --benchmark TYPE           Run benchmark (scaling, distribution, single)" << std::endl;
    std::cout << "  --points N                 Number of points to generate" << std::endl;
    std::cout << "  --distribution DIST        Point distribution type:" << std::endl;
    std::cout << "                             cube, sphere, gaussian, grid, pathological" << std::endl;
    std::cout << "  --min-points N             Minimum points for scaling benchmark" << std::endl;
    std::cout << "  --max-points N             Maximum points for scaling benchmark" << std::endl;
    std::cout << "  --steps N                  Number of steps in scaling benchmark" << std::endl;
    std::cout << "  --runs N                   Number of runs to average" << std::endl;
    std::cout << "  --output FILE              Output CSV file for results" << std::endl;
    std::cout << "  --seed N                   Random seed for reproducible results" << std::endl;
    std::cout << "  --verify-delaunay          Verify Delaunay property (slow for large sets)" << std::endl;
    std::cout << "  --help                     Show this help message" << std::endl;
    std::cout << std::endl;
    std::cout << "EXAMPLES:" << std::endl;
    std::cout << "  # Single benchmark with 1000 random points in cube" << std::endl;
    std::cout << "  " << program_name << " --benchmark single --points 1000 --distribution cube" << std::endl;
    std::cout << std::endl;
    std::cout << "  # Scaling benchmark from 100 to 5000 points, 10 steps" << std::endl;
    std::cout << "  " << program_name << " --benchmark scaling --min-points 100 --max-points 5000 --steps 10" << std::endl;
    std::cout << std::endl;
    std::cout << "  # Distribution comparison with 1000 points each" << std::endl;
    std::cout << "  " << program_name << " --benchmark distribution --points 1000 --output results.csv" << std::endl;
}

PointDistribution parseDistribution(const std::string& dist_str) {
    if (dist_str == "cube") return PointDistribution::UNIFORM_CUBE;
    if (dist_str == "sphere") return PointDistribution::UNIFORM_SPHERE;
    if (dist_str == "gaussian") return PointDistribution::GAUSSIAN_CLUSTER;
    if (dist_str == "grid") return PointDistribution::GRID_WITH_NOISE;
    if (dist_str == "pathological") return PointDistribution::PATHOLOGICAL;
    
    std::cerr << "Unknown distribution: " << dist_str << std::endl;
    std::cerr << "Valid options: cube, sphere, gaussian, grid, pathological" << std::endl;
    exit(1);
}

int main(int argc, char* argv[]) {
    // Default parameters
    std::string benchmark_type = "single";
    size_t points = 1000;
    PointDistribution distribution = PointDistribution::UNIFORM_CUBE;
    size_t min_points = 100;
    size_t max_points = 5000;
    size_t steps = 10;
    size_t runs = 3;
    std::string output_file = "";
    uint32_t seed = std::chrono::steady_clock::now().time_since_epoch().count();
    bool verify_delaunay = false;
    
    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg == "--help") {
            printUsage(argv[0]);
            return 0;
        } else if (arg == "--benchmark" && i + 1 < argc) {
            benchmark_type = argv[++i];
        } else if (arg == "--points" && i + 1 < argc) {
            points = std::stoul(argv[++i]);
        } else if (arg == "--distribution" && i + 1 < argc) {
            distribution = parseDistribution(argv[++i]);
        } else if (arg == "--min-points" && i + 1 < argc) {
            min_points = std::stoul(argv[++i]);
        } else if (arg == "--max-points" && i + 1 < argc) {
            max_points = std::stoul(argv[++i]);
        } else if (arg == "--steps" && i + 1 < argc) {
            steps = std::stoul(argv[++i]);
        } else if (arg == "--runs" && i + 1 < argc) {
            runs = std::stoul(argv[++i]);
        } else if (arg == "--output" && i + 1 < argc) {
            output_file = argv[++i];
        } else if (arg == "--seed" && i + 1 < argc) {
            seed = std::stoul(argv[++i]);
        } else if (arg == "--verify-delaunay") {
            verify_delaunay = true;
        } else {
            std::cerr << "Unknown argument: " << arg << std::endl;
            printUsage(argv[0]);
            return 1;
        }
    }
    
    // Validate parameters
    if (benchmark_type != "single" && benchmark_type != "scaling" && 
        benchmark_type != "distribution") {
        std::cerr << "Invalid benchmark type: " << benchmark_type << std::endl;
        std::cerr << "Valid options: single, scaling, distribution" << std::endl;
        return 1;
    }
    
    if (min_points >= max_points && benchmark_type == "scaling") {
        std::cerr << "min-points must be less than max-points for scaling benchmark" << std::endl;
        return 1;
    }
    
    if (points == 0 || runs == 0 || steps == 0) {
        std::cerr << "Points, runs, and steps must be greater than 0" << std::endl;
        return 1;
    }
    
    try {
        std::cout << "=== 3D Voronoi Diagram Benchmark ===" << std::endl;
        std::cout << "Benchmark type: " << benchmark_type << std::endl;
        std::cout << "Random seed: " << seed << std::endl;
        std::cout << std::endl;
        
        VoronoiBenchmark benchmark(seed);
        std::vector<BenchmarkResult> results;
        
        if (benchmark_type == "single") {
            std::cout << "Running single benchmark:" << std::endl;
            std::cout << "  Points: " << points << std::endl;
            std::cout << "  Distribution: " << VoronoiBenchmark::distributionToString(distribution) << std::endl;
            std::cout << "  Runs: " << runs << std::endl;
            std::cout << std::endl;
            
            // Run multiple times and average
            for (size_t run = 0; run < runs; ++run) {
                std::cout << "Run " << (run + 1) << "/" << runs << "..." << std::endl;
                auto test_points = benchmark.generatePoints(points, distribution);
                auto result = benchmark.runBenchmark(test_points, distribution, verify_delaunay);
                results.push_back(result);
            }
            
        } else if (benchmark_type == "scaling") {
            std::cout << "Running scaling benchmark:" << std::endl;
            std::cout << "  Point range: " << min_points << " to " << max_points << std::endl;
            std::cout << "  Steps: " << steps << std::endl;
            std::cout << "  Distribution: " << VoronoiBenchmark::distributionToString(distribution) << std::endl;
            std::cout << "  Runs per size: " << runs << std::endl;
            std::cout << std::endl;
            
            results = benchmark.runScalingBenchmark(min_points, max_points, steps, 
                                                   distribution, runs);
            
        } else if (benchmark_type == "distribution") {
            std::cout << "Running distribution comparison benchmark:" << std::endl;
            std::cout << "  Points per distribution: " << points << std::endl;
            std::cout << "  Runs per distribution: " << runs << std::endl;
            std::cout << std::endl;
            
            results = benchmark.runDistributionBenchmark(points, runs);
        }
        
        // Print results to console
        benchmark.printResults(results);
        
        // Export to CSV if requested
        if (!output_file.empty()) {
            benchmark.exportResultsToCSV(results, output_file);
        }
        
        std::cout << "Benchmark completed successfully!" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}