#include "visualizer.h"
#include "point_extractor.h"
#include "clustering.h"
#include "lines_detector.h"
#include "main_direction_getter.h"
#include "types.h"

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <vector>
#include <string_view>
#include <map>

#if defined(__GNUC__) && __GNUC__ < 8
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif

static_assert (sizeof (float) == 4, "Float doesn't consist of 4 bytes");

constexpr std::size_t POINT_SIZE = 5u;

constexpr float LINE_VISUALIZATION_RADIUS = 40.f;

namespace {

PointsMap read_point_cloud(const fs::path& point_cloud_path)
{
  PointsMap points_channels_map{};
  std::ifstream point_cloud_file{point_cloud_path, std::ios::in | std::ios::binary};
  while (point_cloud_file) {
    PlainPointXYZI point;
    point_cloud_file.read(reinterpret_cast<char*>(&point.x), sizeof (point.x));
    point_cloud_file.read(reinterpret_cast<char*>(&point.y), sizeof (point.y));
    point_cloud_file.read(reinterpret_cast<char*>(&point.z), sizeof (point.z));
    point_cloud_file.read(reinterpret_cast<char*>(&point.intensity), sizeof (point.intensity));
    float channel;
    point_cloud_file.read(reinterpret_cast<char*>(&channel), sizeof (channel));

    if (point_cloud_file.gcount() > 0)
    {
      points_channels_map[static_cast<int>(channel)].push_back(point);
    }
  }
  return points_channels_map;
}

void save_polynomyals(const fs::path& save_path, const PolynomialsVector& polynomials)
{
  std::ofstream output_file{save_path};
  output_file << std::scientific << std::setprecision(16);
  for (auto iter = polynomials.begin(); iter != polynomials.end(); ++iter) {
    if (iter != polynomials.begin()) {
      output_file << std::endl;
    }
    output_file << iter->coef0 << ';' << iter->coef1 << ';' << iter->coef2 << ';' << iter->coef3;
  }
}

}

int main(int argc, char* argv[])
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] <<
                 " [path to a directory with point clouds]" <<
                 " [path to an output directory to store results]" << std::endl;
    return EXIT_FAILURE;
  }
  fs::path data_dir{argv[1]};
  if (!fs::is_directory(data_dir)) {
    std::cerr << "Provided path " << data_dir << " is not a directory" << std::endl;
    return EXIT_FAILURE;
  }
  fs::path output_dir{argv[2]};
  if (!fs::exists(output_dir)) {
    fs::create_directories(output_dir);
  }
  if (!fs::is_directory(output_dir)) {
    std::cerr << "Provided path " << output_dir << " is not a directory" << std::endl;
    return EXIT_FAILURE;
  }

  vis_utils::Visualizer visualizer{LINE_VISUALIZATION_RADIUS};
  // Iterate over files with point cloud data
  for (auto const& point_cloud_file_entry : fs::directory_iterator{data_dir})
  {
    const auto point_cloud_path = point_cloud_file_entry.path();

    // Check the file
    if (!fs::is_regular_file(point_cloud_path)) {
      continue;
    }
    const auto file_size = fs::file_size(point_cloud_path);
    if (file_size % sizeof (float)) {
      std::cerr << point_cloud_path << " can not be interpreted as an array of floats" << std::endl;
      continue;
    }
    const std::size_t plointcloud_size = file_size / sizeof (float);
    if (plointcloud_size % POINT_SIZE != 0) {
      std::cerr << point_cloud_path << ": the point size is not " << POINT_SIZE << std::endl;
      continue;
    }

    // Read all points into a map
    const auto points_channels_map = read_point_cloud(point_cloud_path);

    // Copy all the points to a single array
    PointsVector all_points;
    all_points.reserve(plointcloud_size / POINT_SIZE);
    for (const auto& [_, point_vector] : points_channels_map) {
      std::copy(point_vector.begin(), point_vector.end(), std::back_inserter(all_points));
    }
    if (all_points.size() != plointcloud_size / POINT_SIZE) {
      std::cerr << point_cloud_path << ": unknown error occurred while reading point cloud" << std::endl;
      continue;
    }

    // Filter out possible marking line points by intensity
    std::vector<PointsVector> possible_marking_line_points;
    possible_marking_line_points.reserve(points_channels_map.size());
    for(auto channel_iter = points_channels_map.rbegin(); channel_iter != points_channels_map.rend(); ++channel_iter) {
      auto outliers = processing_logic::extract_intensity_outliers(channel_iter->second);
      possible_marking_line_points.push_back(std::move(outliers));
    }
    const auto fused_possible_marking_line_points = processing_logic::fuse_points(possible_marking_line_points);

    // Cluster possible marking lines points
    auto all_clusters = processing_logic::cluster_marking_line_points(fused_possible_marking_line_points);
    // Replace clusters by their mass centers
    auto possible_marking_lines_points_mass_centers = processing_logic::get_reduced_point_cloud(all_clusters);

    // Remove outlier regions by second clustering
    possible_marking_lines_points_mass_centers = processing_logic::find_largest_meta_cluster(possible_marking_lines_points_mass_centers);

    // Find main direction using PCA
    const auto main_direction = processing_logic::get_main_direction(possible_marking_lines_points_mass_centers);

    // Detect polynomials
    const auto polynomials = processing_logic::find_lines(possible_marking_lines_points_mass_centers, main_direction);

    // Save the result to a text file
    save_polynomyals(output_dir / point_cloud_path.stem().concat(".txt"), polynomials);

    // Visualize the result
    visualizer.visualize_cloud(all_points, polynomials, point_cloud_path.filename().string());
  }

  return EXIT_SUCCESS;
}
