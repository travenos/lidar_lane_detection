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

// PARAMETERS
constexpr float DISTANCE_IN_CLUSTER2 = 0.4f;
constexpr int MIN_POINTS_PER_CLUSTER2 = 1;
constexpr int MAX_POINTS_PER_CLUSTER2 = 250;

constexpr float METALCUSTER_DISTANCE = 15.f;

namespace {
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

  vis_utils::Visualizer visualizer_main;
  vis_utils::Visualizer visualizer_preprocessed_for_pca;
  vis_utils::Visualizer visualizer_filtered;
  for (auto const& point_cloud_file_entry : fs::directory_iterator{data_dir})
  {
    const auto point_cloud_path = point_cloud_file_entry.path();
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


    // Reading all points into a map
    PointsMap points_channels_map{};
    {
      std::fstream point_cloud_file{point_cloud_path, std::ios::in | std::ios::binary};
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
    }

    PointsVector all_points;
    all_points.reserve(plointcloud_size / POINT_SIZE);
    for (const auto& [_, point_vector] : points_channels_map) {
      std::copy(point_vector.begin(), point_vector.end(), std::back_inserter(all_points));
    }

    assert(all_points.size() == plointcloud_size / POINT_SIZE);

    std::vector<PointsVector> points_vectors;
    points_vectors.reserve(points_channels_map.size());
    std::vector<PointsVector> clusters;
    ChanneledClusteredPointClouds channelled_clouds;
    for(auto channel_iter = points_channels_map.rbegin(); channel_iter != points_channels_map.rend(); ++channel_iter) {
//      visualizer_filtered.visualize_cloud(channel_iter->second, std::to_string(channel_iter->first));
      auto outliers = processing_logic::extract_intensity_outliers(channel_iter->second);
//      visualizer_filtered.visualize_cloud(outliers, std::to_string(channel_iter->first));

//      auto new_clusters = processing_logic::cluster(outliers, DISTANCE_IN_CLUSTER1, MIN_POINTS_PER_CLUSTER1,
//                                                    MAX_POINTS_PER_CLUSTER1);
//      std::copy(new_clusters.begin(), new_clusters.end(), std::back_inserter(clusters));

//      auto& new_beam = channelled_clouds.emplace_back();
//      new_beam = std::move(new_clusters);

      points_vectors.push_back(std::move(outliers));
    }

    const auto fused_filtered_points = processing_logic::fuse_points(points_vectors);

    auto all_clusters = processing_logic::cluster(fused_filtered_points, DISTANCE_IN_CLUSTER2, MIN_POINTS_PER_CLUSTER2,
                                                  MAX_POINTS_PER_CLUSTER2);

    auto fused_filtered_center_mass_points = processing_logic::prepare_cloud_for_pca({all_clusters});

    {
      // Remove outlier regions
      auto meta_clusters = processing_logic::cluster(fused_filtered_center_mass_points, METALCUSTER_DISTANCE, 1,
                                                     std::numeric_limits<int>::max(),
                                                     std::numeric_limits<float>::infinity());
      std::size_t max_cluster_id{0};
      for (std::size_t i{1}; i < meta_clusters.size(); ++i) {
        if (meta_clusters.at(i).size() > meta_clusters.at(max_cluster_id).size()) {
          max_cluster_id = i;
        }
      }
      fused_filtered_center_mass_points = std::move(meta_clusters.at(max_cluster_id));
    }
    const auto main_direction = processing_logic::get_main_direction(fused_filtered_center_mass_points);
//    PolynomialsVector polynomials
//        = {{Polynomial{0.f, 0.f, main_direction.y / main_direction.x, 0.f}}};
    const auto polynomials = processing_logic::find_lines(fused_filtered_center_mass_points, main_direction);
    save_polynomyals(output_dir / point_cloud_path.stem().concat(".txt"), polynomials);

    // Old main logic!!!
    // visualizer_preprocessed_for_pca.visualize_clusters({}, {fused_filtered_center_mass_points}, polynomials, point_cloud_path.filename().string() + "_filtered");
    visualizer_main.visualize_clusters({}, {all_points}, polynomials, point_cloud_path.filename().string());

    //    visualizer_filtered.visualize_clusters(clusters, points_vectors, polynom, point_cloud_path.filename().string() + "_clustered");
    //TODO!!! PCA

//    visualizer_filtered.visualize_cloud(fused_filtered_points, point_cloud_path.filename().string() + "_filtered");
    // TODO!!! Polynomials
    // const auto polynomials = processing_logic::find_lines(channelled_clouds);
    // visualizer_filtered.visualize_clusters(clusters, points_vectors, polynomials, point_cloud_path.filename().string() + "_clustered");
  }

  return EXIT_SUCCESS;
}
