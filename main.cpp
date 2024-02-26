#include "visualizer.h"
#include "point_extractor.h"
#include "clustering.h"
#include "polynomial_creator.h"
#include "main_direction_getter.h"
#include "types.h"

#include <iostream>
#include <fstream>
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
constexpr float DISTANCE_IN_CLUSTER1 = 0.2f;
constexpr int MIN_POINTS_PER_CLUSTER1 = 1;
constexpr int MAX_POINTS_PER_CLUSTER1 = 10;
constexpr float DISTANCE_IN_CLUSTER2 = 0.3f;
constexpr int MIN_POINTS_PER_CLUSTER2 = 1;
constexpr int MAX_POINTS_PER_CLUSTER2 = 50;

auto parse_by_channels(const std::vector<float>& pointcloud_data) {
  PointsMap result;
  for (auto point_iter = pointcloud_data.begin(); point_iter !=  pointcloud_data.end(); point_iter += POINT_SIZE)
  {
    const auto channel = static_cast<int>(*(point_iter + 4));
    const PlainPointXYZI point{*(point_iter + 0), *(point_iter + 1), *(point_iter + 2), *(point_iter + 3)};
    result[channel].push_back(point);
  }
  return result;
}


int main(int argc, char* argv[])
{
  if (argc < 2) {
    std::cerr << "Not enough arguments" << std::endl;
    return 1;
  }
  fs::path data_dir{argv[1]};
  if (!fs::is_directory(data_dir)) {
    std::cerr << "Provided path is not a directory" << std::endl;
    return 1;
  }
  vis_utils::Visualizer visualizer_main;
  vis_utils::Visualizer visualizer_preprocessed_for_pca;
  vis_utils::Visualizer visualizer_filtered;
  for (auto const& point_cloud_dir : fs::directory_iterator{data_dir})
  {
    const auto point_cloud_path = point_cloud_dir.path();
    if (!fs::is_regular_file(point_cloud_path)) {
      continue;
    }
    const auto file_size = fs::file_size(point_cloud_path);
    std::fstream point_cloud_file{point_cloud_path, std::ios::in | std::ios::binary};
    if (file_size % sizeof (float)) {
      std::cerr << point_cloud_path << " can not be interpreted as an array of floats" << std::endl;
      continue;
    }
    const std::size_t plointcloud_size = file_size / sizeof (float);
    if (plointcloud_size % POINT_SIZE != 0) {
      std::cerr << point_cloud_path << ": the point size is not " << POINT_SIZE << std::endl;
      continue;
    }

    std::vector<float> pointcloud_data(plointcloud_size);
    point_cloud_file.read(reinterpret_cast<char*>(pointcloud_data.data()), file_size);

    PointsVector all_points{};
    all_points.reserve(plointcloud_size / POINT_SIZE);
    for (auto point_iter = pointcloud_data.begin(); point_iter !=  pointcloud_data.end(); point_iter += POINT_SIZE)
    {
      const PlainPointXYZI point{*(point_iter + 0), *(point_iter + 1), *(point_iter + 2), *(point_iter + 3)};
      all_points.push_back(point);
    }

    visualizer_main.visualize_cloud(pointcloud_data, POINT_SIZE, point_cloud_path.filename().string());

    const auto channels_map = parse_by_channels(pointcloud_data);
    std::vector<PointsVector> points_vectors;
    points_vectors.reserve(channels_map.size());
    std::vector<PointsVector> clusters;
    ChanneledClusteredPointClouds channelled_clouds;
    for(auto channel_iter = channels_map.rbegin(); channel_iter != channels_map.rend(); ++channel_iter) {
//      visualizer_filtered.visualize_cloud(channel_iter->second, std::to_string(channel_iter->first));
      auto outliers = processing_logic::extract_intensity_outliers(channel_iter->second);
//      visualizer_filtered.visualize_cloud(outliers, std::to_string(channel_iter->first));

      auto new_clusters = processing_logic::cluster(outliers, DISTANCE_IN_CLUSTER1, MIN_POINTS_PER_CLUSTER1,
                                                    MAX_POINTS_PER_CLUSTER1);
      std::copy(new_clusters.begin(), new_clusters.end(), std::back_inserter(clusters));

      auto& new_beam = channelled_clouds.emplace_back();
      new_beam = std::move(new_clusters);

      points_vectors.push_back(std::move(outliers));
    }

    const auto fused_filtered_points = processing_logic::fuse_points(points_vectors);

    //TODO!!! PCA
//    auto all_outliers = processing_logic::extract_intensity_outliers(all_points);
//    auto all_clusters = processing_logic::cluster(all_outliers, 0.3f, 1, 200, 0.3f);

    auto all_clusters = processing_logic::cluster(fused_filtered_points, DISTANCE_IN_CLUSTER2, MIN_POINTS_PER_CLUSTER2,
                                                  MAX_POINTS_PER_CLUSTER2);

    auto fused_filtered_center_mass_points = processing_logic::prepare_cloud_for_pca({all_clusters});

    {
      // Remove outlier regions
      auto meta_clusters = processing_logic::cluster(fused_filtered_center_mass_points, 15.f, 1,
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
    WeightedPolynomialsVector polynom
        = {{WeightedPolynomial{{0.f, 0.f, main_direction.y / main_direction.x, 0.f}, 0.f}}};
    visualizer_preprocessed_for_pca.visualize_clusters(all_clusters, {fused_filtered_center_mass_points}, polynom, point_cloud_path.filename().string() + "_filtered");
//    visualizer_filtered.visualize_clusters(clusters, points_vectors, polynom, point_cloud_path.filename().string() + "_clustered");
    //TODO!!! PCA

//    visualizer_filtered.visualize_cloud(fused_filtered_points, point_cloud_path.filename().string() + "_filtered");
    // TODO!!! Polynomials
    // const auto polynomials = processing_logic::find_lines(channelled_clouds);
    // visualizer_filtered.visualize_clusters(clusters, points_vectors, polynomials, point_cloud_path.filename().string() + "_clustered");
  }

  return 0;
}
