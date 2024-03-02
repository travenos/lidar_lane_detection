#include "clustering.h"

#include "kdtree.h"

#include <cmath>
#include <functional>
#include <limits>

namespace {

// PARAMETERS
constexpr float DISTANCE_IN_CLUSTER = 0.4f;
constexpr int MIN_POINTS_PER_CLUSTER = 1;
constexpr int MAX_POINTS_PER_CLUSTER = 250;

constexpr float METACLUSTER_DISTANCE = 15.f;

std::vector<PointsVector> cluster(const PointsVector& points, float cluster_tolerance,
                                  int min_size, int max_size)
{
  KdTree tree{};
  // May be this part can be improved - add balancing to the filled tree
  for (std::size_t i{0}; i < points.size(); ++i)
  {
    tree.insert(points.at(i), static_cast<int>(i));
  }

  std::vector<bool> processed_points(points.size(), false);

  std::function<void(int, std::vector<int>&)> proximity_func =
      [&] (int point_id, std::vector<int>& cluster)
  {
    const auto& point = points[static_cast<std::size_t>(point_id)];

    processed_points[static_cast<std::size_t>(point_id)] = true;
    cluster.push_back(point_id);
    auto nearby_points = tree.search(point, cluster_tolerance);
    for (int near_point_id: nearby_points)
    {
      if (processed_points[static_cast<std::size_t>(near_point_id)] == false)
      {
        proximity_func(near_point_id, cluster);
      }
    }
  };

  std::vector<PointsVector> clusters;
  for (std::size_t point_id{0}; point_id < points.size(); ++point_id)
  {
    if (processed_points[point_id] == false)
    {
      std::vector<int> new_cluster_ids;
      proximity_func(static_cast<int>(point_id), new_cluster_ids);
      const int cluster_size{static_cast<int>(new_cluster_ids.size())};
      if (cluster_size >= min_size && cluster_size <= max_size)
      {
        PointsVector new_cluster;
        new_cluster.reserve(new_cluster_ids.size());
        for (int new_cluster_id: new_cluster_ids) {
          new_cluster.push_back(points.at(static_cast<std::size_t>(new_cluster_id)));
        }
        clusters.emplace_back(std::move(new_cluster));
      }
    }
  }
  return clusters;
}

}

namespace processing_logic {

std::vector<PointsVector> cluster_marking_line_points(const PointsVector& points)
{
  return cluster(points, DISTANCE_IN_CLUSTER, MIN_POINTS_PER_CLUSTER, MAX_POINTS_PER_CLUSTER);
}

PointsVector find_largest_meta_cluster(const PointsVector& points)
{
  auto meta_clusters = cluster(points, METACLUSTER_DISTANCE, 1, std::numeric_limits<int>::max());
  std::size_t max_cluster_id{0};
  for (std::size_t i{1}; i < meta_clusters.size(); ++i) {
    if (meta_clusters.at(i).size() > meta_clusters.at(max_cluster_id).size()) {
      max_cluster_id = i;
    }
  }
  return meta_clusters.at(max_cluster_id);
}

}
