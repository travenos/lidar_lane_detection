#include "clustering.h"

#include "kdtree.h"

#include <cmath>
#include <functional>

namespace {

bool fits_width(PointsVector& cluster, float squared_max_cluster_width)
{
  if (std::isinf(squared_max_cluster_width)) {
    return true;
  }
  for (const auto& point : cluster) {
    for (const auto& other_point : cluster) {
      const float squared_width = square(point.x - other_point.x) + square(point.y - other_point.y);
      if (squared_width > squared_max_cluster_width) {
        return false;
      }
    }
  }
  return true;
}

}

namespace processing_logic {

std::vector<PointsVector> cluster(const PointsVector& points, float cluster_tolerance,
                                  int min_size, int max_size, float max_cluster_width)
{
  const float squared_max_cluster_width = square(max_cluster_width);
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
        if (fits_width(new_cluster, squared_max_cluster_width)) {
          clusters.emplace_back(std::move(new_cluster));
        }
      }
    }
  }
  return clusters;
}

}
