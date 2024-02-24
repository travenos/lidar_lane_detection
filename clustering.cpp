#include "clustering.h"

#include "kdtree.h"

#include <functional>

namespace processing_logic {

std::vector<PointsVector> cluster(const PointsVector& points, float cluster_tolerance, int min_size)
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
      if (cluster_size >= min_size)
      {
        auto& new_cluster = clusters.emplace_back();
        new_cluster.reserve(new_cluster_ids.size());
        for (int new_cluster_id: new_cluster_ids) {
          new_cluster.push_back(points.at(static_cast<std::size_t>(new_cluster_id)));
        }
      }
    }
  }
  return clusters;
}

}
