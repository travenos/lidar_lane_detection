#pragma once

#include "types.h"

#include <cstddef>
#include <vector>
#include <memory>
#include <string>

namespace pcl::visualization {
class PCLVisualizer;
}

namespace vis_utils {
class Visualizer final {
public:
  Visualizer();
  ~Visualizer();

  void visualize_cloud(const std::vector<float>& pointcloud_data, std::size_t point_size, const std::string& name);
  void visualize_cloud(const PointsVector& pointcloud_data, const std::string& name);
  void visualize_clusters(const std::vector<PointsVector>& clusters, const std::string& name);

private:
  void init_if_needed_();

  std::unique_ptr<pcl::visualization::PCLVisualizer> viewer_;
  bool viewer_initialized_{false};
};
}
