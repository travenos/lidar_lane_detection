#pragma once

#include "types.h"

#include <cstddef>
#include <vector>
#include <memory>
#include <string>

namespace vis_utils {

class ViewerImpl;

class Visualizer final {
public:
  explicit Visualizer(float line_visualization_radius);
  ~Visualizer();

  void visualize_cloud(const PointsVector& pointcloud_data, const PolynomialsVector& polynomials, const std::string& name);
  void visualize_clusters(const std::vector<PointsVector>& clusters, const std::vector<PointsVector>& all_points,
                          const PolynomialsVector& polynomials, const std::string& name); //TODO!!! remove

private:
  void init_if_needed_();

  std::unique_ptr<ViewerImpl> viewer_;
  float line_visualization_radius_;
  bool viewer_initialized_{false};
};
}
