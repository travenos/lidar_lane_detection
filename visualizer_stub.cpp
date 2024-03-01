#include "visualizer.h"

// This is a stub file, which is used in builds where visualization is turned off
// (-DWITH_VISUALIZATION=OFF)

namespace vis_utils {

class ViewerImpl {
};

Visualizer::Visualizer(float line_visualization_radius)
  : line_visualization_radius_{line_visualization_radius}
{
}

Visualizer::~Visualizer() = default;

void Visualizer::visualize_cloud(const PointsVector&, const PolynomialsVector&, const std::string&)
{
}

void Visualizer::visualize_clusters(const std::vector<PointsVector>&, const std::vector<PointsVector>&,
                          const PolynomialsVector&, const std::string&) //TODO!!! remove
{
}

}
