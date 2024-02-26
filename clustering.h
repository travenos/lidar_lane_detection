#pragma once

#include "types.h"

#include <limits>

namespace processing_logic {
std::vector<PointsVector> cluster(const PointsVector& points, float cluster_tolerance,
                                  int min_size, int max_size,
                                  float max_cluster_width = std::numeric_limits<float>::infinity());
}
