#pragma once

#include "types.h"

namespace processing_logic {
std::vector<PointsVector> cluster(const PointsVector& points, float cluster_tolerance, int min_size, int max_size);
}
