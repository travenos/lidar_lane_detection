#pragma once

#include "types.h"

namespace processing_logic {

std::vector<PointsVector> cluster_marking_line_points(const PointsVector& points);
PointsVector find_largest_meta_cluster(const PointsVector& points);

}
