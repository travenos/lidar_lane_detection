#pragma once

#include "types.h"

namespace processing_logic {

PointsVector extract_intensity_outliers(const PointsVector& points);
PointsVector fuse_points(const std::vector<PointsVector>& points_vectors);
}
