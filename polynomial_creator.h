#pragma once

#include "types.h"

namespace processing_logic {
WeightedPolynomialsVector find_lines(const ChanneledClusteredPointClouds& clouds);
WeightedPolynomialsVector find_lines(const PointsVector& cloud, const Vec2D& main_direction);
}
