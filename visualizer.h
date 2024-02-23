#pragma once

#include "types.h"

#include <cstddef>
#include <vector>
#include <string>

namespace vis_utils {
void visualize_cloud(const std::vector<float>& pointcloud_data, std::size_t point_size, const std::string& name);
void visualize_cloud(const std::vector<PlainPointXYZI>& pointcloud_data, const std::string& name);
}
