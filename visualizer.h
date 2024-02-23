#pragma once

#include <cstddef>
#include <vector>

namespace vis_utils {
void visualize_cloud(const std::vector<float>& pointcloud_data, std::size_t point_size);
}
