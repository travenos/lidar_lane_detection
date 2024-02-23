#include "visualizer.h"
#include "types.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string_view>
#include <map>

#if defined(__GNUC__) && __GNUC__ < 8
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif

static_assert (sizeof (float) == 4, "Float doesn't consist of 4 bytes");

constexpr std::size_t POINT_SIZE = 5u;

auto parse_by_channels(const std::vector<float>& pointcloud_data) {
  std::map<int, std::vector<PlainPointXYZI>> result;
  for (auto point_iter = pointcloud_data.begin(); point_iter !=  pointcloud_data.end(); point_iter += POINT_SIZE)
  {
    const auto channel = static_cast<int>(*(point_iter + 4));
    const PlainPointXYZI point{*(point_iter + 0), *(point_iter + 1), *(point_iter + 2), *(point_iter + 3)};
    result[channel].push_back(point);
  }
  return result;
}


int main(int argc, char* argv[])
{
  if (argc < 2) {
    std::cerr << "Not enough arguments" << std::endl;
    return 1;
  }
  fs::path data_dir{argv[1]};
  if (!fs::is_directory(data_dir)) {
    std::cerr << "Provided path is not a directory" << std::endl;
    return 1;
  }
  for (auto const& point_cloud_dir : fs::directory_iterator{data_dir})
  {
    const auto point_cloud_path = point_cloud_dir.path();
    if (!fs::is_regular_file(point_cloud_path)) {
      continue;
    }
    const auto file_size = fs::file_size(point_cloud_path);
    std::fstream point_cloud_file{point_cloud_path, std::ios::in | std::ios::binary};
    if (file_size % sizeof (float)) {
      std::cerr << point_cloud_path << " can not be interpreted as an array of floats" << std::endl;
      continue;
    }
    const std::size_t plointcloud_size = file_size / sizeof (float);
    if (plointcloud_size % POINT_SIZE != 0) {
      std::cerr << point_cloud_path << ": the point size is not " << POINT_SIZE << std::endl;
      continue;
    }

    std::vector<float> pointcloud_data(plointcloud_size);
    point_cloud_file.read(reinterpret_cast<char*>(pointcloud_data.data()), file_size);

    vis_utils::visualize_cloud(pointcloud_data, POINT_SIZE, point_cloud_path.filename().string());

    const auto channels_map = parse_by_channels(pointcloud_data);
    for(auto channel_iter = channels_map.rbegin(); channel_iter != channels_map.rend(); ++channel_iter) {
      vis_utils::visualize_cloud(channel_iter->second, std::to_string(channel_iter->first));
    }
  }

  return 0;
}
