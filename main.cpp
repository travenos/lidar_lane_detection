#include "visualizer.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string_view>

#if defined(__GNUC__) && __GNUC__ < 8
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif

static_assert (sizeof (float) == 4, "Float consists not of 4 bytes");

constexpr std::size_t POINT_SIZE = 5u;

int main(int argc, char* argv[])
{
  if (argc < 2) {
    std::cerr << "Not enough arguments" << std::endl;
    return 1;
  }
  std::string_view point_cloud_path = argv[1];
  const auto file_size = fs::file_size(point_cloud_path);
  std::fstream point_cloud_file{point_cloud_path.data(), std::ios::in | std::ios::binary};
  if (file_size % sizeof (float)) {
    std::cerr << "File can not be interpreted as an array of floats" << std::endl;
    return 1;
  }
  const std::size_t plointcloud_size = file_size / sizeof (float);
  if (plointcloud_size % POINT_SIZE != 0) {
    std::cerr << "The point size is not " << POINT_SIZE << std::endl;
    return 1;
  }

  std::vector<float> pointcloud_data(plointcloud_size);
  point_cloud_file.read(reinterpret_cast<char*>(pointcloud_data.data()), file_size);

  vis_utils::visualize_cloud(pointcloud_data, POINT_SIZE);

  return 0;
}
