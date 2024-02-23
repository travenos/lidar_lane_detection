#include <iostream>
#include <fstream>
#include <memory>
#include <string_view>

#if defined(__GNUC__) && __GNUC__ < 8
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif

// TODO!!! PCL
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
//TODO!!! refactor code style
enum CameraAngle
{
  XY, TopDown, Side, FPS
};
static void initCamera(CameraAngle setAngle, pcl::visualization::PCLVisualizer& viewer)
{

    viewer.setBackgroundColor (0, 0, 0);

    // set camera position and angle
    viewer.initCameraParameters();
    // distance away in meters
    int distance = 16;

    switch(setAngle)
    {
        case XY : viewer.setCameraPosition(-distance, -distance, distance, 1, 1, 0); break;
        case TopDown : viewer.setCameraPosition(0, 0, distance, 1, 0, 1); break;
        case Side : viewer.setCameraPosition(0, -distance, 0, 0, 0, 1); break;
        case FPS : viewer.setCameraPosition(-10, 0, 0, 0, 0, 1);
    }

    if(setAngle!=FPS) {
        viewer.addCoordinateSystem (1.0);
    }
}
// TODO!!! PCL

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

  std::unique_ptr<float[]> pointcloud_data{new float[plointcloud_size]};
  point_cloud_file.read(reinterpret_cast<char*>(pointcloud_data.get()), file_size);
  // TODO!!! PCL
  pcl::visualization::PCLVisualizer viewer{"3D Viewer"};
  initCamera(XY, viewer);
  pcl::PointCloud<pcl::PointXYZI>::Ptr cloud{boost::make_shared<pcl::PointCloud<pcl::PointXYZI>>()};

  cloud->reserve(plointcloud_size / POINT_SIZE);
  const auto* point_ptr = pointcloud_data.get();
  while (point_ptr < pointcloud_data.get() + plointcloud_size)
  {
    pcl::PointXYZI point;
    point.x = point_ptr[0];
    point.y = point_ptr[1];
    point.z = point_ptr[2];
    point.intensity = std::clamp(point_ptr[3] * 10.f, 0.f, 255.f);
    cloud->push_back(point);
    point_ptr += POINT_SIZE;
  }
  pcl::visualization::PointCloudColorHandlerGenericField<pcl::PointXYZI> intensity_distribution(cloud,"intensity");
  viewer.addPointCloud<pcl::PointXYZI>(cloud, intensity_distribution, "Cloud");

  viewer.spin();
  // TODO!!! PCL
  return 0;
}
