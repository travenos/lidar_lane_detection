#include "visualizer.h"

#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>

namespace vis_utils {

constexpr double DISTANCE = 75.;

void visualize_cloud(const std::vector<float>& pointcloud_data, std::size_t point_size, const std::string& name)
{
  if (point_size < 3) {
    throw std::invalid_argument{"Point size should be at least 3"};
  }

  pcl::visualization::PCLVisualizer viewer{name};
  viewer.setBackgroundColor(0, 0, 0);

  pcl::PointCloud<pcl::PointXYZI>::Ptr cloud{boost::make_shared<pcl::PointCloud<pcl::PointXYZI>>()};

  cloud->reserve(pointcloud_data.size() / point_size);
  for (auto point_iter = pointcloud_data.begin(); point_iter !=  pointcloud_data.end(); point_iter += point_size)
  {
    pcl::PointXYZI point;
    point.x = *(point_iter + 0);
    point.y = *(point_iter + 1);
    point.z = *(point_iter + 2);
    if (point_size >= 4) {
      point.intensity = std::clamp(*(point_iter + 3) * 10.f, 0.f, 255.f);
    } else {
      point.intensity = 0;
    }
    cloud->push_back(point);
  }
  pcl::visualization::PointCloudColorHandlerGenericField<pcl::PointXYZI> intensity_distribution(cloud,"intensity");
  viewer.addPointCloud<pcl::PointXYZI>(cloud, intensity_distribution, "Cloud");

  viewer.setCameraPosition(0, 0, DISTANCE, 0, 1, 0);
  viewer.addCoordinateSystem(1.0);
  viewer.spin(); // TODO!!! add logic to avoid multiple windows
}

}
