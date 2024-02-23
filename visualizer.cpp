#include "visualizer.h"

#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>

namespace {
 void visualize(const pcl::PointCloud<pcl::PointXYZI>::Ptr& cloud, const std::string& name)
 {
   constexpr double CAMERA_HEIGHT = 75.;
   pcl::visualization::PCLVisualizer viewer{name};
   viewer.setBackgroundColor(0, 0, 0);
   pcl::visualization::PointCloudColorHandlerGenericField<pcl::PointXYZI> intensity_distribution(cloud,"intensity");
   viewer.addPointCloud<pcl::PointXYZI>(cloud, intensity_distribution, "Cloud");

   viewer.setCameraPosition(0, 0, CAMERA_HEIGHT, 0, 1, 0);
   viewer.addCoordinateSystem(1.0);
   viewer.spin(); // TODO!!! add logic to avoid multiple windows
 }

 float process_intensity(float original_intensity)
 {
   return std::clamp(original_intensity * 10.f, 0.f, 255.f);
 }
}

namespace vis_utils {

void visualize_cloud(const std::vector<float>& pointcloud_data, std::size_t point_size, const std::string& name)
{
  if (point_size < 3) {
    throw std::invalid_argument{"Point size should be at least 3"};
  }

  pcl::PointCloud<pcl::PointXYZI>::Ptr cloud{boost::make_shared<pcl::PointCloud<pcl::PointXYZI>>()};

  cloud->reserve(pointcloud_data.size() / point_size);
  for (auto point_iter = pointcloud_data.begin(); point_iter != pointcloud_data.end(); point_iter += point_size)
  {
    pcl::PointXYZI point;
    point.x = *(point_iter + 0);
    point.y = *(point_iter + 1);
    point.z = *(point_iter + 2);
    if (point_size >= 4) {
      point.intensity = process_intensity(*(point_iter + 3));
    } else {
      point.intensity = 0;
    }
    cloud->push_back(point);
  }
  visualize(cloud, name);
}

void visualize_cloud(const std::vector<PlainPointXYZI>& pointcloud_data, const std::string& name)
{
  pcl::PointCloud<pcl::PointXYZI>::Ptr cloud{boost::make_shared<pcl::PointCloud<pcl::PointXYZI>>()};

  cloud->reserve(pointcloud_data.size());
  for (auto point_iter = pointcloud_data.begin(); point_iter != pointcloud_data.end(); ++point_iter)
  {
    pcl::PointXYZI point;
    point.x = point_iter->x;
    point.y = point_iter->y;
    point.z = point_iter->z;
    point.intensity = process_intensity(point_iter->intensity);
    cloud->push_back(point);
  }
  visualize(cloud, name);
}

}
