#include "visualizer.h"

#include <pcl/common/common.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>

namespace {

constexpr double CAMERA_HEIGHT = 75.;

struct Box
{
  float x_min;
  float y_min;
  float z_min;
  float x_max;
  float y_max;
  float z_max;
};

struct Color
{
  float r;
  float g;
  float b;
};

void render_box(pcl::visualization::PCLVisualizer& viewer, const Box& box, int id, const Color& color, float opacity = 1.f)
{
  std::clamp(opacity, 0.f, 1.f);

  std::string cube = "box" + std::to_string(id);
  viewer.addCube(box.x_min, box.x_max, box.y_min, box.y_max, box.z_min, box.z_max, color.r, color.g, color.b, cube);
  viewer.setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_REPRESENTATION, pcl::visualization::PCL_VISUALIZER_REPRESENTATION_WIREFRAME, cube);
  viewer.setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, color.r, color.g, color.b, cube);
  viewer.setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_OPACITY, opacity, cube);

  std::string cube_fill = "boxFill" + std::to_string(id);
  viewer.addCube(box.x_min, box.x_max, box.y_min, box.y_max, box.z_min, box.z_max, color.r, color.g, color.b, cube_fill);
  viewer.setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_REPRESENTATION, pcl::visualization::PCL_VISUALIZER_REPRESENTATION_SURFACE, cube_fill);
  viewer.setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, color.r, color.g, color.b, cube_fill);
  viewer.setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_OPACITY, opacity * 0.3, cube_fill);
}

void clear(pcl::visualization::PCLVisualizer& viewer)
{
  viewer.removeAllPointClouds();
  viewer.removeAllShapes();
}

void visualize(pcl::visualization::PCLVisualizer& viewer, const pcl::PointCloud<pcl::PointXYZI>::Ptr& cloud)
{
  clear(viewer);
  pcl::visualization::PointCloudColorHandlerGenericField<pcl::PointXYZI> intensity_distribution{cloud, "intensity"};
  viewer.addPointCloud<pcl::PointXYZI>(cloud, intensity_distribution, "Cloud");
  viewer.spin();
}

void visualize(pcl::visualization::PCLVisualizer& viewer,
               const pcl::PointCloud<pcl::PointXYZI>::Ptr& cloud,
               const std::vector<Box>& boxes)
{
  static const std::vector<Color> COLORS = {Color{1,0,0}, Color{1,1,0}, Color{0,0,1}};
  clear(viewer);
  pcl::visualization::PointCloudColorHandlerGenericField<pcl::PointXYZI> intensity_distribution{cloud, "intensity"};
  viewer.addPointCloud<pcl::PointXYZI>(cloud, intensity_distribution, "Cloud");
  for (std::size_t i{0}; i < boxes.size(); ++i) {
    render_box(viewer, boxes.at(i), static_cast<int>(i), COLORS[i % COLORS.size()]);
  }
  viewer.spin();
}

float process_intensity(float original_intensity)
{
  return std::clamp(original_intensity * 10.f, 0.f, 255.f);
}

Box get_bounding_box(const pcl::PointCloud<pcl::PointXYZI>::Ptr& cluster)
{
  // Find bounding box for one of the clusters
  pcl::PointXYZI min_point, max_point;
  pcl::getMinMax3D(*cluster, min_point, max_point);

  Box box;
  box.x_min = min_point.x;
  box.y_min = min_point.y;
  box.z_min = min_point.z;
  box.x_max = max_point.x;
  box.y_max = max_point.y;
  box.z_max = max_point.z;

  return box;
}

pcl::PointXYZI convert_to_pcl(const PlainPointXYZI& input_point)
{
  pcl::PointXYZI point;
  point.x = input_point.x;
  point.y = input_point.y;
  point.z = input_point.z;
  point.intensity = process_intensity(input_point.intensity);
  return point;
}

pcl::PointCloud<pcl::PointXYZI>::Ptr convert_to_cloud(const PointsVector& pointcloud_data)
{
  auto cloud{boost::make_shared<pcl::PointCloud<pcl::PointXYZI>>()};

  cloud->reserve(pointcloud_data.size());
  for (auto point_iter = pointcloud_data.begin(); point_iter != pointcloud_data.end(); ++point_iter) {
    cloud->push_back(convert_to_pcl(*point_iter));
  }
  return cloud;
}

}

namespace vis_utils {

Visualizer::Visualizer()
  : viewer_{std::make_unique<pcl::visualization::PCLVisualizer>("", false)}
{
}

Visualizer::~Visualizer() = default;

void Visualizer::visualize_cloud(const std::vector<float>& pointcloud_data, std::size_t point_size, const std::string& name)
{
  if (point_size < 3) {
    throw std::invalid_argument{"Point size should be at least 3"};
  }

  auto cloud{boost::make_shared<pcl::PointCloud<pcl::PointXYZI>>()};

  cloud->reserve(pointcloud_data.size() / point_size);
  for (auto point_iter = pointcloud_data.begin(); point_iter != pointcloud_data.end(); point_iter += point_size) {
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
  init_if_needed_();
  viewer_->setWindowName(name);
  visualize(*viewer_, cloud);
}

void Visualizer::visualize_cloud(const PointsVector& pointcloud_data, const std::string& name)
{
  auto cloud = convert_to_cloud(pointcloud_data);
  init_if_needed_();
  viewer_->setWindowName(name);
  visualize(*viewer_, cloud);
}

void Visualizer::visualize_clusters(const std::vector<PointsVector>& clusters, const std::string& name)
{
  auto cloud{boost::make_shared<pcl::PointCloud<pcl::PointXYZI>>()};
  std::vector<Box> boxes;
  boxes.reserve(clusters.size());
  for (const auto& cluster: clusters) {
    const auto new_cloud = convert_to_cloud(cluster);
    *cloud += *new_cloud;
    boxes.emplace_back(get_bounding_box(new_cloud));
  }
  init_if_needed_();
  viewer_->setWindowName(name);
  visualize(*viewer_, cloud, boxes);
}

void Visualizer::init_if_needed_()
{
  if (!viewer_initialized_) {
    viewer_->createInteractor();
    viewer_->setBackgroundColor(0, 0, 0);
    viewer_->setCameraPosition(0, 0, CAMERA_HEIGHT, 0, 1, 0);
    viewer_->addCoordinateSystem(1.0);
    viewer_initialized_ = true;
  }
}

}
