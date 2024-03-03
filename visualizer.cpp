#include "visualizer.h"

#include <pcl/common/common.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>

namespace {

constexpr double CAMERA_HEIGHT = 60.;

struct Color
{
  float r;
  float g;
  float b;
};


void clear(pcl::visualization::PCLVisualizer& viewer)
{
  viewer.removeAllPointClouds();
  viewer.removeAllShapes();
}

void visualize(pcl::visualization::PCLVisualizer& viewer,
               const pcl::PointCloud<pcl::PointXYZI>::Ptr& cloud,
               const std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr>& curves)
{  
  constexpr Color LINE_COLOR{255.f, 0.f, 255.f};
  clear(viewer);
  pcl::visualization::PointCloudColorHandlerGenericField<pcl::PointXYZI> intensity_distribution{cloud, "intensity"};
  viewer.addPointCloud<pcl::PointXYZI>(cloud, intensity_distribution, "Cloud");

  for (std::size_t i{0}; i < curves.size(); ++i) {
    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> curve_color(curves.at(i), LINE_COLOR.r, LINE_COLOR.g, LINE_COLOR.b); // Red color
    viewer.addPointCloud<pcl::PointXYZ>(curves.at(i), curve_color, "PolynomialCurve" + std::to_string(i));
  }
  viewer.spin();
}

float process_intensity(float original_intensity)
{
  const float intensity = std::clamp(original_intensity * 10.f, 0.f, 255.f);
  return intensity;
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
  pcl::PointCloud<pcl::PointXYZI>::Ptr cloud{new pcl::PointCloud<pcl::PointXYZI>};

  cloud->reserve(pointcloud_data.size());
  for (auto point_iter = pointcloud_data.begin(); point_iter != pointcloud_data.end(); ++point_iter) {
    cloud->push_back(convert_to_pcl(*point_iter));
  }
  return cloud;
}

auto create_polynomials_representation(const PolynomialsVector& polynomials,
                                       float line_visualization_radius)
{
  constexpr float X_STEP = 0.01f;

  const float squared_radius = square(line_visualization_radius);

  std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> curves;
  curves.reserve(polynomials.size());

  for (const auto& polynomial : polynomials) {
    pcl::PointCloud<pcl::PointXYZ>::Ptr curve{new pcl::PointCloud<pcl::PointXYZ>};
    // Generate points along the polynomial curve
    float x = -line_visualization_radius;
    while (x <= line_visualization_radius) {
      const auto x_square = square(x);
      const auto a = polynomial.coef0;
      const auto b = polynomial.coef1;
      const auto c = polynomial.coef2;
      const auto d = polynomial.coef3;
      float y = a * x_square * x + b * x_square + c * x + d;
      x += X_STEP;
      if (square(x) + square(y) < squared_radius) {
        curve->points.emplace_back(x, y, 0.f);
      }
    }
    curves.emplace_back(std::move(curve));
  }
  return curves;
}

}

namespace vis_utils {

class ViewerImpl : public pcl::visualization::PCLVisualizer {
  using pcl::visualization::PCLVisualizer::PCLVisualizer;
};

Visualizer::Visualizer(float line_visualization_radius)
  : viewer_{std::make_unique<ViewerImpl>("", false)}
  , line_visualization_radius_{line_visualization_radius}
{
}

Visualizer::~Visualizer() = default;

void Visualizer::visualize_cloud(const PointsVector& pointcloud_data,
                                 const PolynomialsVector& polynomials,
                                 const std::string& name)
{
  auto cloud = convert_to_cloud(pointcloud_data);
  init_if_needed_();
  const auto curves = create_polynomials_representation(polynomials, line_visualization_radius_);
  viewer_->setWindowName(name);
  visualize(*viewer_, cloud, curves);
}

void Visualizer::init_if_needed_()
{
  if (!viewer_initialized_) {
    viewer_->createInteractor();
    viewer_->setBackgroundColor(0, 0, 0);
    viewer_->addCoordinateSystem(1.0);
    viewer_->setCameraPosition(0, 0, CAMERA_HEIGHT, 0, 1, 0);
    viewer_initialized_ = true;
  }
}

}
