#include "point_extractor.h"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace {

// PARAMETERS
constexpr float EXTRACTION_MULTIPLIER = 1.5f;
constexpr float MIN_DIV = 1.f;
constexpr float ROI_RADIUS = 55.f;
constexpr float HEIGHT_TOLERANCE = 1.5f;

float get_median_intensity(PointsVector points)
{
  auto m = points.begin() + points.size() / 2;
  std::nth_element(points.begin(), m, points.end(),
                   [](const auto& first, const auto& second)
                    { return first.intensity < second.intensity; }
  );
  return points[points.size() / 2].intensity;
}

float get_median_height(PointsVector points)
{
  auto m = points.begin() + points.size() / 2;
  std::nth_element(points.begin(), m, points.end(),
                   [](const auto& first, const auto& second)
                    { return first.z < second.z; }
  );
  return points[points.size() / 2].z;
}

float get_deviation(const PointsVector& points, float mean)
{
  float div{};
  for (const auto& point: points) {
    div += square(point.intensity - mean);
  }
  if (!points.empty()) {
    div /= static_cast<float>(points.size());
    div = std::sqrt(div);
  }
  return div;
}

PlainPointXYZI get_mass_center(const PointsVector& points)
{
  PlainPointXYZI result{};
  for (const auto& point : points) {
    result.x += point.x;
    result.y += point.y;
    result.z += point.z;
    result.intensity += point.intensity;
  }
  if (!points.empty()) {
    result.x /= static_cast<float>(points.size());
    result.y /= static_cast<float>(points.size());
    result.z /= static_cast<float>(points.size());
    result.intensity /= static_cast<float>(points.size());
  }
  return result;
}

}

namespace processing_logic {

PointsVector extract_intensity_outliers(const PointsVector& points)
{
  PointsVector filtered;
  const auto median = get_median_intensity(points);
  const auto deviation = std::max(get_deviation(points, median), MIN_DIV);
  for (const auto& point : points) {
    if (point.intensity - median > EXTRACTION_MULTIPLIER * deviation) {
      filtered.push_back(point);
    }
  }
  return filtered;
}

PointsVector fuse_points(const std::vector<PointsVector>& points_vectors)
{
  PointsVector fused;
  const auto size_summator = [](std::size_t sum, const PointsVector& vec) { return sum + vec.size(); };
  std::size_t num_points = std::accumulate(points_vectors.begin(), points_vectors.end(), std::size_t{},
                                           size_summator);
  fused.reserve(num_points);
  for (const auto& points_vector: points_vectors) {
    std::copy(points_vector.begin(), points_vector.end(), std::back_inserter(fused));
  }
  return fused;
}

PointsVector get_reduced_point_cloud(const std::vector<PointsVector>& clusters)
{
  constexpr float ROI_RADIUS_SQUARED = square(ROI_RADIUS);
  PointsVector clustered_cloud;

  // Replace points by cluster centers and filter by ROI
  // Note: this part may be parallelized in the future
  for (const auto& cluster : clusters) {
    const auto mass_center = get_mass_center(cluster);
    if (square(mass_center.x) + square(mass_center.y) < ROI_RADIUS_SQUARED) {
      clustered_cloud.push_back(mass_center);
    }
  }

  // Filter by height
  PointsVector result_cloud;
  result_cloud.reserve(clustered_cloud.size());
  const float median_height = get_median_height(clustered_cloud);
  for (const auto& point : clustered_cloud) {
    if (std::fabs(point.z - median_height) < HEIGHT_TOLERANCE) {
      result_cloud.push_back(point);
    }
  }
  return result_cloud;
}

}
