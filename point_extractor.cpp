#include "point_extractor.h"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace {

// PARAMETERS
constexpr int EXTRACTION_MULTIPLIER = 2;
constexpr float MIN_DIV = 1.;
constexpr bool USE_MEDIAN = true;

float square(float val) { return val * val; }

float get_mean_intensity(const PointsVector& points)
{
  float sum{0.f};
  for (const auto& point : points) {
    sum += point.intensity;
  }
  return sum / points.size();
}

float get_median_intensity(PointsVector points)
{
  auto m = points.begin() + points.size() / 2;
  std::nth_element(points.begin(), m, points.end(),
                   [](const auto& first, const auto& second)
                    { return first.intensity < second.intensity; }
  );
  return points[points.size() / 2].intensity;
}

float get_deviation(const PointsVector& points, float mean)
{
  float div{};
  for (const auto& point: points) {
    div += square(point.intensity - mean);
  }
  if (!points.empty()) {
    div /= points.size();
    div = std::sqrt(div);
  }
  return div;
}
}

namespace processing_logic {

PointsVector extract_intensity_outliers(const PointsVector& points)
{
  PointsVector filtered;
  const auto mean = USE_MEDIAN ? get_median_intensity(points) : get_mean_intensity(points);
  const auto deviation = std::max(get_deviation(points, mean), MIN_DIV);
  for (const auto& point : points) {
    if (point.intensity - mean > EXTRACTION_MULTIPLIER * deviation) {
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

}
