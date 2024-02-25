#include "polynomial_creator.h"

#include <Eigen/Core>
#include <Eigen/LU>

#include <array>
#include <cassert>
#include <cmath>
#include <set>
#include <random>
#include <numeric>

namespace {

constexpr int RANDOM_SEED = 42;
constexpr int POINTS_PER_SEGMENT = 2;
constexpr int POLYNOMIAL_ORDER = 3;

constexpr float MAX_WEIGHT_PART = 0.15f;
constexpr float MAX_WEIGHT = 1.f;
constexpr float MIN_WEIGHT = 0.3f;
constexpr int RANSAC_ITERATIONS = 10000;
constexpr double GEOMETRIC_DISTRIBUTION_PARAM = 0.04;
constexpr std::size_t BEST_POLYNOMIALS_PER_PART_NUMBER = 10u;

// Assume, that z-coordinate is 0
using Cluster = std::vector<Eigen::Vector2f>;
using MetaCluster = std::vector<Cluster>;

struct Segments {
  enum Segment : std::size_t {
    UpRight, UpLeft, DownLeft, DownRight, Count
  };
};

constexpr auto square(float value)
{
  return value * value;
}

Eigen::Vector2f get_mass_center(const PointsVector& points)
{
  Eigen::Vector2f result{Eigen::Vector2f::Zero()};
  for (const auto& point : points) {
    result.x() += point.x;
    result.y() += point.y;
  }
  if (!points.empty()) {
    result /= static_cast<float>(points.size());
  }
  return result;
}

bool is_left_to_line(const Eigen::Vector2f& point, const Eigen::Vector2f& line_direction)
{
  return (line_direction.x() * point.y() - line_direction.y() * point.x()) >= 0;
}

auto segment(const ChanneledClusteredPointClouds& clouds,
             const Eigen::Vector2f& main_direction)
{
  Eigen::Vector2f traverse_direction{-main_direction.y(), main_direction.x()};
  std::array<MetaCluster, Segments::Count> result{{}};
  // Closer beams have higher values, reordering beams in reverse order
  for (const auto& beam: clouds) {
    for (auto& segment : result) {
      segment.emplace_back();
    }
    for (const auto& cluster : beam) {
      const auto mass_center = get_mass_center(cluster);
      const bool is_left_to_main_direction = is_left_to_line(mass_center, main_direction);
      const bool is_left_to_traverse_direction = is_left_to_line(mass_center, traverse_direction);

      std::size_t segment_id;
      if (is_left_to_main_direction) {
        segment_id = is_left_to_traverse_direction ? Segments::UpLeft : Segments::UpRight;
      } else {
        segment_id = is_left_to_traverse_direction ? Segments::DownLeft : Segments::DownRight;
      }

      auto& result_cluster = result[segment_id].back();
      result_cluster.emplace_back(mass_center);
    }
  }
  return result;
}

std::size_t calculate_not_empty_beams(const MetaCluster& beams)
{
  std::size_t result{};
  for (const auto& beam : beams) {
    result += static_cast<std::size_t>(!beam.empty());
  }
  return result;
}

WeightedPolynomialsVector select_results(const std::array<std::map<float, Eigen::Vector4f>, 2>& polynomials_map) {
  WeightedPolynomialsVector result;
  result.reserve(BEST_POLYNOMIALS_PER_PART_NUMBER * 2);
  for (const auto& road_part : polynomials_map) {
    std::size_t counter{0u};
    float min_loss{1.f};
    if (!road_part.empty()) {
      min_loss = road_part.begin()->first;
    }

    for (const auto& polynomial : road_part) {
      WeightedPolynomial result_entry;
      result_entry.coef0 = polynomial.second(0);
      result_entry.coef1 = polynomial.second(1);
      result_entry.coef2 = polynomial.second(2);
      result_entry.coef3 = polynomial.second(3);
      result_entry.weight = min_loss / polynomial.first;

      result.push_back(result_entry);
      ++counter;
      if (counter >= BEST_POLYNOMIALS_PER_PART_NUMBER) {
        break;
      }
    }
  }
  return result;
}

namespace ransac {

Eigen::Vector4f create_polynomial(const std::array<Eigen::Vector2f, POLYNOMIAL_ORDER + 1>& points)
{
  Eigen::Matrix4f A;
  Eigen::Vector4f b;
  for (int i{0}; i < A.rows(); ++i) {
    const auto x = points[static_cast<std::size_t>(i)].x();
    A(i, 1) = square(x);
    A(i, 0) = A(i, 1) * x; // x^3 = x^2 * x
    A(i, 2) = x;
    A(i, 3) = 1.;
  }
  for (int i{0}; i < b.rows(); ++i) {
    const auto y = points[static_cast<std::size_t>(i)].y();
    b(i) = y;
  }
  Eigen::Vector4f polynomial = A.lu().solve(b);
  return polynomial;
}

float get_weight(std::size_t beam_id, std::size_t beams_number)
{
  const auto part = static_cast<float>(beam_id) / static_cast<float>(beams_number);
  if (part > MAX_WEIGHT_PART) {
    const auto position = (part - MAX_WEIGHT_PART) / (1.f - MAX_WEIGHT_PART);
    return MAX_WEIGHT * (1.f - position) + MIN_WEIGHT * position;
  } else {
    return MAX_WEIGHT;
  }
}

float calculate_y_by_polynomial(const Eigen::Vector4f& polynomial, float x)
{
  const float x_square = square(x);
  return polynomial(0) * x_square * x + polynomial(1) * x_square + polynomial(2) * x + polynomial(3);
}

float evaluate_polynomial(const Eigen::Vector4f& polynomial,
                          const MetaCluster& left_cluster, const MetaCluster& right_cluster)
{
  float loss{0.};
  for (const auto* cluster : {&left_cluster, &right_cluster}) {
    for (std::size_t beam_id{0u}; beam_id < cluster->size(); ++beam_id) {
      if (cluster->at(beam_id).empty()) {
        continue;
      }
      float min_y_diff{INFINITY};
      for (const auto& point : cluster->at(beam_id)) {
        const float y_approx = calculate_y_by_polynomial(polynomial, point.x());
        const float diff = std::fabs(y_approx - point.y());
        if (diff < min_y_diff) {
          min_y_diff = diff;
        }
      }
      assert(!std::isinf(min_y_diff));

      loss += square(min_y_diff) * get_weight(beam_id, cluster->size());
    }
  }
  return loss;
}

auto perform_ransac(const std::array<MetaCluster, Segments::Count>& segmented_clusters)
{
  std::mt19937 gen{RANDOM_SEED};
  std::geometric_distribution<> segment_distrib{GEOMETRIC_DISTRIBUTION_PARAM};

  const std::array<std::array<std::size_t, 2>, 2> segment_ids
      = {{{{Segments::UpLeft, Segments::UpRight}}, {{Segments::DownLeft, Segments::DownRight}}}};
  // Fist iteration is up, the second one is down
  std::array<std::map<float, Eigen::Vector4f>, 2> result; // TODO!!! think about other data structure
  for (std::size_t part_id{0}; part_id < segment_ids.size(); ++part_id) {
    // Select two random beams from each segment. Select random cluster from each beam
    // Select beams
    for (int iter_counter{}; iter_counter < RANSAC_ITERATIONS; ++iter_counter)
    {
      std::set<std::size_t> left_beams_ids{};
      std::set<std::size_t> right_beams_ids{};
      while (left_beams_ids.size() < POINTS_PER_SEGMENT) {
        auto beam_id = static_cast<std::size_t>(std::round(segment_distrib(gen)));
        const auto segment_id = segment_ids[part_id][0];
        if (beam_id < segmented_clusters[segment_id].size()
            && !segmented_clusters[segment_id].at(beam_id).empty()) {
          left_beams_ids.emplace(beam_id);
        }
      }
      while (right_beams_ids.size() < POINTS_PER_SEGMENT) {
        auto beam_id = static_cast<std::size_t>(std::round(segment_distrib(gen)));
        const auto segment_id = segment_ids[part_id][1];
        if (beam_id < segmented_clusters[segment_id].size()
            && !segmented_clusters[segment_id].at(beam_id).empty()) {
          right_beams_ids.emplace(beam_id);
        }
      }
      // Select clusters per each beam
      std::array<Eigen::Vector2f, POLYNOMIAL_ORDER + 1> selected_points;
      std::size_t points_selected{0};
      for (const std::size_t beam_id : left_beams_ids) {
        const auto segment_id = segment_ids[part_id][0];
        std::uniform_int_distribution<std::size_t>
            point_distrib{0, segmented_clusters[segment_id].at(beam_id).size() - 1};
        const std::size_t point_id = point_distrib(gen);
        selected_points.at(points_selected) = segmented_clusters[segment_id].at(beam_id).at(point_id);
        ++points_selected;
      }
      for (const std::size_t beam_id : right_beams_ids) {
        const auto segment_id = segment_ids[part_id][1];
        std::uniform_int_distribution<std::size_t>
            point_distrib{0, segmented_clusters[segment_id].at(beam_id).size() - 1};
        const std::size_t point_id = point_distrib(gen);
        selected_points.at(points_selected) = segmented_clusters[segment_id].at(beam_id).at(point_id);
        ++points_selected;
      }

      // Create polynomial by selected points
      const auto polynomial = create_polynomial(selected_points);

      // Evaluate polynomial
      const auto left_id = segment_ids[part_id][0];
      const auto right_id = segment_ids[part_id][0];
      const auto score
          = evaluate_polynomial(polynomial, segmented_clusters[left_id], segmented_clusters[right_id]);

      // Record result
      result.at(part_id)[score] = polynomial;
    }
  }
  return result;
}

}

}

namespace processing_logic {

WeightedPolynomialsVector find_lines(const ChanneledClusteredPointClouds& clouds)
{
  const Eigen::Vector2f MAIN_DIRECTION{1.f, 0.f}; // TODO!!! calculate using PCA
  const auto segments = segment(clouds, MAIN_DIRECTION);
  for (std::size_t i{0u}; i < segments.size(); ++i) {
    if (calculate_not_empty_beams(segments.at(i)) < POINTS_PER_SEGMENT) {
      throw std::invalid_argument{"Not enough data in segment number " + std::to_string(i)};
    }
  }

  const auto results_map = ransac::perform_ransac(segments);
  const auto selected_results = select_results(results_map);
  return selected_results;
}

}
