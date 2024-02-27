#include "polynomial_creator.h"

#include <Eigen/Core>
#include <Eigen/LU>

// TODO!!! clear includes
#include <array>
#include <cassert>
#include <cmath>
#include <set>
#include <random>
#include <numeric>
#include <list>

namespace {

constexpr int RANDOM_SEED = 43;
constexpr int POINTS_PER_SEGMENT = 1;

constexpr float NEAR_Y_DIST = 0.3f;
constexpr float MAX_WEIGHT_PART = 0.1f;
constexpr float MAX_WEIGHT = 1.f;
constexpr float MIN_WEIGHT = 0.1f;
constexpr int RANSAC_ITERATIONS = 1000;
constexpr double GEOMETRIC_DISTRIBUTION_PARAM = 0.05;
constexpr std::size_t BEST_POLYNOMIALS_PER_PART_NUMBER = 2u;

// Assume, that z-coordinate is 0
using Cluster = std::vector<Eigen::Vector2f>;
using MetaCluster = std::vector<Cluster>;

struct Segments {
  enum Segment : std::size_t {
    UpRight, UpLeft, DownLeft, DownRight, Count
  };
};

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

WeightedPolynomialsVector select_results(const std::array<std::map<float, Eigen::Vector2f>, 2>& polynomials_map) {
  WeightedPolynomialsVector result;
  std::map<float, bool> lines{{{2, false}, {-2, false}}};
  result.reserve(lines.size());
  for (const auto& road_part : polynomials_map) {
    std::size_t counter{0u};
    float min_loss{1.f};
    if (!road_part.empty()) {
      min_loss = road_part.begin()->first;
    }

    for (const auto& polynomial : road_part) {
      float b = polynomial.second(1);
      bool line_is_needed{false};
      for (auto& line: lines) {
        if (!line.second && std::fabs(line.first - b) < 1.) {
          line_is_needed = true;
          line.second = true;
          break;
        }
      }
      if (line_is_needed) {
        WeightedPolynomial result_entry;
        result_entry.coef0 = 0.f;
        result_entry.coef1 = 0.f;
        result_entry.coef2 = polynomial.second(0);
        result_entry.coef3 = polynomial.second(1);
        result_entry.weight = min_loss / polynomial.first;

        result.push_back(result_entry);
        ++counter;
        if (counter >= lines.size()) {
          break;
        }
      }
    }
  }
  return result;
}

namespace ransac {

Eigen::Vector2f create_polynomial(const std::pair<Eigen::Vector2f, Eigen::Vector2f>& points)
{
  Eigen::Matrix2f A;
  Eigen::Vector2f b;
  A(0, 0) = points.first.x();
  A(0, 1) = 1.;
  A(1, 0) = points.second.x();
  A(1, 1) = 1.;
  b(0) = points.first.y();
  b(1) = points.second.y();
  Eigen::Vector2f polynomial = A.lu().solve(b);
  return polynomial;
}

float get_weight(std::size_t beam_id, std::size_t beams_number)
{
  const auto part = static_cast<float>(beam_id) / static_cast<float>(beams_number);
  if (part > MAX_WEIGHT_PART) {
    return 0.f;
    const auto position = (part - MAX_WEIGHT_PART) / (1.f - MAX_WEIGHT_PART);
    return MAX_WEIGHT * (1.f - position) + MIN_WEIGHT * position;
  } else {
    return MAX_WEIGHT;
  }
}

float calculate_y_by_polynomial(const Eigen::Vector2f& polynomial, float x)
{
  return polynomial(0) * x + polynomial(1);
}

float evaluate_polynomial(const Eigen::Vector2f& polynomial,
                          const MetaCluster& left_cluster, const MetaCluster& right_cluster)
{
  float angle = std::atan(polynomial(0));
  std::array<float,2> points{{0., 0.}};
  std::size_t i{0u};
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

      points[i] += (min_y_diff < NEAR_Y_DIST) * get_weight(beam_id, cluster->size());
    }
    ++i;
  }
  //return -(points[0] + points[1] - std::fabs(points[0] - points[1]) - std::fabs(angle) * 50);
  return -(points[0] + points[1]);
}

auto perform_ransac(const std::array<MetaCluster, Segments::Count>& segmented_clusters)
{
  std::mt19937 gen{RANDOM_SEED};
  std::geometric_distribution<> segment_distrib{GEOMETRIC_DISTRIBUTION_PARAM};

  const std::array<std::array<std::size_t, 2>, 2> segment_ids
      = {{{{Segments::UpLeft, Segments::UpRight}}, {{Segments::DownLeft, Segments::DownRight}}}};
  // Fist iteration is up, the second one is down
  std::array<std::map<float, Eigen::Vector2f>, 2> result; // TODO!!! think about other data structure
  for (std::size_t part_id{0}; part_id < segment_ids.size(); ++part_id) {
    // Select two random beams from each segment. Select random cluster from each beam
    // Select beams
    for (int iter_counter{}; iter_counter < RANSAC_ITERATIONS; ++iter_counter)
    {
      std::optional<std::size_t> left_beam_id{};
      std::optional<std::size_t> right_beam_id{};
      while (!left_beam_id.has_value()) {
        auto beam_id = static_cast<std::size_t>(std::round(segment_distrib(gen)));
        const auto segment_id = segment_ids[part_id][0];
        if (beam_id < segmented_clusters[segment_id].size()
            && !segmented_clusters[segment_id].at(beam_id).empty()) {
          left_beam_id.emplace(beam_id);
        }
      }
      while (!right_beam_id.has_value()) {
        auto beam_id = static_cast<std::size_t>(std::round(segment_distrib(gen)));
        const auto segment_id = segment_ids[part_id][1];
        if (beam_id < segmented_clusters[segment_id].size()
            && !segmented_clusters[segment_id].at(beam_id).empty()) {
          right_beam_id.emplace(beam_id);
        }
      }
      // Select clusters per each beam
      std::pair<Eigen::Vector2f, Eigen::Vector2f> selected_points;
      {
        const auto segment_id = segment_ids[part_id][0];
        std::uniform_int_distribution<std::size_t>
            point_distrib{0, segmented_clusters[segment_id].at(left_beam_id.value()).size() - 1};
        const std::size_t point_id = point_distrib(gen);
        selected_points.first = segmented_clusters[segment_id].at(left_beam_id.value()).at(point_id);
      }
      {
        const auto segment_id = segment_ids[part_id][1];
        std::uniform_int_distribution<std::size_t>
            point_distrib{0, segmented_clusters[segment_id].at(right_beam_id.value()).size() - 1};
        const std::size_t point_id = point_distrib(gen);
        selected_points.second = segmented_clusters[segment_id].at(right_beam_id.value()).at(point_id);
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

//==================================================================================

namespace {

constexpr float deg_2_rad(float deg) { return deg * M_PIf32 / 180.f; }

constexpr float MAX_ANGLE_DEVIATION = deg_2_rad(10.f);
constexpr float MAIN_PART_RELATIVE_SIZE = 0.6f;
constexpr float MAIN_PART_SEGMENT_WIDTH = 7.f;
constexpr float NEAR_POINT_DIST = 0.1f;
constexpr float MIN_DIST_FROM_ORIGIN = 1.f;
constexpr float MAX_DIST_FROM_ORIGIN = 4.5f;
constexpr std::size_t MAX_BEST_ELEMENTS_COUNT = 400u;

const float MAX_ANGLE_DEVIATION_TAN = std::tan(MAX_ANGLE_DEVIATION);



template <typename T>
class PriorityList {
public:
  template <typename FCN>
  PriorityList(std::size_t max_size, FCN&& comparer)
    : comparer_{std::forward<FCN>(comparer)}
    , max_size_{max_size}
  {
  }

  const std::list<T>& get_data() const { return data_; }

  std::list<T> obtain_list() { return std::move(data_); }

  template <typename... Args>
  bool emplace(Args&&... args)
  {
    T element{std::forward<Args>(args)...};
    bool inserted{false};
    for (auto iter = data_.begin(); iter != data_.end(); ++iter) {
      // If the new element is greater than the element under iterator
      if (comparer_(*iter, element)) {
        data_.insert(iter, std::move(element));
        inserted = true;
        break;
      }
    }
    if (!inserted && data_.size() < max_size_) {
      data_.push_back(std::move(element));
    }
    if (data_.size() > max_size_) {
      data_.pop_back();
    }
    return inserted;
  }

private:
  std::list<T> data_;
  std::function<bool(const T&, const T&)> comparer_;
  std::size_t max_size_;
};

Eigen::MatrixXf point_cloud_to_matrix(const PointsVector& cloud)
{
  Eigen::MatrixXf points{static_cast<int>(cloud.size()), 2};
  for (std::size_t i{0}; i < cloud.size(); ++i) {
    const auto i_int = static_cast<int>(i);
    points.row(i_int).x() = cloud[i].x;
    points.row(i_int).y() = cloud[i].y;
  }
  return points;
}

void rotate_matrix(Eigen::MatrixXf& matrix, float angle)
{
  Eigen::Matrix2f rotation_matrix;
  const float angle_sin = std::sin(angle);
  const float angle_cos = std::cos(angle);
  rotation_matrix << angle_cos, angle_sin,
                    -angle_sin, angle_cos;

  matrix *= rotation_matrix;
}

void get_min_max(const Eigen::MatrixXf& matrix, Eigen::Vector2f& min, Eigen::Vector2f& max)
{
  min = matrix.colwise().minCoeff();
  max = matrix.colwise().maxCoeff();
}

std::array<Eigen::MatrixXf, 2> get_left_right_segments(const Eigen::MatrixXf& matrix,
                                                       const Eigen::Vector2f& left_down,
                                                       const Eigen::Vector2f& right_up)
{
  std::array<Eigen::MatrixXf, 2> result;
  // segment_id == 0 is left part, segment_id == 1 is right part
  for (std::size_t segment_id{0u}; segment_id < result.size(); ++segment_id) {
    std::vector<int> indices;
    for (int i{0}; i < matrix.rows(); ++i) {
      const auto& vec = matrix.row(i);
      if ((vec.transpose().array() > left_down.array()).all() && (vec.transpose().array() < right_up.array()).all()) {
        if (segment_id == 0 && vec.x() < 0.f) {
          indices.push_back(i);
        } else if (segment_id == 1 && vec.x() >= 0.f) {
          indices.push_back(i);
        }
      }
    }
    result[segment_id] = Eigen::MatrixXf{indices.size(), 2};
    int point_counter{0};
    for (int index : indices) {
      result[segment_id].row(point_counter) = matrix.row(index);
      ++point_counter;
    }
  }
  return result;
}

bool check_skew( const Eigen::Vector2f& left_pt, const Eigen::Vector2f& right_pt)
{
  const auto vec = right_pt - left_pt;
  return std::fabs(vec.y() / vec.x()) < MAX_ANGLE_DEVIATION_TAN;
}

Eigen::Vector3f get_line(const Eigen::Vector2f& left_pt, const Eigen::Vector2f& right_pt)
{
  Eigen::Vector3f line_coefs;
  line_coefs(0) = right_pt.y() - left_pt.y();
  line_coefs(1) = left_pt.x() - right_pt.x();
  line_coefs(2) = right_pt.x() * left_pt.y() - left_pt.x() * right_pt.y();

  if (line_coefs(0) < 0.f) {
    line_coefs *= -1;
  }

  return line_coefs;
}

float get_dist(const Eigen::Vector2f& point, const Eigen::Vector3f& line)
{
  float distance = std::fabs(Eigen::Vector3f{point.x(), point.y(), 1.}.dot(line)) / line.topRows<2>().norm();
  return distance;
}

bool point_fits_line(const Eigen::Vector2f& point, const Eigen::Vector3f& line)
{
  return get_dist(point, line) < NEAR_POINT_DIST;
}

float score_line(const Eigen::Vector3f& line, const std::array<Eigen::MatrixXf, 2>& point_cloud_segments)
{
  std::size_t fitting_points{0u};
  for (const auto& segment : point_cloud_segments) {
    for (int i{0}; i < segment.rows(); ++i) {
      fitting_points += point_fits_line(segment.row(i), line);
    }
  }

  return static_cast<float>(fitting_points);
}

auto find_best_lines(const std::array<Eigen::MatrixXf, 2>& point_cloud_segments)
{
  using StoredElement = std::pair<float, Eigen::Vector3f>;
  auto comparer = [](const StoredElement& first_val, const StoredElement& second_val) {
    return first_val.first < second_val.first;
  };
  PriorityList<StoredElement> result{MAX_BEST_ELEMENTS_COUNT, comparer};

  const auto& right_part = point_cloud_segments[1];
  const auto& left_part = point_cloud_segments[0];
  // Make all possible pairs between points from letf and right.
  // No need to ese RANSAC because of fiew points number.
  for (int right_point_id{0}; right_point_id < right_part.rows(); ++right_point_id) {
    const auto& right_point = right_part.row(right_point_id);
    for (int left_point_id{0}; left_point_id < left_part.rows(); ++left_point_id) {
      const auto& left_point = left_part.row(left_point_id);
      if (check_skew(left_point, right_point)) {
        const auto line = get_line(left_point, right_point);
        const float dist_from_origin = get_dist(Eigen::Vector2f::Zero(), line);
        if (dist_from_origin > MIN_DIST_FROM_ORIGIN && dist_from_origin < MAX_DIST_FROM_ORIGIN) {
          const float line_score = score_line(line, point_cloud_segments);
          result.emplace(line_score, line);
        }
      }
    }
  }
  return result.obtain_list();
}

auto choose_best_lines_pair(std::array<std::list<std::pair<float, Eigen::Vector3f>>, 2>& lines)
{
  if (lines[0].empty()) {
    throw std::invalid_argument{"No lines found in the upper part"};
  }
  if (lines[1].empty()) {
    throw std::invalid_argument{"No lines found in the down part"};
  }
  float min_cos{INFINITY};
  const auto* best_up_line = &(lines[0].front().second);
  const auto* best_down_line = &(lines[1].front().second);
  for (const auto& scored_up_line: lines[0]) {
    for (const auto& scored_down_line: lines[1]) {
      Eigen::Vector2f norm_vec_up = scored_up_line.second.topRows<2>().normalized();
      Eigen::Vector2f norm_vec_down = scored_down_line.second.topRows<2>().normalized();
      const float lines_cos = norm_vec_up.dot(norm_vec_down);
      if (lines_cos < min_cos) {
        best_up_line = &(scored_up_line.second);
        best_down_line = &(scored_down_line.second);
        min_cos = lines_cos;
      }
    }
  }
  return std::array{*best_up_line, *best_down_line};
}

}

#include <iostream> // TODO!!!

namespace processing_logic {

WeightedPolynomialsVector find_lines(const PointsVector& cloud, const Vec2D& main_direction)
{
  Eigen::MatrixXf matrix_cloud = point_cloud_to_matrix(cloud);
  // Rotate the point cloud to put it along the main adirection
  rotate_matrix(matrix_cloud, -std::atan2(main_direction.y, main_direction.x));
  // Find size of the cloud
  Eigen::Vector2f cloud_down_left;
  Eigen::Vector2f cloud_up_right;
  get_min_max(matrix_cloud, cloud_down_left, cloud_up_right);
  // Find length of the main part
  const float main_part_lenght = (cloud_up_right.x() - cloud_down_left.x()) * MAIN_PART_RELATIVE_SIZE;
  // segment_id == 0 is upper part, segment_id == 1 is down part
  std::array<std::list<std::pair<float, Eigen::Vector3f>>, 2> results;
  for (std::size_t segment_id{0u}; segment_id < 2u; ++segment_id) {
    Eigen::Vector2f left_down;
    Eigen::Vector2f right_up;
    const auto main_part_segment_width = MAIN_PART_SEGMENT_WIDTH;//MAX_DIST_FROM_ORIGIN + main_part_lenght / 2.f * MAX_ANGLE_DEVIATION_TAN;
    if (segment_id == 0u) { // upper part
      left_down = Eigen::Vector2f{-main_part_lenght / 2.f, 0.f};
      right_up = Eigen::Vector2f{main_part_lenght / 2.f, main_part_segment_width};
    } else { // down part
      left_down = Eigen::Vector2f{-main_part_lenght / 2.f, -main_part_segment_width};
      right_up = Eigen::Vector2f{main_part_lenght / 2.f, 0.f};
    }
    auto segments = get_left_right_segments(matrix_cloud, left_down, right_up);
    results[segment_id] = find_best_lines(segments);
  }
  auto best_pair = choose_best_lines_pair(results);

  //TODO!!! temporary
  std::cout << best_pair.front().transpose() << std::endl;
  std::cout << best_pair.back().transpose() << std::endl;

  WeightedPolynomialsVector result{};
  for (const auto& line : best_pair) {
    Eigen::Vector2f point{0, -line(2) / line(1)};
    float angle = std::atan2(main_direction.y, main_direction.x);
    Eigen::Matrix2f rotation_matrix;
    const float angle_sin = std::sin(angle);
    const float angle_cos = std::cos(angle);
    rotation_matrix << angle_cos, -angle_sin,
                      angle_sin, angle_cos;

    Eigen::Vector2f point2 = rotation_matrix * point;
    Eigen::Vector2f vector_normal2 = rotation_matrix * line.topRows<2>();
    float c2 = -point2.dot(vector_normal2);

    WeightedPolynomial polynomial;
    polynomial.coef2 = - vector_normal2(0) / vector_normal2(1);
    polynomial.coef3 = - c2 / vector_normal2(1);
    result.push_back(polynomial);
    result.push_back(polynomial);
  }
  return result;
  //TODO!!! temporary
}

}
