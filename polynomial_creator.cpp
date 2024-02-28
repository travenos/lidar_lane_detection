#include "polynomial_creator.h"

#include <Eigen/Core>

#include <array>
#include <cmath>
#include <list>

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

  const std::list<T>& get_list() const { return data_; }

  std::list<T> obtain_list() { return std::move(data_); }

  template <typename U>
  bool push(U&& element)
  {
    bool inserted{false};
    for (auto iter = data_.begin(); iter != data_.end(); ++iter) {
      // If the new element is greater or equal to the element under iterator
      if (comparer_(*iter, element)) {
        data_.insert(iter, std::forward<U>(element));
        inserted = true;
        break;
      }
    }
    if (!inserted && data_.size() < max_size_) {
      data_.push_back(std::forward<U>(element));
    }
    if (data_.size() > max_size_) {
      data_.pop_back();
    }
    return inserted;
  }

  template <typename... Args>
  bool emplace(Args&&... args)
  {
    return push(T{std::forward<Args>(args)...});
  }

  void merge(const PriorityList<T>& other)
  {
    for (const auto& other_element : other.get_list()) {
      push(other);
    }
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
//  std::cout << best_pair.front().transpose() << std::endl;
//  std::cout << best_pair.back().transpose() << std::endl;

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
