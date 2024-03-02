#include "lines_detector.h"

#include <Eigen/Core>
#include <Eigen/Cholesky>

#include <array>
#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <list>
#include <vector>

namespace {

constexpr float deg_2_rad(float deg) { return deg * M_PIf32 / 180.f; }

constexpr float MAX_ANGLE_DEVIATION = deg_2_rad(10.f);
constexpr float MAIN_PART_RELATIVE_SIZE = 0.6f;
constexpr float MAIN_PART_SEGMENT_WIDTH = 6.f;
constexpr float NEAR_POINT_DIST_IN_SCORING = 0.1f;
constexpr float MIN_DIST_FROM_ORIGIN = 1.f;
constexpr float MAX_DIST_FROM_ORIGIN = 4.5f;
constexpr std::size_t MAX_BEST_ELEMENTS_COUNT = 150u;

constexpr float ADDITIONAL_PART_OVERLAP_RELATIVE_SIZE = 0.5f;
constexpr float ADDITIONAL_SEGMENT_ANGLE = deg_2_rad(45.f);
constexpr float MAX_ANGLE_DEV_IN_SEGMENT = deg_2_rad(35.f);
constexpr float HOW_SEGMENT_IS_BETTER_COEF = 2.f;

constexpr float NEAR_POINT_DIST_IN_POINT_SELECTION = 0.3f;

const float MAX_ANGLE_DEVIATION_TAN = std::tan(MAX_ANGLE_DEVIATION);
const float MIN_ANGLE_DEV_IN_SEGMENT_COS = std::cos(MAX_ANGLE_DEV_IN_SEGMENT);

struct LeftRightSegment {
  enum Id : std::size_t {
    Left = 0u, Right, Size
  };
};

struct UpDownSegment {
  enum Id : std::size_t {
    Up = 0u, Down, Size
  };
};

struct Line {
  Line() = default;
  Line(const Eigen::Vector3f& line)
    : k{-line(0) / line(1)}
    , b{-line(2) / line(1)}
  {}

  float get_y(float x) const { return k * x + b; }

  float k{};
  float b{};
};

using PointCloudSegments = std::array<Eigen::MatrixXf, LeftRightSegment::Size>;
using WeightedLinesSet = std::array<std::list<std::pair<float, Eigen::Vector3f>>, UpDownSegment::Size>;
using WeightedLinesPair = std::array<std::pair<float, Eigen::Vector3f>, UpDownSegment::Size>;

template <typename T>
class PriorityList final {
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
      push(other_element);
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

PointCloudSegments get_main_left_right_segments(const Eigen::MatrixXf& matrix,
                                                const Eigen::Vector2f& left_down,
                                                const Eigen::Vector2f& right_up)
{
  PointCloudSegments result;
  // segment_id == 0 is left part, segment_id == 1 is right part
  for (std::size_t segment_id{0u}; segment_id < LeftRightSegment::Size; ++segment_id) {
    std::vector<int> indices;
    for (int i{0}; i < matrix.rows(); ++i) {
      const auto& point = matrix.row(i);
      if ((point.transpose().array() > left_down.array()).all() && (point.transpose().array() < right_up.array()).all()) {
        if (segment_id == LeftRightSegment::Left && point.x() < 0.f) {
          indices.push_back(i);
        } else if (segment_id == LeftRightSegment::Right && point.x() >= 0.f) {
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

bool point_fits_line(const Eigen::Vector2f& point, const Eigen::Vector3f& line, float tolerance)
{
  return get_dist(point, line) < tolerance;
}

float score_line(const Eigen::Vector3f& line, const PointCloudSegments& point_cloud_segments)
{
  std::size_t fitting_points{0u};
  for (const auto& segment : point_cloud_segments) {
    for (int i{0}; i < segment.rows(); ++i) {
      fitting_points += point_fits_line(segment.row(i), line, NEAR_POINT_DIST_IN_SCORING);
    }
  }

  return static_cast<float>(fitting_points);
}

auto find_best_lines(const PointCloudSegments& point_cloud_segments)
{
  using StoredElement = std::pair<float, Eigen::Vector3f>;
  auto score_comparer = [](const StoredElement& first_val, const StoredElement& second_val) {
    return first_val.first < second_val.first;
  };
  PriorityList<StoredElement> result{MAX_BEST_ELEMENTS_COUNT, score_comparer};

  const auto& right_part = point_cloud_segments[LeftRightSegment::Right];
  const auto& left_part = point_cloud_segments[LeftRightSegment::Left];
  // Make all possible pairs between points from letf and right.
  // No need to ese RANSAC because of fiew points number.
  // Note: this process can be parallelized
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
  return result;
}

bool are_lines_distant(const Eigen::Vector3f& line1, const Eigen::Vector3f& line2)
{
  return std::fabs(-line1(2) / line1(1) - (-line2(2) / line2(1))) >= MIN_DIST_FROM_ORIGIN;
}

auto choose_best_lines_pair(const WeightedLinesSet& lines)
{
  if (lines[UpDownSegment::Up].empty()) {
    throw std::invalid_argument{"No lines found in the upper part"};
  }
  if (lines[UpDownSegment::Down].empty()) {
    throw std::invalid_argument{"No lines found in the down part"};
  }
  float max_cos{0.};
  const auto* best_up_line = &(lines[UpDownSegment::Up].front());
  const auto* best_down_line = &(lines[UpDownSegment::Down].front());
  for (const auto& scored_up_line: lines[UpDownSegment::Up]) {
    for (const auto& scored_down_line: lines[UpDownSegment::Down]) {
      const auto& up_line = scored_up_line.second;
      const auto& down_line = scored_down_line.second;
      // Filter out too close lines
      if (!are_lines_distant(up_line, down_line)) {
        continue;
      }
      // Filter out pairs where up lines are lower than down lines
      if (-up_line(2) / up_line(1) < -down_line(2) / down_line(1)) {
        continue;
      }
      Eigen::Vector2f norm_vec_up = up_line.topRows<2>().normalized();
      Eigen::Vector2f norm_vec_down = down_line.topRows<2>().normalized();
      const float lines_cos = std::fabs(norm_vec_up.dot(norm_vec_down));
      if (lines_cos > max_cos) {
        best_up_line = &(scored_up_line);
        best_down_line = &(scored_down_line);
        max_cos = lines_cos;
      }
    }
  }
  return std::array{*best_up_line, *best_down_line};
}

Line get_segment_line(const Line& main_line, float x_start, float line_angle)
{
  float y_start = main_line.get_y(x_start);
  Line segment_line;
  segment_line.k = std::tan(line_angle);
  segment_line.b = y_start - segment_line.k * x_start;
  return segment_line;
}

template <typename Comp>
PointCloudSegments get_additional_left_right_segments(const Eigen::MatrixXf& matrix,
                                                      const Line& main_line,
                                                      const Line& segment_line,
                                                      float x_start, float x_middle, Comp comparer)
{
  PointCloudSegments result;
  // segment_id == 0 is left part, segment_id == 1 is right part
  for (std::size_t segment_id{0u}; segment_id < result.size(); ++segment_id) {
    std::vector<int> indices;
    for (int i{0}; i < matrix.rows(); ++i) {
      const auto& point = matrix.row(i);
      if (comparer(point.x(), x_start)) { // point.x() < x_start for left part and point.x() > x_start for right part
        const float y_main_line = main_line.get_y(point.x());
        const float y_segment_line = segment_line.get_y(point.x());
        const auto min_max_y = std::minmax(y_main_line, y_segment_line);
        if (point.y() > min_max_y.first && point.y() <= min_max_y.second) {
          if (segment_id == LeftRightSegment::Left && point.x() < x_middle) {
            indices.push_back(i);
          } else if (segment_id == LeftRightSegment::Right && point.x() >= x_middle) {
            indices.push_back(i);
          }
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

float get_intersection_x(const Eigen::Vector3f& line1, const Eigen::Vector3f& line2)
{
  const float num = line2(2) / line2(1) - line1(2) / line1(1);
  const float den = line1(0) / line1(1) - line2(0) / line2(1);
  return num / den;
}

bool angle_is_near(const Eigen::Vector3f& line1, const Eigen::Vector3f& line2)
{
  Eigen::Vector2f norm_vect1 = line1.topRows<2>().normalized();
  Eigen::Vector2f norm_vect2 = line2.topRows<2>().normalized();
  return std::fabs(norm_vect1.dot(norm_vect2)) > MIN_ANGLE_DEV_IN_SEGMENT_COS;
}

template <typename Comp>
auto find_best_lines_in_segment(const PointCloudSegments& point_cloud_segments,
                                const Eigen::Vector3f& main_line, float x_middle, Comp comparer)
{
  using StoredElement = std::pair<float, Eigen::Vector3f>;
  auto score_comparer = [](const StoredElement& first_val, const StoredElement& second_val) {
    return first_val.first < second_val.first;
  };
  PriorityList<StoredElement> result{MAX_BEST_ELEMENTS_COUNT, score_comparer};

  const auto& right_part = point_cloud_segments[LeftRightSegment::Right];
  const auto& left_part = point_cloud_segments[LeftRightSegment::Left];
  // Make all possible pairs between points from letf and right.
  // No need to ese RANSAC because of fiew points number.
  for (int right_point_id{0}; right_point_id < right_part.rows(); ++right_point_id) {
    const auto& right_point = right_part.row(right_point_id);
    for (int left_point_id{0}; left_point_id < left_part.rows(); ++left_point_id) {
      const auto& left_point = left_part.row(left_point_id);
      const auto line = get_line(left_point, right_point);
      const float intersection_with_main_x = get_intersection_x(main_line, line);
      // x_middle < intersection_with_main_x for left part and x_middle > intersection_with_main_x for right part
      if (comparer(x_middle, intersection_with_main_x)
          && !comparer(-x_middle, intersection_with_main_x)
          && angle_is_near(main_line, line)) {
        const float line_score = score_line(line, point_cloud_segments);
        result.emplace(line_score, line);
      }
    }
  }
  return result;
}

float get_pair_score(const WeightedLinesPair& lines)
{
  return lines[UpDownSegment::Up].first + lines[UpDownSegment::Down].first;
}

WeightedLinesPair get_best_segment_pair(const WeightedLinesPair& best_main_pair,
                                        const std::array<float, UpDownSegment::Size>& main_pair_scores,
                                        const WeightedLinesSet& segment_results)
{
  WeightedLinesPair best_segment_pair;
  try {
    best_segment_pair = choose_best_lines_pair(segment_results);
    const auto min_score = HOW_SEGMENT_IS_BETTER_COEF * (main_pair_scores[UpDownSegment::Down] + main_pair_scores[UpDownSegment::Up]);
    if (get_pair_score(best_segment_pair) >= min_score) {
      return best_segment_pair;
    }
  } catch (const std::invalid_argument&) {
  }
  best_segment_pair = best_main_pair;
  best_segment_pair[UpDownSegment::Down].first = main_pair_scores[UpDownSegment::Down];
  best_segment_pair[UpDownSegment::Up].first = main_pair_scores[UpDownSegment::Up];
  return best_segment_pair;
}

template <typename Comp>
Eigen::MatrixXf fuse_segments(const PointCloudSegments& left_right_segments,
                              float x_start, Comp comparer)
{
  std::array<std::vector<int>, 2> indices;
  for (std::size_t segment_id{0u}; segment_id < LeftRightSegment::Size; ++segment_id) {
    for (int id{0}; id < left_right_segments[segment_id].rows(); ++id)
    {
      const float point_x = left_right_segments[segment_id].row(id).x();
      if (comparer(point_x, x_start)) { // point_x < x_start for left part and point_x > x_start for right part
        indices[segment_id].push_back(id);
      }
    }
  }
  const std::size_t total_indices = indices[LeftRightSegment::Left].size() + indices[LeftRightSegment::Right].size();
  Eigen::MatrixXf result{total_indices, 2};
  int result_id{0};
  for (std::size_t segment_id{0u}; segment_id < LeftRightSegment::Size; ++segment_id) {
    for (int id : indices[segment_id])
    {
      result.row(result_id) = left_right_segments[segment_id].row(id);
      ++result_id;
    }
  }
  return result;
}

Eigen::MatrixXf get_near_points_in_segment(const Eigen::Vector3f& line, const Eigen::MatrixXf& points)
{
  std::vector<int> indices;
  for (int id{0}; id < points.rows(); ++id) {
    if (point_fits_line(points.row(id), line, NEAR_POINT_DIST_IN_POINT_SELECTION)) {
      indices.push_back(id);
    }
  }
  Eigen::MatrixXf result_points{indices.size(), 2};
  int result_id{0};
  for (int id : indices) {
    result_points.row(result_id) = points.row(id);
    ++result_id;
  }
  return result_points;
}

Eigen::MatrixXf concatenate_points(const std::initializer_list<const Eigen::MatrixXf*>& point_matrices)
{
  Eigen::MatrixXf result;
  if (point_matrices.size() == 0u) {
    return result;
  }
  int rows{0};
  for (const auto* point_matrix : point_matrices) {
    rows += point_matrix->rows();
  }
  const auto cols = (*point_matrices.begin())->cols();
  result.conservativeResize(rows, cols);
  int start_row{0};
  for (const auto* point_matrix : point_matrices) {
    result.block(start_row, 0, point_matrix->rows(), cols) = *point_matrix;
    start_row += point_matrix->rows();
  }
  return result;
}

Eigen::Vector4f approxymate_points_by_polinomials(const Eigen::MatrixXf& points)
{
  Eigen::MatrixXf A{points.rows(), 4};
  for (int i{0}; i < points.rows(); ++i) {
    const float x = points.row(i).x();
    const float x_square = square(x);
    A(i, 0) = x * x_square; // x^3
    A(i, 1) = x_square;     // x^2
    A(i, 2) = x;            // x
    A(i, 3) = 1.f;          // constant
  }

  // Solve the normal equation A^TAx = A^Ty for x using the least squares method
  Eigen::Vector4f coefficients = (A.transpose() * A).ldlt().solve(A.transpose() * points.col(1));
  return coefficients;
}


}

namespace processing_logic {

PolynomialsVector find_lines(const PointsVector& cloud, const Vec2D& main_direction)
{
  Eigen::MatrixXf matrix_cloud = point_cloud_to_matrix(cloud);
  // Rotate the point cloud to put it along the main adirection
  const float main_direction_angle = std::atan2(main_direction.y, main_direction.x);
  rotate_matrix(matrix_cloud, -main_direction_angle);
  // Find size of the cloud
  Eigen::Vector2f cloud_down_left;
  Eigen::Vector2f cloud_up_right;
  get_min_max(matrix_cloud, cloud_down_left, cloud_up_right);
  const float cloud_length = cloud_up_right.x() - cloud_down_left.x();
  // Find length of the main part
  const float main_part_lenght = cloud_length * MAIN_PART_RELATIVE_SIZE;
  // Find best lines in the main part
  // segment_id == 0 is upper part, segment_id == 1 is down part
  WeightedLinesSet main_part_results;
  std::array<Eigen::MatrixXf, UpDownSegment::Size> main_segment_points;
  for (std::size_t segment_id{0u}; segment_id < UpDownSegment::Size; ++segment_id) {
    Eigen::Vector2f left_down;
    Eigen::Vector2f right_up;
    if (segment_id == UpDownSegment::Up) { // upper part
      left_down = Eigen::Vector2f{-main_part_lenght / 2.f, 0.f};
      right_up = Eigen::Vector2f{main_part_lenght / 2.f, MAIN_PART_SEGMENT_WIDTH};
    } else { // down part
      left_down = Eigen::Vector2f{-main_part_lenght / 2.f, -MAIN_PART_SEGMENT_WIDTH};
      right_up = Eigen::Vector2f{main_part_lenght / 2.f, 0.f};
    }
    auto left_right_segments = get_main_left_right_segments(matrix_cloud, left_down, right_up);
    main_part_results[segment_id] = find_best_lines(left_right_segments).obtain_list();
    main_segment_points[segment_id] = concatenate_points({&(left_right_segments[LeftRightSegment::Left]),
                                                          &(left_right_segments[LeftRightSegment::Right])});
  }
  const auto best_main_pair = choose_best_lines_pair(main_part_results);
  std::array<float, UpDownSegment::Size> main_pair_scores{{0.f, 0.f}};

  // Find best lines in left additional segments
  WeightedLinesSet left_segment_results;
  std::array<Eigen::MatrixXf, UpDownSegment::Size> left_segment_points;
  for (std::size_t segment_id{0u}; segment_id < UpDownSegment::Size; ++segment_id) {
    const auto& main_line_coefs = best_main_pair[segment_id].second;
    const Line main_line{main_line_coefs};
    const float segment_x_start = -cloud_length / 2.f * (MAIN_PART_RELATIVE_SIZE - ADDITIONAL_PART_OVERLAP_RELATIVE_SIZE);
    const float segment_x_middle = (segment_x_start + cloud_down_left.x()) / 2.f;
    const auto up_segment_line = get_segment_line(main_line, segment_x_start, -ADDITIONAL_SEGMENT_ANGLE);
    auto up_segments = get_additional_left_right_segments(matrix_cloud,
                                                          main_line,
                                                          up_segment_line,
                                                          segment_x_start,
                                                          segment_x_middle,
                                                          std::less<float>{});
    auto best_lines = find_best_lines_in_segment(up_segments, main_line_coefs, segment_x_middle, std::less<float>{});
    const auto down_segment_line = get_segment_line(main_line, segment_x_start, ADDITIONAL_SEGMENT_ANGLE);
    auto down_segments = get_additional_left_right_segments(matrix_cloud,
                                                            main_line,
                                                            down_segment_line,
                                                            segment_x_start,
                                                            segment_x_middle,
                                                            std::less<float>{});
    best_lines.merge(find_best_lines_in_segment(down_segments, main_line_coefs, segment_x_middle, std::less<float>{}));
    main_pair_scores[segment_id] = score_line(main_line_coefs, up_segments) + score_line(main_line_coefs, down_segments);
    left_segment_results[segment_id] = best_lines.obtain_list();

    const float points_segment_x_start = -cloud_length / 2.f * MAIN_PART_RELATIVE_SIZE;
    Eigen::MatrixXf up_segment_points = fuse_segments(up_segments, points_segment_x_start, std::less<float>{});
    Eigen::MatrixXf down_segment_points = fuse_segments(down_segments, points_segment_x_start, std::less<float>{});
    left_segment_points[segment_id] = concatenate_points({&up_segment_points, &down_segment_points});
  }
  const auto best_left_segment_pair = get_best_segment_pair(best_main_pair, main_pair_scores, left_segment_results);

  // Find best lines in right additional segments
  WeightedLinesSet right_segment_results;
  std::array<Eigen::MatrixXf, UpDownSegment::Size> right_segment_points;
  for (std::size_t segment_id{0u}; segment_id < UpDownSegment::Size; ++segment_id) {
    const auto& main_line_coefs = best_main_pair[segment_id].second;
    const Line main_line{main_line_coefs};
    const float segment_x_start = cloud_length / 2.f * (MAIN_PART_RELATIVE_SIZE - ADDITIONAL_PART_OVERLAP_RELATIVE_SIZE);
    const float segment_x_middle = (cloud_up_right.x() + segment_x_start) / 2.f;
    const auto up_segment_line = get_segment_line(main_line, segment_x_start, ADDITIONAL_SEGMENT_ANGLE);
    auto up_segments = get_additional_left_right_segments(matrix_cloud,
                                                          main_line,
                                                          up_segment_line,
                                                          segment_x_start,
                                                          segment_x_middle,
                                                          std::greater<float>{});
    auto best_lines = find_best_lines_in_segment(up_segments, main_line_coefs, segment_x_middle, std::greater<float>{});
    const auto down_segment_line = get_segment_line(main_line, segment_x_start, -ADDITIONAL_SEGMENT_ANGLE);
    auto down_segments = get_additional_left_right_segments(matrix_cloud,
                                                            main_line,
                                                            down_segment_line,
                                                            segment_x_start,
                                                            segment_x_middle,
                                                            std::greater<float>{});
    best_lines.merge(find_best_lines_in_segment(down_segments, main_line_coefs, segment_x_middle, std::greater<float>{}));
    main_pair_scores[segment_id] = score_line(main_line_coefs, up_segments) + score_line(main_line_coefs, down_segments);
    right_segment_results[segment_id] = best_lines.obtain_list();

    const float points_segment_x_start = cloud_length / 2.f * MAIN_PART_RELATIVE_SIZE;
    Eigen::MatrixXf up_segment_points = fuse_segments(up_segments, points_segment_x_start, std::greater<float>{});
    Eigen::MatrixXf down_segment_points = fuse_segments(down_segments, points_segment_x_start, std::greater<float>{});
    right_segment_points[segment_id] = concatenate_points({&up_segment_points, &down_segment_points});
  }
  const auto best_right_segment_pair = get_best_segment_pair(best_main_pair, main_pair_scores, right_segment_results);

  PolynomialsVector result{};
  for (std::size_t segment_id{0u}; segment_id < UpDownSegment::Size; ++segment_id) {
    // Find points which lay near lines in each segment
    auto left_segment_line_points = get_near_points_in_segment(best_left_segment_pair[segment_id].second, left_segment_points[segment_id]);
    auto main_segment_line_points = get_near_points_in_segment(best_main_pair[segment_id].second, main_segment_points[segment_id]);
    auto right_segment_line_points = get_near_points_in_segment(best_right_segment_pair[segment_id].second, right_segment_points[segment_id]);
    auto line_points = concatenate_points({&left_segment_line_points, &main_segment_line_points, &right_segment_line_points});
    // Rotate points to the original coordinales system
    rotate_matrix(line_points, main_direction_angle);
    // Approximate points by a polynomial
    auto polynomial_vec = approxymate_points_by_polinomials(line_points);
    Polynomial polynomial{polynomial_vec(0), polynomial_vec(1), polynomial_vec(2), polynomial_vec(3)};
    result.push_back(polynomial);
  }
  return result;
}

}
