#include "main_direction_getter.h"
#include <Eigen/Eigenvalues>
#include <stdexcept>

constexpr int DIMENSIONS = 2;

namespace processing_logic {


Vec2D get_main_direction(const PointsVector& points)
{
  Eigen::MatrixXf data{points.size(), DIMENSIONS};
  for (std::size_t i{0u}; i < points.size(); ++i) {
    auto& point = points.at(i);
    data.row(static_cast<int>(i)).x() = point.x;
    data.row(static_cast<int>(i)).y() = point.y;
  }

  Eigen::RowVectorXf mean = data.colwise().mean();
  data = data.rowwise() - mean;

  Eigen::MatrixXf cov = (data.adjoint() * data / static_cast<float>(data.rows() - 1));
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eigensolver{cov};

  if (eigensolver.info() != Eigen::Success) {
    throw std::logic_error{"Failed to compute PCA: " + std::to_string(eigensolver.info())};
  }

  Eigen::VectorXf eigenvalues = eigensolver.eigenvalues();
  Eigen::MatrixXf eigenvectors = eigensolver.eigenvectors();

  Eigen::VectorXf::Index max_index;
  eigenvalues.maxCoeff(&max_index);

  auto main_direction_vec = eigenvectors.row(max_index);
  if (main_direction_vec.x() < 0) {
    main_direction_vec *= (-1);
  }
  return Vec2D{main_direction_vec.x(), main_direction_vec.y()};
}

}
