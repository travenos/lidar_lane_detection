#include "kdtree.h"
#include <array>
#include <cmath>
#include <cstddef>
#include <numeric>

using PointCoords = std::array<float, 3>;

// Structure to represent node of kd tree
struct Node
{
  PointCoords point;
  int id;
  Node* left;
  Node* right;

  Node(const PointCoords& coords, int setId)
    :	point{coords}, id{setId}, left{nullptr}, right{nullptr}
  {}

  ~Node()
  {
    delete left;
    delete right;
  }
private:
  Node(const Node&) = delete;
  Node& operator=(const Node&) = delete;
};

namespace  {

void insert(const PointCoords& point, int id, std::size_t k, Node*& node)
{
  if (node == nullptr)
  {
    node = new Node{point, id};
    return;
  }

  const std::size_t next_k = (k + 1) % point.size();
  if (point.at(k) < node->point.at(k))
  {
    ::insert(point, id, next_k, node->left);
  }
  else
  {
    ::insert(point, id, next_k, node->right);
  }
}

bool is_fitting_tolerance(const PointCoords& target, const PointCoords& other_node, float distance_tol)
{
  bool is_near{false};
  bool diff_not_zeros{false};
  PointCoords diffs{};
  for (std::size_t i{0}; i < diffs.size(); ++i)
  {
    diffs[i] = std::fabs(target[i] - other_node[i]);
    is_near |= (diffs[i] <= distance_tol);
    diff_not_zeros |= (diffs[i] > 0.f);
  }
  if (!diff_not_zeros) {  // Filter duplicated points
    return false;
  }
  if (is_near)
  {
    const float square_sum = std::accumulate(diffs.begin(), diffs.end(), 0.f,
                                             [](float a, float b) { return a + b * b; });
    const float square_tol = distance_tol * distance_tol;
    is_near = (square_sum <= square_tol);
  }
  return is_near;
}

void search(const PointCoords& target, float distance_tol, std::size_t k, const Node* node, std::vector<int>& ids)
{
  if (node == nullptr)
  {
    return;
  }

  if (::is_fitting_tolerance(target, node->point, distance_tol))
  {
    ids.push_back(node->id);
  }

  const std::size_t nextK = (k + 1) % target.size();
  if (target[k] - distance_tol < node->point[k])
  {
    ::search(target, distance_tol, nextK, node->left, ids);
  }
  if (target[k] + distance_tol > node->point[k])
  {
    ::search(target, distance_tol, nextK, node->right, ids);
  }
}

}

KdTree::~KdTree()
{
  delete root_;
}

void KdTree::insert(const PlainPointXYZI& point, int id)
{
  PointCoords point_coords{{point.x, point.y, point.z}};
  ::insert(point_coords, id, 0, root_);
}

// return a list of point ids in the tree that are within distance of target
std::vector<int> KdTree::search(const PlainPointXYZI& target, float distance_tol) const
{
  std::vector<int> ids;
  PointCoords target_coords{{target.x, target.y, target.z}};
  ::search(target_coords, distance_tol, 0, root_, ids);
  return ids;
}



