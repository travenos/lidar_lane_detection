#include "kdtree.h"
#include <cassert>
#include <cmath>
#include <cstddef>
#include <numeric>


// Structure to represent node of kd tree
struct Node
{
  std::vector<float> point;
  int id;
  Node* left;
  Node* right;

  Node(std::vector<float> arr, int setId)
    :	point{std::move(arr)}, id{setId}, left{nullptr}, right{nullptr}
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

void insert(std::vector<float> point, int id, std::size_t k, Node*& node)
{
  if (node == nullptr)
  {
    node = new Node{std::move(point), id};
    return;
  }

  assert(point.size() == node->point.size());

  const std::size_t next_k = (k + 1) % point.size();
  if (point.at(k) < node->point.at(k))
  {
    ::insert(std::move(point), id, next_k, node->left);
  }
  else
  {
    ::insert(std::move(point), id, next_k, node->right);
  }
}

bool is_fitting_tolerance(const std::vector<float>& target, const std::vector<float>& other_node, float distance_tol)
{
  assert(target.size() == other_node.size());

  bool isNear{false};
  std::vector<float> diffs(target.size());
  for (std::size_t i{0}; i < target.size(); ++i)
  {
    diffs[i] = std::fabs(target[i] - other_node[i]);
    isNear |= (diffs[i] <= distance_tol);
  }
  if (isNear)
  {
    const float square_sum = std::accumulate(diffs.begin(), diffs.end(), 0, [](float a, float b) {return a + b * b;});
    const float square_tol = distance_tol * distance_tol;
    isNear = (square_sum <= square_tol);
  }
  return isNear;
}

void search(const std::vector<float>& target, float distance_tol, std::size_t k, const Node* node, std::vector<int>& ids)
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

void KdTree::insert(std::vector<float> point, int id)
{
  ::insert(std::move(point), id, 0, root_);
}

// return a list of point ids in the tree that are within distance of target
std::vector<int> KdTree::search(const std::vector<float>& target, float distance_tol) const
{
  std::vector<int> ids;
  ::search(target, distance_tol, 0, root_, ids);
  return ids;
}



