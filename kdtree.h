#pragma once

#include <vector>

struct Node;

struct KdTree
{
  KdTree() = default;
  ~KdTree();
  void insert(std::vector<float> point, int id);
  // return a list of point ids in the tree that are within distance of target
  std::vector<int> search(const std::vector<float>& target, float distance_tol) const;

private:
  Node* root_{nullptr};

  KdTree(const KdTree&) = delete;
  KdTree& operator=(const KdTree&) = delete;
};

