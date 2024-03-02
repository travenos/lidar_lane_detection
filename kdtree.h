#pragma once

#include "types.h"

#include <vector>

struct Node;

struct KdTree
{
  KdTree() = default;
  ~KdTree();
  void insert(const PlainPointXYZI& point, int id);
  // Return a list of point ids in the tree that are within distance of target
  std::vector<int> search(const PlainPointXYZI& target, float distance_tol) const;

private:
  Node* root_{nullptr};

  KdTree(const KdTree&) = delete;
  KdTree& operator=(const KdTree&) = delete;
};

