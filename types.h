#pragma once

#include <map>
#include <vector>

struct PlainPointXYZI
{
  float x{};
  float y{};
  float z{};
  float intensity{};
};

using PointsVector = std::vector<PlainPointXYZI>;
using PointsMap = std::map<int, PointsVector>;
