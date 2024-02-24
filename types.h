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
using ClusteredPointClouds = std::vector<PointsVector>;

struct Polynom
{
  float coef0{};
  float coef1{};
  float coef2{};
  float coef3{};
};

struct WeightedPolynom: Polynom
{
  float weight{};
};

using WeightedPolynomsVector = std::vector<WeightedPolynom>;
