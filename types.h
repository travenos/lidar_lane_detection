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
using ChanneledClusteredPointClouds = std::vector<ClusteredPointClouds>;

struct Polynomial
{
  float coef0{};
  float coef1{};
  float coef2{};
  float coef3{};
};

struct WeightedPolynomial: Polynomial
{
  float weight{};
};

using WeightedPolynomialsVector = std::vector<WeightedPolynomial>;
