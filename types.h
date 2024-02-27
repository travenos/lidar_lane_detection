#pragma once

#include <map>
#include <vector>

constexpr inline float square(float val) { return val * val; }

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
using ChanneledClusteredPointClouds = std::vector<ClusteredPointClouds>; // TODO!!! May be remove

struct Polynomial
{
  float coef0{};
  float coef1{};
  float coef2{};
  float coef3{};
};

struct WeightedPolynomial: Polynomial // TODO!!!
{
  float weight{};
};

using WeightedPolynomialsVector = std::vector<WeightedPolynomial>; // TODO!!!
using PolynomialsVector = std::vector<Polynomial>;

struct Vec2D {
  float x{};
  float y{};
};
