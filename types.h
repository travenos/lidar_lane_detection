#pragma once

#include <map>
#include <vector>

constexpr inline float square(float val) { return val * val; }

struct PlainPointXYZI
{
  float x;
  float y;
  float z;
  float intensity;
};

using PointsVector = std::vector<PlainPointXYZI>;
using PointsMap = std::map<int, PointsVector>;

struct Polynomial
{
  float coef0;
  float coef1;
  float coef2;
  float coef3;
};

using PolynomialsVector = std::vector<Polynomial>;

struct Vec2D {
  float x;
  float y;
};
