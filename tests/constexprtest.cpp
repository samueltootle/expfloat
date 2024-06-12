#include <cmath>
#include <stdio.h>
#include "../include/expansion_math.h"

static inline void TEST_64bit_product() {
  constexpr double a = 1.0 / 6.0;
  constexpr double b = 3.0e-6;
  constexpr double c = M_PI;
  constexpr expansion_math::float2<float> sa = expansion_math::split<float>(a);
  constexpr expansion_math::float2<float> sb = expansion_math::split<float>(b);
  constexpr expansion_math::float2<float> sc = expansion_math::split<float>(c);

  constexpr expansion_math::float2<float> prod = sa * (sb * sc);
  constexpr double res_d = a * b * c;
  constexpr double res   = expansion_math::recast_sum<double>(prod);
  static_assert(std::abs(1. - res / (res_d)) <= 1e-13);
}

static inline void TEST_64bit_sum() {
  constexpr double a = 1.0 / 6.0;
  constexpr double b = 3.0e-6;
  constexpr double c = M_PI;
  constexpr expansion_math::float2<float> sa = expansion_math::split<float>(a);
  constexpr expansion_math::float2<float> sb = expansion_math::split<float>(b);
  constexpr expansion_math::float2<float> sc = expansion_math::split<float>(c);

  constexpr expansion_math::float2<float> sum = sa + sb + sc;

  constexpr double res_d = a + b + c;
  constexpr double res   = expansion_math::recast_sum<double>(sum);
  static_assert(std::abs(1. - res / (res_d)) <= 1e-13);
}

static inline void TEST_64bit_division() {
  constexpr double b = 6.0;
  constexpr double c = M_PI;
  constexpr expansion_math::float2<float> sb = expansion_math::split<float>(b);
  constexpr expansion_math::float2<float> sc = expansion_math::split<float>(c);

  constexpr expansion_math::float2 res = sc / sb;
  constexpr double res_d = c / b;

  static_assert(std::abs(1. - expansion_math::recast_sum<double>(res) / (res_d)) <= 1e-13);
}

static inline void TEST_negative_float2() {
  static constexpr double a = 1.0 / 6.0;
  static constexpr double res_d = -1.0 / 6.0;
  static constexpr expansion_math::float2<float> sa = expansion_math::split<float>(a);
  static constexpr expansion_math::float2<float> res = -sa;

  static_assert(std::abs(1. - expansion_math::recast_sum<double>(res) / (res_d)) <= 1e-13);
}

static inline void TEST_64bit_pow() {
  static constexpr double a = M_PI;
  static constexpr double res_d = std::pow(a, 5);
  static constexpr expansion_math::float2<float> sa = expansion_math::split<float>(a);
  static constexpr expansion_math::float2<float> res = expansion_math::pow_expansion(sa, 5);
  static_assert(std::abs(1. - expansion_math::recast_sum<double>(res) / (res_d)) <= 1e-12);
}

static inline void TEST_64bit_sqrt() {
  static constexpr double a = 1.0 / M_PI;
  static constexpr double res_d = std::sqrt(a);
  static constexpr expansion_math::float2<float> sa = expansion_math::split<float>(a);
  static constexpr expansion_math::float2<float> res = expansion_math::sqrt_expansion(sa);
  static_assert(std::abs(1. - expansion_math::recast_sum<double>(res) / (res_d)) <= 1e-12);
}

static inline void TEST_64bit_exp() {
  static constexpr double a = 1.0 / M_PI;
  static constexpr double res_d = std::exp(a);
  static constexpr expansion_math::float2<float> sa = expansion_math::split<float>(a);
  static constexpr expansion_math::float2<float> res = expansion_math::exp_TAYLOR(sa);
  static_assert(std::abs(1. - expansion_math::recast_sum<double>(res) / (res_d)) <= 1e-12);
}

int main() {
  return 0;
}
