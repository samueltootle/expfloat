#include <cmath>
#include "../include/expansion_math.h"

constexpr static inline void TEST_64bit_product() {
  constexpr double a = 1.0 / 6.0;
  constexpr double b = 3.0e-6;
  constexpr double c = M_PI;
  constexpr expansion_math::float2<float> sa = expansion_math::split<float>(a);
  constexpr expansion_math::float2<float> sb = expansion_math::split<float>(b);
  constexpr expansion_math::float2<float> sc = expansion_math::split<float>(c);

  constexpr expansion_math::float2<float> prod = expansion_math::scale_expansion(
    sa, expansion_math::scale_expansion(sb, sc)
  );
  constexpr double res_d = a * b * c;
  constexpr double res   = expansion_math::recast_sum<double>(prod);
  static_assert(std::abs(1. - res / (res_d)) <= 1e-14);
}

static inline void TEST_64bit_sum() {
  constexpr double a = 1.0 / 6.0;
  constexpr double b = 3.0e-6;
  constexpr double c = M_PI;
  constexpr expansion_math::float2<float> sa = expansion_math::split<float>(a);
  constexpr expansion_math::float2<float> sb = expansion_math::split<float>(b);
  constexpr expansion_math::float2<float> sc = expansion_math::split<float>(c);

  constexpr expansion_math::float2<float> sum = expansion_math::grow_expansion(
    sa, expansion_math::grow_expansion(sb, sc)
  );

  constexpr double res_d = a + b + c;
  constexpr double res   = expansion_math::recast_sum<double>(sum);
  static_assert(std::abs(1. - res / (res_d)) <= 1e-14);
}

static inline void TEST_64bit_division() {
  constexpr double b = 6.0;
  constexpr double c = M_PI;
  constexpr expansion_math::float2<float> sb = expansion_math::split<float>(b);
  constexpr expansion_math::float2<float> sc = expansion_math::split<float>(c);

  constexpr expansion_math::float2 res = division_expansion(sc, sb);
  constexpr double res_d = c / b;

  static_assert(std::abs(1. - expansion_math::recast_sum<double>(res) / (res_d)) <= 1e-14);
}

static inline void TEST_negative_float2() {
  static constexpr double a = 1.0 / 6.0;
  static constexpr double res_d = -1.0 / 6.0;
  static constexpr expansion_math::float2<float> sa = expansion_math::split<float>(a);
  static constexpr expansion_math::float2<float> res = -sa;

  static_assert(std::abs(1. - expansion_math::recast_sum<double>(res) / (res_d)) <= 1e-14);
 
}

int main() {
  TEST_64bit_product();
  TEST_64bit_sum();
  TEST_64bit_division();
  TEST_negative_float2();

  return 0;
}
