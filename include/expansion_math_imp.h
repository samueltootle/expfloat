//
// Created by hari on 8/14/15.
// Rewritten by Samuel Tootle: June, 2024
//
#ifndef EXPFLOAT_EXPANSION_MATH_IMP_H
#define EXPFLOAT_EXPANSION_MATH_IMP_H

#include "expansion_math.h"

namespace expansion_math {

// Paper: https://andrewthall.org/papers/df64_qf128.pdf
// In principle float2<storage_prec_t> could be replaced by 
// Float<storage_prec_t, N> where N are the number of values that need
// to be stored
template<class storage_prec_t = float>
struct MY_ALIGN(sizeof(storage_prec_t) * 2U) float2 {
  FUNCTION_DECORATORS_HD
  constexpr float2() : value(0.0), remainder(0.0) {} 
  
  FUNCTION_DECORATORS_HD
  constexpr float2(storage_prec_t const & value_, storage_prec_t const & remainder_) 
      : value(value_), remainder(remainder_) {}  
  
  FUNCTION_DECORATORS_HD
  constexpr float2<storage_prec_t>& operator=(float2<storage_prec_t> const & o) {
    if(this == &o)
      return *this;
    
    value = o.value;
    remainder = o.remainder;
    return *this;
  }

  FUNCTION_DECORATORS_HD
  constexpr float2<storage_prec_t>& operator+=( float2<storage_prec_t> const & rhs) {
    value += rhs.value;
    remainder += rhs.remainder;
    return *this;
  }

  FUNCTION_DECORATORS_HD
  constexpr float2<storage_prec_t>& operator-=( float2<storage_prec_t> const & rhs) {
    value -= rhs.value;
    remainder -= rhs.remainder;
    return *this;
  }

  FUNCTION_DECORATORS_HD
  friend constexpr float2<storage_prec_t> operator+( float2<storage_prec_t> const & lhs, float2<storage_prec_t> const & rhs) {
    return grow_expansion(lhs, rhs);
  }

  FUNCTION_DECORATORS_HD
  friend constexpr float2<storage_prec_t> operator+( float const & c, float2<storage_prec_t> const & rhs) {
    float2<storage_prec_t> lhs(c, 0);
    return lhs + rhs;
  }

  FUNCTION_DECORATORS_HD
  friend constexpr float2<storage_prec_t> operator+( float2<storage_prec_t> const & lhs, float const & c) {
    return c + lhs;
  }

  FUNCTION_DECORATORS_HD
  friend constexpr float2<storage_prec_t> operator+( double const & c, float2<storage_prec_t> const & rhs) {
    float2<storage_prec_t> lhs = split<float, double>(c);
    return lhs + rhs;
  }

  FUNCTION_DECORATORS_HD
  friend constexpr float2<storage_prec_t> operator+( float2<storage_prec_t> const & lhs, double const & c) {
    return c + lhs;
  }
  
  FUNCTION_DECORATORS_HD
  friend constexpr float2<storage_prec_t> operator-(float2<storage_prec_t> const & rhs) {
    float2<storage_prec_t> lhs(-rhs.value, -rhs.remainder);
    return lhs;
  }

  FUNCTION_DECORATORS_HD
  friend constexpr float2<storage_prec_t> operator-( float2<storage_prec_t> const & lhs, float2<storage_prec_t> const & rhs) {
    return grow_expansion(lhs, -rhs);
  }

  FUNCTION_DECORATORS_HD
  friend constexpr float2<storage_prec_t> operator/( float2<storage_prec_t> const & num, float2<storage_prec_t> const & den) {
    return division_expansion(num, den);
  }
  
  FUNCTION_DECORATORS_HD
  friend constexpr float2<storage_prec_t> operator*( float2<storage_prec_t> const & p1, float2<storage_prec_t> const & p2) {
    return scale_expansion(p1, p2);
  }

  FUNCTION_DECORATORS_HD
  friend constexpr float2<storage_prec_t> operator*( float2<storage_prec_t> p1, float const & p2_) {
    float2<storage_prec_t> p2(p2_, 0.0);
    return p1 * p2;
  }

  FUNCTION_DECORATORS_HD
  friend constexpr float2<storage_prec_t> operator*( float const & p1_, float2<storage_prec_t> p2) {
    return p1_ * p2;
  }
  
  // Declarations to avoid implicit conversions
  template<class T>
  FUNCTION_DECORATORS_HD
  friend constexpr float2<storage_prec_t> operator*( float2<storage_prec_t> p1, T const & p2_) = delete;

  template<class T>
  FUNCTION_DECORATORS_HD
  friend constexpr float2<storage_prec_t> operator*(T const & p1_, float2<storage_prec_t> p2) = delete;

  template<class T>
  FUNCTION_DECORATORS_HD
  friend constexpr float2<storage_prec_t> operator+( float2<storage_prec_t> p1, T const & p2_) = delete;

  template<class T>
  FUNCTION_DECORATORS_HD
  friend constexpr float2<storage_prec_t> operator+(T const & p1_, float2<storage_prec_t> p2) = delete;
  // END deletion of implicit conversions

  // Stored values
  storage_prec_t value;
  storage_prec_t remainder;

};

template<class storage_prec_t = float, class input_t = float>
FUNCTION_DECORATORS_HD
constexpr float2<storage_prec_t> split(input_t const a) {
  // FIXME verify UL is ok
  constexpr unsigned int S = sizeof(storage_prec_t) * 4;
  float2<storage_prec_t> res;
  input_t c = ((1 << S) + 1) * a;
  input_t ab = c - a;
  res.value = storage_prec_t(c - ab);
  res.remainder = storage_prec_t(a - res.value);
  return res;
}

FUNCTION_DECORATORS_HD
constexpr inline float2<float> fast_two_sum(float const a, float const b) {
  float v = 0.0;
  float2<float> result;

  result.value = a+b;
  if ((a>b) == (a>-b)) {
      v  = result.value - a;
      result.remainder = b - v;
  } else {
      v  = result.value - b;
      result.remainder = a - v;
  }
  return result;
}

FUNCTION_DECORATORS_HD
constexpr inline float2<float> two_sum(float const a, float const b) {
  float2<float> result;

  result.value = a + b;

  float bv = result.value - a;
  float av = result.value - bv;

  float br = b - bv;
  float ar = a - av;

  result.remainder = ar + br;

  return result;
}

FUNCTION_DECORATORS_HD
constexpr inline float2<float> grow_expansion (float const val_, float const rem_, float b) {
  float2<float> tmp_res = two_sum(b, rem_);
  return two_sum(tmp_res.value, val_);
}

FUNCTION_DECORATORS_HD
constexpr inline float2<float> grow_expansion (float2<float> const in, float b) {
  return grow_expansion(in.value, in.remainder, b);
}

// Implementation of df64_add from ref above
FUNCTION_DECORATORS_HD
constexpr inline float2<float> grow_expansion (float2<float> const in, float2<float> const b) {
  float2<float> s = two_sum(in.value, b.value);
  float2<float> t = two_sum(in.remainder, b.remainder);
  s.remainder += t.value;
  s = fast_two_sum(s.value, s.remainder);
  s.remainder += t.remainder;
  s = fast_two_sum(s.value, s.remainder);
  return s;
}

FUNCTION_DECORATORS_HD
constexpr inline float2<float> fast_expansion_sum (float const val_, float const rem_, float f1, float f2) {
  float2<float> res;
  float2<float> tmp_res = two_sum(f2, rem_);
  if (val_ < f1) { // want branch to fail
      tmp_res = two_sum(tmp_res.value, val_);
      res     = two_sum(tmp_res.value, f1);
  } else {
      tmp_res = two_sum(tmp_res.value, f1);
      res     = two_sum(tmp_res.value, val_);
  }
  return res;
}

FUNCTION_DECORATORS_HD
constexpr inline float2<float> two_product (float const a, float const b) {

  float2<float> res;
  res.value = a * b;

  float2<float> sa = split<float, float>(a);
  float2<float> sb = split<float, float>(b);

  float err1 = res.value - (sa.value * sb.value);
  float err2 = err1 - (sa.remainder * sb.value);
  err1 = err2 - (sa.value * sb.remainder);

  res.remainder = (sa.remainder * sb.remainder) - err1;
  return res;
}

// scale (e1,e2) by a
FUNCTION_DECORATORS_HD
constexpr inline float2<float> scale_expansion(float const val_, float const rem_, float const a) {
  float2<float> tmp_res1 = two_product(rem_, a);
  float2<float> tmp_res2 = two_product(val_, a);
  
  tmp_res1 = two_sum(tmp_res1.value, tmp_res2.remainder);
  return two_sum(tmp_res1.value, tmp_res2.value);
}

// Implementation of df64_prod from ref above
FUNCTION_DECORATORS_HD
constexpr inline float2<float> scale_expansion(float2<float> in, float2<float> const a) {
  float2<float> p = two_product(in.value, a.value);
  p.remainder += in.value * a.remainder;
  p.remainder += in.remainder * a.value;
  p = fast_two_sum(p.value, p.remainder);
  return p;
}

// Implementation of df64_div from ref above
FUNCTION_DECORATORS_HD
constexpr inline float2<float> division_expansion(float2<float> numerator, float2<float> const denominator) {
  // FIXME need to throw if denominator is zero
  float xn = 1.0f / denominator.value;
  float2<float> yn(numerator.value * xn, 0.0f);
  float diff = 0.0;
  {
    float2<float> tmp = scale_expansion(denominator, yn);
    diff = grow_expansion(numerator, float2<float>(-tmp.value, -tmp.remainder)).value;
  }

  float2<float> p = two_product(xn, diff);
  return grow_expansion(yn, p);
}

// Implementation of pow for expansion
FUNCTION_DECORATORS_HD
constexpr inline float2<float> pow_expansion(float2<float> var, const int exp) {
  float2<float> res;
  if(exp == 0) {
    res = float2<float>(1.0, 0);
  } else if(exp == 1) {
    res = var;
  } else {
    res = var;
    for(int p = 2; p <= exp; ++p) {
      res = scale_expansion(res, var);
    }
  } 
  return res;  
}

// scale (e1,e2) by a
FUNCTION_DECORATORS_HD
constexpr inline float2<float> scale_expansion(float2<float> in, float const a) {
  return scale_expansion(in.value, in.remainder, a);
}

FUNCTION_DECORATORS_HD
constexpr inline float2<float> daxpy (float const val_, float const rem_, float scale_factor, float grow_factor) {
  float2<float> tmp_res = scale_expansion(val_, rem_, scale_factor);
  return grow_expansion(tmp_res.value, tmp_res.remainder, grow_factor);
}

FUNCTION_DECORATORS_HD
constexpr inline float2<float> sqrt_expansion(float2<float> x) {
  float xn = 1.0 / std::sqrt(x.value);
  float2<float> yn(x.value * xn, 0.0);
  float2<float> ynsq = pow_expansion(yn, 2);

  float diff = (x - ynsq).value;
  float2<float> prod = two_product(xn, diff) / float2<float>(2.0, 0.0);

  return yn + prod;
}

// computes e^x using taylor expansion form - see ref, Fig 3
FUNCTION_DECORATORS_HD
constexpr inline float2<float> exp_TAYLOR(float2<float> x) {
  const float thresh = 1.0e-20 * std::exp(x.value);
  float2<float> s = float2<float>(1.0f, 0) + x;
  float2<float> p = x * x;
  float m = 2.0f;
  float2<float> f = float2<float>(2.0f, 0);
  float2<float> t = p / f;

  while(std::fabs(t.value) > thresh) {
    s = s + t;
    p = p * x;
    m += 1.0f;
    f = f * float2<float>(m, 0.0f);
    t = p / f;
  }
  return s + t;
}

// Take a given float2, cast the stored value and remainder to cast_t
// before adding them together.
template<class cast_t, class float2_t>
FUNCTION_DECORATORS_HD
constexpr cast_t float2_sum(float2_t& f) {
  return static_cast<cast_t>(f.value) + static_cast<cast_t>(f.remainder);
}

// For an arbitrary input list of float2 arguments,
// recast and sum their components using float2_sum
// and sum them all together.
template<class cast_t, class... floats_t>
FUNCTION_DECORATORS_HD
constexpr cast_t recast_sum(floats_t... input) {
  return (float2_sum<cast_t>(input) + ...);
}

// todo Add code for dot product, matvec, dgemm, rk4

// FUNCTION_DECORATORS_HD
// inline void exp_dot (float* a, float* b, unsigned int n, float& r1, float& r2) {
//     float x=0.0, y=0.0;
//     r1=0.0; r2=0.0;
//     for (int i = 0; i < n; ++i) {
//         // two_product(a[i], b[i], x, y);
//         // fast_expansion_sum(r1, r2, x, y);
//         grow_expansion(r1, r2, a[i]*b[i]);
//     }
// }

// template <typename T>
// FUNCTION_DECORATORS_HD
// inline T dot (T* a, T* b, int n) {
//     T res = 0.0;
//     for (int i = 0; i < n; ++i) {
//         res += a[i]*b[i];
//     }
//     return res;
// }

}

#endif //EXPFLOAT_EXPANSION_MATH_H