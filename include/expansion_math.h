//
// Created by hari on 8/14/15.
// Rewritten by Samuel Tootle: June, 2024
//

#ifndef EXPFLOAT_EXPANSION_MATH_H
#define EXPFLOAT_EXPANSION_MATH_H

#if __NVCC__
    #define FUNCTION_DECORATORS_HD __host__ __device__
    #define MY_ALIGN(n) __align__(n)
#else
    #define FUNCTION_DECORATORS_HD
    #define MY_ALIGN(n) alignas(n)
#endif

namespace expansion_math {
// #define S 16

// @brief computes the sum of two single precision floating point values and computes
// their double precision sum in the expansion form
template<class storage_prec_t>
struct MY_ALIGN(sizeof(storage_prec_t) * 2U) float2;

template<class storage_prec_t, class input_t>
FUNCTION_DECORATORS_HD
constexpr float2<storage_prec_t> split(input_t const a);

FUNCTION_DECORATORS_HD
constexpr inline float2<float> fast_two_sum(float const a, float const b);

FUNCTION_DECORATORS_HD
constexpr inline float2<float> two_sum(float const a, float const b);

FUNCTION_DECORATORS_HD
constexpr inline float2<float> grow_expansion (float const val_, float const rem_, float b);

FUNCTION_DECORATORS_HD
constexpr inline float2<float> grow_expansion (float2<float> const in, float b);

// Implementation of df64_add from ref above
FUNCTION_DECORATORS_HD
constexpr inline float2<float> grow_expansion (float2<float> const in, float2<float> const b);

FUNCTION_DECORATORS_HD
constexpr inline float2<float> fast_expansion_sum (float const val_, float const rem_, float f1, float f2);

FUNCTION_DECORATORS_HD
constexpr inline float2<float> two_product (float const a, float const b);

// scale (e1,e2) by a
FUNCTION_DECORATORS_HD
constexpr inline float2<float> scale_expansion(float const val_, float const rem_, float const a);

// Implementation of df64_prod from ref above
FUNCTION_DECORATORS_HD
constexpr inline float2<float> scale_expansion(float2<float> in, float2<float> const a);

// Implementation of df64_div from ref above
FUNCTION_DECORATORS_HD
constexpr inline float2<float> division_expansion(float2<float> numerator, float2<float> const denominator);

// Implementation of pow for expansion
FUNCTION_DECORATORS_HD
constexpr inline float2<float> pow_expansion(float2<float> var, const int exp);

// scale (e1,e2) by a
FUNCTION_DECORATORS_HD
constexpr inline float2<float> scale_expansion(float2<float> in, float const a);

FUNCTION_DECORATORS_HD
constexpr inline float2<float> daxpy (float const val_, float const rem_, float scale_factor, float grow_factor);

// Take a given float2, cast the stored value and remainder to cast_t
// before adding them together.
template<class cast_t, class float2_t>
FUNCTION_DECORATORS_HD
constexpr cast_t float2_sum(float2_t& f);

// For an arbitrary input list of float2 arguments,
// recast and sum their components using float2_sum
// and sum them all together.
template<class cast_t, class... floats_t>
FUNCTION_DECORATORS_HD
constexpr cast_t recast_sum(floats_t... input);

}
#include "expansion_math_imp.h"
#endif //EXPFLOAT_EXPANSION_MATH_H
