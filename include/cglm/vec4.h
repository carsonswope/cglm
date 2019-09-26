/*
 * Copyright (c), Recep Aslantas.
 *
 * MIT License (MIT), http://opensource.org/licenses/MIT
 * Full license can be found in the LICENSE file
 */

/*
 Macros:
   GLM_VEC4_ONE_INIT
   GLM_VEC4_BLACK_INIT
   GLM_VEC4_ZERO_INIT
   GLM_VEC4_ONE
   GLM_VEC4_BLACK
   GLM_VEC4_ZERO

 Functions:
   CGLM_INLINE void  glm_vec4(vec3 v3, float last, vec4 dest);
   CGLM_INLINE void  glm_vec4_copy3(vec4 a, vec3 dest);
   CGLM_INLINE void  glm_vec4_copy(vec4 v, vec4 dest);
   CGLM_INLINE void  glm_vec4_ucopy(vec4 v, vec4 dest);
   CGLM_INLINE float glm_vec4_dot(vec4 a, vec4 b);
   CGLM_INLINE float glm_vec4_norm2(vec4 v);
   CGLM_INLINE float glm_vec4_norm(vec4 v);
   CGLM_INLINE void  glm_vec4_add(vec4 a, vec4 b, vec4 dest);
   CGLM_INLINE void  glm_vec4_adds(vec4 v, float s, vec4 dest);
   CGLM_INLINE void  glm_vec4_sub(vec4 a, vec4 b, vec4 dest);
   CGLM_INLINE void  glm_vec4_subs(vec4 v, float s, vec4 dest);
   CGLM_INLINE void  glm_vec4_mul(vec4 a, vec4 b, vec4 dest);
   CGLM_INLINE void  glm_vec4_scale(vec4 v, float s, vec4 dest);
   CGLM_INLINE void  glm_vec4_scale_as(vec4 v, float s, vec4 dest);
   CGLM_INLINE void  glm_vec4_div(vec4 a, vec4 b, vec4 dest);
   CGLM_INLINE void  glm_vec4_divs(vec4 v, float s, vec4 dest);
   CGLM_INLINE void  glm_vec4_addadd(vec4 a, vec4 b, vec4 dest);
   CGLM_INLINE void  glm_vec4_subadd(vec4 a, vec4 b, vec4 dest);
   CGLM_INLINE void  glm_vec4_muladd(vec4 a, vec4 b, vec4 dest);
   CGLM_INLINE void  glm_vec4_muladds(vec4 a, float s, vec4 dest);
   CGLM_INLINE void  glm_vec4_maxadd(vec4 a, vec4 b, vec4 dest);
   CGLM_INLINE void  glm_vec4_minadd(vec4 a, vec4 b, vec4 dest);
   CGLM_INLINE void  glm_vec4_negate(vec4 v);
   CGLM_INLINE void  glm_vec4_inv(vec4 v);
   CGLM_INLINE void  glm_vec4_inv_to(vec4 v, vec4 dest);
   CGLM_INLINE void  glm_vec4_normalize(vec4 v);
   CGLM_INLINE void  glm_vec4_normalize_to(vec4 vec, vec4 dest);
   CGLM_INLINE float glm_vec4_distance(vec4 a, vec4 b);
   CGLM_INLINE void  glm_vec4_maxv(vec4 a, vec4 b, vec4 dest);
   CGLM_INLINE void  glm_vec4_minv(vec4 a, vec4 b, vec4 dest);
   CGLM_INLINE void  glm_vec4_clamp(vec4 v, float minVal, float maxVal);
   CGLM_INLINE void  glm_vec4_lerp(vec4 from, vec4 to, float t, vec4 dest)

 DEPRECATED:
   glm_vec4_dup
   glm_vec4_flipsign
   glm_vec4_flipsign_to
   glm_vec4_inv
   glm_vec4_inv_to
   glm_vec4_mulv
 */

#ifndef cglm_vec4_h
#define cglm_vec4_h

#include "common.h"
#include "vec4-ext.h"
#include "util.h"

/* DEPRECATED! functions */
#define glm_vec4_dup3(v, dest)         glm_vec4_copy3(v, dest)
#define glm_vec4_dup(v, dest)          glm_vec4_copy(v, dest)
#define glm_vec4_flipsign(v)           glm_vec4_negate(v)
#define glm_vec4_flipsign_to(v, dest)  glm_vec4_negate_to(v, dest)
#define glm_vec4_inv(v)                glm_vec4_negate(v)
#define glm_vec4_inv_to(v, dest)       glm_vec4_negate_to(v, dest)
#define glm_vec4_mulv(a, b, d)         glm_vec4_mul(a, b, d)

#define GLM_VEC4_ONE_INIT   {1.0f, 1.0f, 1.0f, 1.0f}
#define GLM_VEC4_BLACK_INIT {0.0f, 0.0f, 0.0f, 1.0f}
#define GLM_VEC4_ZERO_INIT  {0.0f, 0.0f, 0.0f, 0.0f}

#define GLM_VEC4_ONE        ((vec4)GLM_VEC4_ONE_INIT)
#define GLM_VEC4_BLACK      ((vec4)GLM_VEC4_BLACK_INIT)
#define GLM_VEC4_ZERO       ((vec4)GLM_VEC4_ZERO_INIT)

// /*!
//  * @brief init vec4 using vec3
//  *
//  * @param[in]  v3   vector3
//  * @param[in]  last last item
//  * @param[out] dest destination
//  */
// CGLM_INLINE
// void
// glm_vec4(vec3 v3, float last, vec4 dest) {
//   (*dest).x = (*v3).x;
//   (*dest).y = (*v3).y;
//   (*dest).z = (*v3).z;
//   (*dest).w = last;
// }

// /*!
//  * @brief copy first 3 members of [a] to [dest]
//  *
//  * @param[in]  a    source
//  * @param[out] dest destination
//  */
// CGLM_INLINE
// void
// glm_vec4_copy3(vec4 a, vec3 dest) {
//   (*dest).x = (*a).x;
//   (*dest).y = (*a).y;
//   (*dest).z = (*a).z;
// }

/*!
 * @brief copy all members of [a] to [dest]
 *
 * @param[in]  v    source
 * @param[out] dest destination
 */
CGLM_INLINE
void
glm_vec4_copy(vec4 v, vec4 dest) {
  (*dest).x = (*v).x;
  (*dest).y = (*v).y;
  (*dest).z = (*v).z;
  (*dest).w = (*v).w;
}

// /*!
//  * @brief copy all members of [a] to [dest]
//  *
//  * alignment is not required
//  *
//  * @param[in]  v    source
//  * @param[out] dest destination
//  */
// CGLM_INLINE
// void
// glm_vec4_ucopy(vec4 v, vec4 dest) {
//   (*dest).x = (*v).x;
//   (*dest).y = (*v).y;
//   (*dest).z = (*v).z;
//   (*dest).w = (*v).w;
// }

// /*!
//  * @brief make vector zero
//  *
//  * @param[in, out]  v vector
//  */
// CGLM_INLINE
// void
// glm_vec4_zero(vec4 v) {
// #if defined( __SSE__ ) || defined( __SSE2__ )
//   glmm_store(v, _mm_setzero_ps());
// #elif defined(CGLM_NEON_FP)
//   vst1q_f32(v, vdupq_n_f32(0.0f));
// #else
//   (*v).x = 0.0f;
//   (*v).y = 0.0f;
//   (*v).z = 0.0f;
//   (*v).w = 0.0f;
// #endif
// }

// /*!
//  * @brief make vector one
//  *
//  * @param[in, out]  v vector
//  */
// CGLM_INLINE
// void
// glm_vec4_one(vec4 v) {
// #if defined( __SSE__ ) || defined( __SSE2__ )
//   glmm_store(v, _mm_set1_ps(1.0f));
// #elif defined(CGLM_NEON_FP)
//   vst1q_f32(v, vdupq_n_f32(1.0f));
// #else
//   (*v).x = 1.0f;
//   (*v).y = 1.0f;
//   (*v).z = 1.0f;
//   (*v).w = 1.0f;
// #endif
// }

// /*!
//  * @brief vec4 dot product
//  *
//  * @param[in] a vector1
//  * @param[in] b vector2
//  *
//  * @return dot product
//  */
// CGLM_INLINE
// float
// glm_vec4_dot(vec4 a, vec4 b) {
// #if defined(CGLM_SIMD)
//   return glmm_dot(glmm_load(a), glmm_load(b));
// #else
//   return (*a).x * (*b).x + (*a).y * (*b).y + (*a).z * (*b).z + (*a).w * (*b).w;
// #endif
// }

// /*!
//  * @brief norm * norm (magnitude) of vec
//  *
//  * we can use this func instead of calling norm * norm, because it would call
//  * sqrtf fuction twice but with this func we can avoid func call, maybe this is
//  * not good name for this func
//  *
//  * @param[in] v vec4
//  *
//  * @return norm * norm
//  */
// CGLM_INLINE
// float
// glm_vec4_norm2(vec4 v) {
//   return glm_vec4_dot(v, v);
// }

// /*!
//  * @brief norm (magnitude) of vec4
//  *
//  * @param[in] v vector
//  *
//  * @return norm
//  */
// CGLM_INLINE
// float
// glm_vec4_norm(vec4 v) {
// #if defined(CGLM_SIMD)
//   return glmm_norm(glmm_load(v));
// #else
//   return sqrtf(glm_vec4_dot(v, v));
// #endif
// }

/*!
 * @brief add b vector to a vector store result in dest
 *
 * @param[in]  a    vector1
 * @param[in]  b    vector2
 * @param[out] dest destination vector
 */
CGLM_INLINE
void
glm_vec4_add(vec4 a, vec4 b, vec4 dest) {
  (*dest).x = (*a).x + (*b).x;
  (*dest).y = (*a).y + (*b).y;
  (*dest).z = (*a).z + (*b).z;
  (*dest).w = (*a).w + (*b).w;
}

// /*!
//  * @brief add scalar to v vector store result in dest (d = v + vec(s))
//  *
//  * @param[in]  v    vector
//  * @param[in]  s    scalar
//  * @param[out] dest destination vector
//  */
// CGLM_INLINE
// void
// glm_vec4_adds(vec4 v, float s, vec4 dest) {
// #if defined( __SSE__ ) || defined( __SSE2__ )
//   glmm_store(dest, _mm_add_ps(glmm_load(v), _mm_set1_ps(s)));
// #elif defined(CGLM_NEON_FP)
//   vst1q_f32(dest, vaddq_f32(vld1q_f32(v), vdupq_n_f32(s)));
// #else
//   (*dest).x = (*v).x + s;
//   (*dest).y = (*v).y + s;
//   (*dest).z = (*v).z + s;
//   (*dest).w = (*v).w + s;
// #endif
// }

// /*!
//  * @brief subtract b vector from a vector store result in dest (d = a - b)
//  *
//  * @param[in]  a    vector1
//  * @param[in]  b    vector2
//  * @param[out] dest destination vector
//  */
// CGLM_INLINE
// void
// glm_vec4_sub(vec4 a, vec4 b, vec4 dest) {
// #if defined( __SSE__ ) || defined( __SSE2__ )
//   glmm_store(dest, _mm_sub_ps(glmm_load(a), glmm_load(b)));
// #elif defined(CGLM_NEON_FP)
//   vst1q_f32(dest, vsubq_f32(vld1q_f32(a), vld1q_f32(b)));
// #else
//   (*dest).x = (*a).x - (*b).x;
//   (*dest).y = (*a).y - (*b).y;
//   (*dest).z = (*a).z - (*b).z;
//   (*dest).w = (*a).w - (*b).w;
// #endif
// }

// /*!
//  * @brief subtract scalar from v vector store result in dest (d = v - vec(s))
//  *
//  * @param[in]  v    vector
//  * @param[in]  s    scalar
//  * @param[out] dest destination vector
//  */
// CGLM_INLINE
// void
// glm_vec4_subs(vec4 v, float s, vec4 dest) {
// #if defined( __SSE__ ) || defined( __SSE2__ )
//   glmm_store(dest, _mm_sub_ps(glmm_load(v), _mm_set1_ps(s)));
// #elif defined(CGLM_NEON_FP)
//   vst1q_f32(dest, vsubq_f32(vld1q_f32(v), vdupq_n_f32(s)));
// #else
//   (*dest).x = (*v).x - s;
//   (*dest).y = (*v).y - s;
//   (*dest).z = (*v).z - s;
//   (*dest).w = (*v).w - s;
// #endif
// }

// /*!
//  * @brief multiply two vector (component-wise multiplication)
//  *
//  * @param a    vector1
//  * @param b    vector2
//  * @param dest dest = ((*a).x * (*b).x, (*a).y * (*b).y, (*a).z * (*b).z, (*a).w * (*b).w)
//  */
// CGLM_INLINE
// void
// glm_vec4_mul(vec4 a, vec4 b, vec4 dest) {
// #if defined( __SSE__ ) || defined( __SSE2__ )
//   glmm_store(dest, _mm_mul_ps(glmm_load(a), glmm_load(b)));
// #elif defined(CGLM_NEON_FP)
//   vst1q_f32(dest, vmulq_f32(vld1q_f32(a), vld1q_f32(b)));
// #else
//   (*dest).x = (*a).x * (*b).x;
//   (*dest).y = (*a).y * (*b).y;
//   (*dest).z = (*a).z * (*b).z;
//   (*dest).w = (*a).w * (*b).w;
// #endif
// }

// /*!
//  * @brief multiply/scale vec4 vector with scalar: result = v * s
//  *
//  * @param[in]  v    vector
//  * @param[in]  s    scalar
//  * @param[out] dest destination vector
//  */
// CGLM_INLINE
// void
// glm_vec4_scale(vec4 v, float s, vec4 dest) {
// #if defined( __SSE__ ) || defined( __SSE2__ )
//   glmm_store(dest, _mm_mul_ps(glmm_load(v), _mm_set1_ps(s)));
// #elif defined(CGLM_NEON_FP)
//   vst1q_f32(dest, vmulq_f32(vld1q_f32(v), vdupq_n_f32(s)));
// #else
//   (*dest).x = (*v).x * s;
//   (*dest).y = (*v).y * s;
//   (*dest).z = (*v).z * s;
//   (*dest).w = (*v).w * s;
// #endif
// }

// /*!
//  * @brief make vec4 vector scale as specified: result = unit(v) * s
//  *
//  * @param[in]  v    vector
//  * @param[in]  s    scalar
//  * @param[out] dest destination vector
//  */
// CGLM_INLINE
// void
// glm_vec4_scale_as(vec4 v, float s, vec4 dest) {
//   float norm;
//   norm = glm_vec4_norm(v);

//   if (norm == 0.0f) {
//     glm_vec4_zero(dest);
//     return;
//   }

//   glm_vec4_scale(v, s / norm, dest);
// }

// /*!
//  * @brief div vector with another component-wise division: d = a / b
//  *
//  * @param[in]  a    vector 1
//  * @param[in]  b    vector 2
//  * @param[out] dest result = ((*a).x/(*b).x, (*a).y/(*b).y, (*a).z/(*b).z, (*a).w/(*b).w)
//  */
// CGLM_INLINE
// void
// glm_vec4_div(vec4 a, vec4 b, vec4 dest) {
// #if defined( __SSE__ ) || defined( __SSE2__ )
//   glmm_store(dest, _mm_div_ps(glmm_load(a), glmm_load(b)));
// #else
//   (*dest).x = (*a).x / (*b).x;
//   (*dest).y = (*a).y / (*b).y;
//   (*dest).z = (*a).z / (*b).z;
//   (*dest).w = (*a).w / (*b).w;
// #endif
// }

// /*!
//  * @brief div vec4 vector with scalar: d = v / s
//  *
//  * @param[in]  v    vector
//  * @param[in]  s    scalar
//  * @param[out] dest destination vector
//  */
// CGLM_INLINE
// void
// glm_vec4_divs(vec4 v, float s, vec4 dest) {
// #if defined( __SSE__ ) || defined( __SSE2__ )
//   glmm_store(dest, _mm_div_ps(glmm_load(v), _mm_set1_ps(s)));
// #else
//   glm_vec4_scale(v, 1.0f / s, dest);
// #endif
// }

// /*!
//  * @brief add two vectors and add result to sum
//  *
//  * it applies += operator so dest must be initialized
//  *
//  * @param[in]  a    vector 1
//  * @param[in]  b    vector 2
//  * @param[out] dest dest += (a + b)
//  */
// CGLM_INLINE
// void
// glm_vec4_addadd(vec4 a, vec4 b, vec4 dest) {
// #if defined( __SSE__ ) || defined( __SSE2__ )
//   glmm_store(dest, _mm_add_ps(glmm_load(dest),
//                               _mm_add_ps(glmm_load(a),
//                                          glmm_load(b))));
// #elif defined(CGLM_NEON_FP)
//   vst1q_f32(dest, vaddq_f32(vld1q_f32(dest),
//                             vaddq_f32(vld1q_f32(a),
//                                       vld1q_f32(b))));
// #else
//   (*dest).x += (*a).x + (*b).x;
//   (*dest).y += (*a).y + (*b).y;
//   (*dest).z += (*a).z + (*b).z;
//   (*dest).w += (*a).w + (*b).w;
// #endif
// }

// /*!
//  * @brief sub two vectors and add result to dest
//  *
//  * it applies += operator so dest must be initialized
//  *
//  * @param[in]  a    vector 1
//  * @param[in]  b    vector 2
//  * @param[out] dest dest += (a - b)
//  */
// CGLM_INLINE
// void
// glm_vec4_subadd(vec4 a, vec4 b, vec4 dest) {
// #if defined( __SSE__ ) || defined( __SSE2__ )
//   glmm_store(dest, _mm_add_ps(glmm_load(dest),
//                               _mm_sub_ps(glmm_load(a),
//                                          glmm_load(b))));
// #elif defined(CGLM_NEON_FP)
//   vst1q_f32(dest, vaddq_f32(vld1q_f32(dest),
//                             vsubq_f32(vld1q_f32(a),
//                                       vld1q_f32(b))));
// #else
//   (*dest).x += (*a).x - (*b).x;
//   (*dest).y += (*a).y - (*b).y;
//   (*dest).z += (*a).z - (*b).z;
//   (*dest).w += (*a).w - (*b).w;
// #endif
// }

// /*!
//  * @brief mul two vectors and add result to dest
//  *
//  * it applies += operator so dest must be initialized
//  *
//  * @param[in]  a    vector 1
//  * @param[in]  b    vector 2
//  * @param[out] dest dest += (a * b)
//  */
// CGLM_INLINE
// void
// glm_vec4_muladd(vec4 a, vec4 b, vec4 dest) {
// #if defined( __SSE__ ) || defined( __SSE2__ )
//   glmm_store(dest, _mm_add_ps(glmm_load(dest),
//                               _mm_mul_ps(glmm_load(a),
//                                          glmm_load(b))));
// #elif defined(CGLM_NEON_FP)
//   vst1q_f32(dest, vaddq_f32(vld1q_f32(dest),
//                             vmulq_f32(vld1q_f32(a),
//                                       vld1q_f32(b))));
// #else
//   (*dest).x += (*a).x * (*b).x;
//   (*dest).y += (*a).y * (*b).y;
//   (*dest).z += (*a).z * (*b).z;
//   (*dest).w += (*a).w * (*b).w;
// #endif
// }

// /*!
//  * @brief mul vector with scalar and add result to sum
//  *
//  * it applies += operator so dest must be initialized
//  *
//  * @param[in]  a    vector
//  * @param[in]  s    scalar
//  * @param[out] dest dest += (a * b)
//  */
// CGLM_INLINE
// void
// glm_vec4_muladds(vec4 a, float s, vec4 dest) {
// #if defined( __SSE__ ) || defined( __SSE2__ )
//   glmm_store(dest, _mm_add_ps(glmm_load(dest),
//                               _mm_mul_ps(glmm_load(a),
//                                          _mm_set1_ps(s))));
// #elif defined(CGLM_NEON_FP)
//   vst1q_f32(dest, vaddq_f32(vld1q_f32(dest),
//                             vsubq_f32(vld1q_f32(a),
//                                       vdupq_n_f32(s))));
// #else
//   (*dest).x += (*a).x * s;
//   (*dest).y += (*a).y * s;
//   (*dest).z += (*a).z * s;
//   (*dest).w += (*a).w * s;
// #endif
// }

// /*!
//  * @brief add max of two vector to result/dest
//  *
//  * it applies += operator so dest must be initialized
//  *
//  * @param[in]  a    vector 1
//  * @param[in]  b    vector 2
//  * @param[out] dest dest += max(a, b)
//  */
// CGLM_INLINE
// void
// glm_vec4_maxadd(vec4 a, vec4 b, vec4 dest) {
// #if defined( __SSE__ ) || defined( __SSE2__ )
//   glmm_store(dest, _mm_add_ps(glmm_load(dest),
//                               _mm_max_ps(glmm_load(a),
//                                          glmm_load(b))));
// #elif defined(CGLM_NEON_FP)
//   vst1q_f32(dest, vaddq_f32(vld1q_f32(dest),
//                             vmaxq_f32(vld1q_f32(a),
//                                       vld1q_f32(b))));
// #else
//   (*dest).x += glm_max((*a).x, (*b).x);
//   (*dest).y += glm_max((*a).y, (*b).y);
//   (*dest).z += glm_max((*a).z, (*b).z);
//   (*dest).w += glm_max((*a).w, (*b).w);
// #endif
// }

// /*!
//  * @brief add min of two vector to result/dest
//  *
//  * it applies += operator so dest must be initialized
//  *
//  * @param[in]  a    vector 1
//  * @param[in]  b    vector 2
//  * @param[out] dest dest += min(a, b)
//  */
// CGLM_INLINE
// void
// glm_vec4_minadd(vec4 a, vec4 b, vec4 dest) {
// #if defined( __SSE__ ) || defined( __SSE2__ )
//   glmm_store(dest, _mm_add_ps(glmm_load(dest),
//                               _mm_min_ps(glmm_load(a),
//                                          glmm_load(b))));
// #elif defined(CGLM_NEON_FP)
//   vst1q_f32(dest, vaddq_f32(vld1q_f32(dest),
//                             vminq_f32(vld1q_f32(a),
//                                       vld1q_f32(b))));
// #else
//   (*dest).x += glm_min((*a).x, (*b).x);
//   (*dest).y += glm_min((*a).y, (*b).y);
//   (*dest).z += glm_min((*a).z, (*b).z);
//   (*dest).w += glm_min((*a).w, (*b).w);
// #endif
// }

// /*!
//  * @brief negate vector components and store result in dest
//  *
//  * @param[in]  v     vector
//  * @param[out] dest  result vector
//  */
// CGLM_INLINE
// void
// glm_vec4_negate_to(vec4 v, vec4 dest) {
// #if defined( __SSE__ ) || defined( __SSE2__ )
//   glmm_store(dest, _mm_xor_ps(glmm_load(v), _mm_set1_ps(-0.0f)));
// #elif defined(CGLM_NEON_FP)
//   vst1q_f32(dest, veorq_s32(vld1q_f32(v), vdupq_n_f32(-0.0f)));
// #else
//   (*dest).x = -(*v).x;
//   (*dest).y = -(*v).y;
//   (*dest).z = -(*v).z;
//   (*dest).w = -(*v).w;
// #endif
// }

// /*!
//  * @brief flip sign of all vec4 members
//  *
//  * @param[in, out]  v  vector
//  */
// CGLM_INLINE
// void
// glm_vec4_negate(vec4 v) {
//   glm_vec4_negate_to(v, v);
// }

// /*!
//  * @brief normalize vec4 to dest
//  *
//  * @param[in]  v    source
//  * @param[out] dest destination
//  */
// CGLM_INLINE
// void
// glm_vec4_normalize_to(vec4 v, vec4 dest) {
// #if defined( __SSE__ ) || defined( __SSE2__ )
//   __m128 xdot, x0;
//   float  dot;

//   x0   = glmm_load(v);
//   xdot = glmm_vdot(x0, x0);
//   dot  = _mm_cvtss_f32(xdot);

//   if (dot == 0.0f) {
//     glmm_store(dest, _mm_setzero_ps());
//     return;
//   }

//   glmm_store(dest, _mm_div_ps(x0, _mm_sqrt_ps(xdot)));
// #else
//   float norm;

//   norm = glm_vec4_norm(v);

//   if (norm == 0.0f) {
//     glm_vec4_zero(dest);
//     return;
//   }

//   glm_vec4_scale(v, 1.0f / norm, dest);
// #endif
// }

// /*!
//  * @brief normalize vec4 and store result in same vec
//  *
//  * @param[in, out] v vector
//  */
// CGLM_INLINE
// void
// glm_vec4_normalize(vec4 v) {
//   glm_vec4_normalize_to(v, v);
// }

// /**
//  * @brief distance between two vectors
//  *
//  * @param[in] a vector1
//  * @param[in] b vector2
//  * @return returns distance
//  */
// CGLM_INLINE
// float
// glm_vec4_distance(vec4 a, vec4 b) {
// #if defined( __SSE__ ) || defined( __SSE2__ )
//   return glmm_norm(_mm_sub_ps(glmm_load(b), glmm_load(a)));
// #elif defined(CGLM_NEON_FP)
//   return glmm_norm(vsubq_f32(glmm_load(a), glmm_load(b)));
// #else
//   return sqrtf(glm_pow2((*b).x - (*a).x)
//              + glm_pow2((*b).y - (*a).y)
//              + glm_pow2((*b).z - (*a).z)
//              + glm_pow2((*b).w - (*a).w));
// #endif
// }

// /*!
//  * @brief max values of vectors
//  *
//  * @param[in]  a    vector1
//  * @param[in]  b    vector2
//  * @param[out] dest destination
//  */
// CGLM_INLINE
// void
// glm_vec4_maxv(vec4 a, vec4 b, vec4 dest) {
// #if defined( __SSE__ ) || defined( __SSE2__ )
//   glmm_store(dest, _mm_max_ps(glmm_load(a), glmm_load(b)));
// #elif defined(CGLM_NEON_FP)
//   vst1q_f32(dest, vmaxq_f32(vld1q_f32(a), vld1q_f32(b)));
// #else
//   (*dest).x = glm_max((*a).x, (*b).x);
//   (*dest).y = glm_max((*a).y, (*b).y);
//   (*dest).z = glm_max((*a).z, (*b).z);
//   (*dest).w = glm_max((*a).w, (*b).w);
// #endif
// }

// /*!
//  * @brief min values of vectors
//  *
//  * @param[in]  a    vector1
//  * @param[in]  b    vector2
//  * @param[out] dest destination
//  */
// CGLM_INLINE
// void
// glm_vec4_minv(vec4 a, vec4 b, vec4 dest) {
// #if defined( __SSE__ ) || defined( __SSE2__ )
//   glmm_store(dest, _mm_min_ps(glmm_load(a), glmm_load(b)));
// #elif defined(CGLM_NEON_FP)
//   vst1q_f32(dest, vminq_f32(vld1q_f32(a), vld1q_f32(b)));
// #else
//   (*dest).x = glm_min((*a).x, (*b).x);
//   (*dest).y = glm_min((*a).y, (*b).y);
//   (*dest).z = glm_min((*a).z, (*b).z);
//   (*dest).w = glm_min((*a).w, (*b).w);
// #endif
// }

// /*!
//  * @brief clamp vector's individual members between min and max values
//  *
//  * @param[in, out]  v      vector
//  * @param[in]       minVal minimum value
//  * @param[in]       maxVal maximum value
//  */
// CGLM_INLINE
// void
// glm_vec4_clamp(vec4 v, float minVal, float maxVal) {
// #if defined( __SSE__ ) || defined( __SSE2__ )
//   glmm_store(v, _mm_min_ps(_mm_max_ps(glmm_load(v), _mm_set1_ps(minVal)),
//                            _mm_set1_ps(maxVal)));
// #elif defined(CGLM_NEON_FP)
//   vst1q_f32(v, vminq_f32(vmaxq_f32(vld1q_f32(v), vdupq_n_f32(minVal)),
//                          vdupq_n_f32(maxVal)));
// #else
//   (*v).x = glm_clamp((*v).x, minVal, maxVal);
//   (*v).y = glm_clamp((*v).y, minVal, maxVal);
//   (*v).z = glm_clamp((*v).z, minVal, maxVal);
//   (*v).w = glm_clamp((*v).w, minVal, maxVal);
// #endif
// }

// /*!
//  * @brief linear interpolation between two vector
//  *
//  * formula:  from + s * (to - from)
//  *
//  * @param[in]   from from value
//  * @param[in]   to   to value
//  * @param[in]   t    interpolant (amount) clamped between 0 and 1
//  * @param[out]  dest destination
//  */
// CGLM_INLINE
// void
// glm_vec4_lerp(vec4 from, vec4 to, float t, vec4 dest) {
//   float4 s, v;

//   /* from + s * (to - from) */
//   glm_vec4_broadcast(glm_clamp_zo(t), s);
//   glm_vec4_sub(to, from, &v);
//   glm_vec4_mul(&s, &v, &v);
//   glm_vec4_add(from, &v, dest);
// }

// /*!
//  * @brief helper to fill vec4 as [S^3, S^2, S, 1]
//  *
//  * @param[in]   s    parameter
//  * @param[out]  dest destination
//  */
// CGLM_INLINE
// void
// glm_vec4_cubic(float s, vec4 dest) {
//   float ss;

//   ss = s * s;

//   (*dest).x = ss * s;
//   (*dest).y = ss;
//   (*dest).z = s;
//   (*dest).w = 1.0f;
// }

#endif /* cglm_vec4_h */
