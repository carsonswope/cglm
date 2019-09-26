/*
 * Copyright (c), Recep Aslantas.
 *
 * MIT License (MIT), http://opensource.org/licenses/MIT
 * Full license can be found in the LICENSE file
 */

/*
 Macros:
   GLM_VEC3_ONE_INIT
   GLM_VEC3_ZERO_INIT
   GLM_VEC3_ONE
   GLM_VEC3_ZERO
   GLM_YUP
   GLM_ZUP
   GLM_XUP

 Functions:
   CGLM_INLINE void  glm_vec3(vec4 v4, vec3 dest);
   CGLM_INLINE void  glm_vec3_copy(vec3 a, vec3 dest);
   CGLM_INLINE void  glm_vec3_zero(vec3 v);
   CGLM_INLINE void  glm_vec3_one(vec3 v);
   CGLM_INLINE float glm_vec3_dot(vec3 a, vec3 b);
   CGLM_INLINE float glm_vec3_norm2(vec3 v);
   CGLM_INLINE float glm_vec3_norm(vec3 v);
   CGLM_INLINE void  glm_vec3_add(vec3 a, vec3 b, vec3 dest);
   CGLM_INLINE void  glm_vec3_adds(vec3 a, float s, vec3 dest);
   CGLM_INLINE void  glm_vec3_sub(vec3 a, vec3 b, vec3 dest);
   CGLM_INLINE void  glm_vec3_subs(vec3 a, float s, vec3 dest);
   CGLM_INLINE void  glm_vec3_mul(vec3 a, vec3 b, vec3 dest);
   CGLM_INLINE void  glm_vec3_scale(vec3 v, float s, vec3 dest);
   CGLM_INLINE void  glm_vec3_scale_as(vec3 v, float s, vec3 dest);
   CGLM_INLINE void  glm_vec3_div(vec3 a, vec3 b, vec3 dest);
   CGLM_INLINE void  glm_vec3_divs(vec3 a, float s, vec3 dest);
   CGLM_INLINE void  glm_vec3_addadd(vec3 a, vec3 b, vec3 dest);
   CGLM_INLINE void  glm_vec3_subadd(vec3 a, vec3 b, vec3 dest);
   CGLM_INLINE void  glm_vec3_muladd(vec3 a, vec3 b, vec3 dest);
   CGLM_INLINE void  glm_vec3_muladds(vec3 a, float s, vec3 dest);
   CGLM_INLINE void  glm_vec3_maxadd(vec3 a, vec3 b, vec3 dest);
   CGLM_INLINE void  glm_vec3_minadd(vec3 a, vec3 b, vec3 dest);
   CGLM_INLINE void  glm_vec3_flipsign(vec3 v);
   CGLM_INLINE void  glm_vec3_flipsign_to(vec3 v, vec3 dest);
   CGLM_INLINE void  glm_vec3_negate_to(vec3 v, vec3 dest);
   CGLM_INLINE void  glm_vec3_negate(vec3 v);
   CGLM_INLINE void  glm_vec3_inv(vec3 v);
   CGLM_INLINE void  glm_vec3_inv_to(vec3 v, vec3 dest);
   CGLM_INLINE void  glm_vec3_normalize(vec3 v);
   CGLM_INLINE void  glm_vec3_normalize_to(vec3 v, vec3 dest);
   CGLM_INLINE void  glm_vec3_cross(vec3 a, vec3 b, vec3 d);
   CGLM_INLINE void  glm_vec3_crossn(vec3 a, vec3 b, vec3 dest);
   CGLM_INLINE float glm_vec3_distance(vec3 a, vec3 b);
   CGLM_INLINE float glm_vec3_angle(vec3 a, vec3 b);
   CGLM_INLINE void  glm_vec3_rotate(vec3 v, float angle, vec3 axis);
   CGLM_INLINE void  glm_vec3_rotate_m4(mat4 m, vec3 v, vec3 dest);
   CGLM_INLINE void  glm_vec3_rotate_m3(mat3 m, vec3 v, vec3 dest);
   CGLM_INLINE void  glm_vec3_proj(vec3 a, vec3 b, vec3 dest);
   CGLM_INLINE void  glm_vec3_center(vec3 a, vec3 b, vec3 dest);
   CGLM_INLINE float glm_vec3_distance2(vec3 a, vec3 b);
   CGLM_INLINE void  glm_vec3_maxv(vec3 a, vec3 b, vec3 dest);
   CGLM_INLINE void  glm_vec3_minv(vec3 a, vec3 b, vec3 dest);
   CGLM_INLINE void  glm_vec3_ortho(vec3 v, vec3 dest);
   CGLM_INLINE void  glm_vec3_clamp(vec3 v, float minVal, float maxVal);
   CGLM_INLINE void  glm_vec3_lerp(vec3 from, vec3 to, float t, vec3 dest);

 Convenient:
   CGLM_INLINE void  glm_cross(vec3 a, vec3 b, vec3 d);
   CGLM_INLINE float glm_dot(vec3 a, vec3 b);
   CGLM_INLINE void  glm_normalize(vec3 v);
   CGLM_INLINE void  glm_normalize_to(vec3 v, vec3 dest);

 DEPRECATED:
   glm_vec3_dup
   glm_vec3_flipsign
   glm_vec3_flipsign_to
   glm_vec3_inv
   glm_vec3_inv_to
   glm_vec3_mulv
 */

#ifndef cglm_vec3_h
#define cglm_vec3_h

#include "common.h"
#include "vec4.h"
#include "vec3-ext.h"
#include "util.h"

/* DEPRECATED! use _copy, _ucopy versions */
#define glm_vec3_dup(v, dest)         glm_vec3_copy(v, dest)
#define glm_vec3_flipsign(v)          glm_vec3_negate(v)
#define glm_vec3_flipsign_to(v, dest) glm_vec3_negate_to(v, dest)
#define glm_vec3_inv(v)               glm_vec3_negate(v)
#define glm_vec3_inv_to(v, dest)      glm_vec3_negate_to(v, dest)
#define glm_vec3_mulv(a, b, d)        glm_vec3_mul(a, b, d)

#define GLM_VEC3_ONE_INIT   {1.0f, 1.0f, 1.0f}
#define GLM_VEC3_ZERO_INIT  {0.0f, 0.0f, 0.0f}

#define GLM_VEC3_ONE  ((vec3)GLM_VEC3_ONE_INIT)
#define GLM_VEC3_ZERO ((vec3)GLM_VEC3_ZERO_INIT)

#define GLM_YUP  ((vec3){0.0f, 1.0f, 0.0f})
#define GLM_ZUP  ((vec3){0.0f, 0.0f, 1.0f})
#define GLM_XUP  ((vec3){1.0f, 0.0f, 0.0f})

// /*!
//  * @brief init vec3 using vec4
//  *
//  * @param[in]  v4   vector4
//  * @param[out] dest destination
//  */
// CGLM_INLINE
// void
// glm_vec3(vec4 v, vec3 dest) {
//   (*dest).x = (*v).x;
//   (*dest).y = (*v).y;
//   (*dest).z = (*v).z;
// }

/*!
 * @brief copy all members of [a] to [dest]
 *
 * @param[in]  a    source
 * @param[out] dest destination
 */
CGLM_INLINE
void
glm_vec3_copy(vec3 a, vec3 dest) {
  (*dest).x = (*a).x;
  (*dest).y = (*a).y;
  (*dest).z = (*a).z;
}

// /*!
//  * @brief make vector zero
//  *
//  * @param[in, out]  v vector
//  */
// CGLM_INLINE
// void
// glm_vec3_zero(vec3 v) {
//   (*v).x = (*v).y = (*v).z = 0.0f;
// }

// /*!
//  * @brief make vector one
//  *
//  * @param[in, out]  v vector
//  */
// CGLM_INLINE
// void
// glm_vec3_one(vec3 v) {
//   (*v).x = (*v).y = (*v).z = 1.0f;
// }

/*!
 * @brief vec3 dot product
 *
 * @param[in] a vector1
 * @param[in] b vector2
 *
 * @return dot product
 */
CGLM_INLINE
float glm_vec3_dot(vec3 a, vec3 b) {
  return (*a).x * (*b).x + (*a).y * (*b).y + (*a).z * (*b).z;
}

/*!
 * @brief norm * norm (magnitude) of vec
 *
 * we can use this func instead of calling norm * norm, because it would call
 * sqrtf fuction twice but with this func we can avoid func call, maybe this is
 * not good name for this func
 *
 * @param[in] v vector
 *
 * @return norm * norm
 */
CGLM_INLINE
float
glm_vec3_norm2(vec3 v) {
  return glm_vec3_dot(v, v);
}

/*!
 * @brief norm (magnitude) of vec3
 *
 * @param[in] v vector
 *
 * @return norm
 */
CGLM_INLINE
float
glm_vec3_norm(vec3 v) {
  return sqrtf(glm_vec3_norm2(v));
}

// /*!
//  * @brief add a vector to b vector store result in dest
//  *
//  * @param[in]  a    vector1
//  * @param[in]  b    vector2
//  * @param[out] dest destination vector
//  */
// CGLM_INLINE
// void
// glm_vec3_add(vec3 a, vec3 b, vec3 dest) {
//   (*dest).x = (*a).x + (*b).x;
//   (*dest).y = (*a).y + (*b).y;
//   (*dest).z = (*a).z + (*b).z;
// }

// /*!
//  * @brief add scalar to v vector store result in dest (d = v + s)
//  *
//  * @param[in]  v    vector
//  * @param[in]  s    scalar
//  * @param[out] dest destination vector
//  */
// CGLM_INLINE
// void
// glm_vec3_adds(vec3 v, float s, vec3 dest) {
//   (*dest).x = (*v).x + s;
//   (*dest).y = (*v).y + s;
//   (*dest).z = (*v).z + s;
// }

// /*!
//  * @brief subtract b vector from a vector store result in dest
//  *
//  * @param[in]  a    vector1
//  * @param[in]  b    vector2
//  * @param[out] dest destination vector
//  */
// CGLM_INLINE
// void
// glm_vec3_sub(vec3 a, vec3 b, vec3 dest) {
//   (*dest).x = (*a).x - (*b).x;
//   (*dest).y = (*a).y - (*b).y;
//   (*dest).z = (*a).z - (*b).z;
// }

// /*!
//  * @brief subtract scalar from v vector store result in dest (d = v - s)
//  *
//  * @param[in]  v    vector
//  * @param[in]  s    scalar
//  * @param[out] dest destination vector
//  */
// CGLM_INLINE
// void
// glm_vec3_subs(vec3 v, float s, vec3 dest) {
//   (*dest).x = (*v).x - s;
//   (*dest).y = (*v).y - s;
//   (*dest).z = (*v).z - s;
// }

// /*!
//  * @brief multiply two vector (component-wise multiplication)
//  *
//  * @param a    vector1
//  * @param b    vector2
//  * @param dest v3 = ((*a).x * (*b).x, (*a).y * (*b).y, (*a).z * (*b).z)
//  */
// CGLM_INLINE
// void
// glm_vec3_mul(vec3 a, vec3 b, vec3 dest) {
//   (*dest).x = (*a).x * (*b).x;
//   (*dest).y = (*a).y * (*b).y;
//   (*dest).z = (*a).z * (*b).z;
// }

/*!
 * @brief multiply/scale vec3 vector with scalar: result = v * s
 *
 * @param[in]  v    vector
 * @param[in]  s    scalar
 * @param[out] dest destination vector
 */
CGLM_INLINE
void
glm_vec3_scale(vec3 v, float s, vec3 dest) {
  (*dest).x = (*v).x * s;
  (*dest).y = (*v).y * s;
  (*dest).z = (*v).z * s;
}

// /*!
//  * @brief make vec3 vector scale as specified: result = unit(v) * s
//  *
//  * @param[in]  v    vector
//  * @param[in]  s    scalar
//  * @param[out] dest destination vector
//  */
// CGLM_INLINE
// void
// glm_vec3_scale_as(vec3 v, float s, vec3 dest) {
//   float norm;
//   norm = glm_vec3_norm(v);

//   if (norm == 0.0f) {
//     glm_vec3_zero(dest);
//     return;
//   }

//   glm_vec3_scale(v, s / norm, dest);
// }

// /*!
//  * @brief div vector with another component-wise division: d = a / b
//  *
//  * @param[in]  a    vector 1
//  * @param[in]  b    vector 2
//  * @param[out] dest result = ((*a).x/(*b).x, (*a).y/(*b).y, (*a).z/(*b).z)
//  */
// CGLM_INLINE
// void
// glm_vec3_div(vec3 a, vec3 b, vec3 dest) {
//   (*dest).x = (*a).x / (*b).x;
//   (*dest).y = (*a).y / (*b).y;
//   (*dest).z = (*a).z / (*b).z;
// }

// /*!
//  * @brief div vector with scalar: d = v / s
//  *
//  * @param[in]  v    vector
//  * @param[in]  s    scalar
//  * @param[out] dest result = ((*a).x/s, (*a).y/s, (*a).z/s)
//  */
// CGLM_INLINE
// void
// glm_vec3_divs(vec3 v, float s, vec3 dest) {
//   (*dest).x = (*v).x / s;
//   (*dest).y = (*v).y / s;
//   (*dest).z = (*v).z / s;
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
// glm_vec3_addadd(vec3 a, vec3 b, vec3 dest) {
//   (*dest).x += (*a).x + (*b).x;
//   (*dest).y += (*a).y + (*b).y;
//   (*dest).z += (*a).z + (*b).z;
// }

// /*!
//  * @brief sub two vectors and add result to dest
//  *
//  * it applies += operator so dest must be initialized
//  *
//  * @param[in]  a    vector 1
//  * @param[in]  b    vector 2
//  * @param[out] dest dest += (a + b)
//  */
// CGLM_INLINE
// void
// glm_vec3_subadd(vec3 a, vec3 b, vec3 dest) {
//   (*dest).x += (*a).x - (*b).x;
//   (*dest).y += (*a).y - (*b).y;
//   (*dest).z += (*a).z - (*b).z;
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
// glm_vec3_muladd(vec3 a, vec3 b, vec3 dest) {
//   (*dest).x += (*a).x * (*b).x;
//   (*dest).y += (*a).y * (*b).y;
//   (*dest).z += (*a).z * (*b).z;
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
// glm_vec3_muladds(vec3 a, float s, vec3 dest) {
//   (*dest).x += (*a).x * s;
//   (*dest).y += (*a).y * s;
//   (*dest).z += (*a).z * s;
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
// glm_vec3_maxadd(vec3 a, vec3 b, vec3 dest) {
//   (*dest).x += glm_max((*a).x, (*b).x);
//   (*dest).y += glm_max((*a).y, (*b).y);
//   (*dest).z += glm_max((*a).z, (*b).z);
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
// glm_vec3_minadd(vec3 a, vec3 b, vec3 dest) {
//   (*dest).x += glm_min((*a).x, (*b).x);
//   (*dest).y += glm_min((*a).y, (*b).y);
//   (*dest).z += glm_min((*a).z, (*b).z);
// }

// /*!
//  * @brief negate vector components and store result in dest
//  *
//  * @param[in]   v     vector
//  * @param[out]  dest  result vector
//  */
// CGLM_INLINE
// void
// glm_vec3_negate_to(vec3 v, vec3 dest) {
//   (*dest).x = -(*v).x;
//   (*dest).y = -(*v).y;
//   (*dest).z = -(*v).z;
// }

// /*!
//  * @brief negate vector components
//  *
//  * @param[in, out]  v  vector
//  */
// CGLM_INLINE
// void
// glm_vec3_negate(vec3 v) {
//   glm_vec3_negate_to(v, v);
// }

/*!
 * @brief normalize vec3 and store result in same vec
 *
 * @param[in, out] v vector
 */
CGLM_INLINE
void
glm_vec3_normalize(vec3 v) {
  float norm;

  norm = glm_vec3_norm(v);

  if (norm == 0.0f) {
    (*v).x = (*v).y = (*v).z = 0.0f;
    return;
  }

  glm_vec3_scale(v, 1.0f / norm, v);
}

// /*!
//  * @brief normalize vec3 to dest
//  *
//  * @param[in]  v    source
//  * @param[out] dest destination
//  */
// CGLM_INLINE
// void
// glm_vec3_normalize_to(vec3 v, vec3 dest) {
//   float norm;

//   norm = glm_vec3_norm(v);

//   if (norm == 0.0f) {
//     glm_vec3_zero(dest);
//     return;
//   }

//   glm_vec3_scale(v, 1.0f / norm, dest);
// }

/*!
 * @brief cross product of two vector (RH)
 *
 * @param[in]  a    vector 1
 * @param[in]  b    vector 2
 * @param[out] dest destination
 */
CGLM_INLINE
void
glm_vec3_cross(vec3 a, vec3 b, vec3 dest) {
  /* (u2.v3 - u3.v2, u3.v1 - u1.v3, u1.v2 - u2.v1) */
  (*dest).x = (*a).y * (*b).z - (*a).z * (*b).y;
  (*dest).y = (*a).z * (*b).x - (*a).x * (*b).z;
  (*dest).z = (*a).x * (*b).y - (*a).y * (*b).x;
}

/*!
 * @brief cross product of two vector (RH) and normalize the result
 *
 * @param[in]  a    vector 1
 * @param[in]  b    vector 2
 * @param[out] dest destination
 */
CGLM_INLINE
void
glm_vec3_crossn(vec3 a, vec3 b, vec3 dest) {
  glm_vec3_cross(a, b, dest);
  glm_vec3_normalize(dest);
}

// /*!
//  * @brief angle betwen two vector
//  *
//  * @param[in] a  vector1
//  * @param[in] b  vector2
//  *
//  * @return angle as radians
//  */
// CGLM_INLINE
// float
// glm_vec3_angle(vec3 a, vec3 b) {
//   float norm, dot;

//   /* maybe compiler generate approximation instruction (rcp) */
//   norm = 1.0f / (glm_vec3_norm(a) * glm_vec3_norm(b));
//   dot  = glm_vec3_dot(a, b) * norm;

//   if (dot > 1.0f)
//     return 0.0f;
//   else if (dot < -1.0f)
//     return CGLM_PI;

//   return acosf(dot);
// }

// /*!
//  * @brief rotate vec3 around axis by angle using Rodrigues' rotation formula
//  *
//  * @param[in, out] v     vector
//  * @param[in]      axis  axis vector (must be unit vector)
//  * @param[in]      angle angle by radians
//  */
// CGLM_INLINE
// void
// glm_vec3_rotate(vec3 v, float angle, vec3 axis) {
//   vec3   v1, v2, k;
//   float  c, s;

//   c = cosf(angle);
//   s = sinf(angle);

//   glm_vec3_normalize_to(axis, k);

//   /* Right Hand, Rodrigues' rotation formula:
//         v = v*cos(t) + (kxv)sin(t) + k*(k.v)(1 - cos(t))
//    */
//   glm_vec3_scale(v, c, v1);

//   glm_vec3_cross(k, v, v2);
//   glm_vec3_scale(v2, s, v2);

//   glm_vec3_add(v1, v2, v1);

//   glm_vec3_scale(k, glm_vec3_dot(k, v) * (1.0f - c), v2);
//   glm_vec3_add(v1, v2, v);
// }

// /*!
//  * @brief apply rotation matrix to vector
//  *
//  *  matrix format should be (no perspective):
//  *   a  b  c  x
//  *   e  f  g  y
//  *   i  j  k  z
//  *   0  0  0  w
//  *
//  * @param[in]  m    affine matrix or rot matrix
//  * @param[in]  v    vector
//  * @param[out] dest rotated vector
//  */
// CGLM_INLINE
// void
// glm_vec3_rotate_m4(mat4 m, vec3 v, vec3 dest) {
//   vec4 x, y, z, res;

//   glm_vec4_normalize_to(m->m + 0, x);
//   glm_vec4_normalize_to(m->m + 1, y);
//   glm_vec4_normalize_to(m->m + 2, z);

//   glm_vec4_scale(x,   (*v).x, res);
//   glm_vec4_muladds(y, (*v).y, res);
//   glm_vec4_muladds(z, (*v).z, res);

//   glm_vec3(res, dest);
// }

// /*!
//  * @brief apply rotation matrix to vector
//  *
//  * @param[in]  m    affine matrix or rot matrix
//  * @param[in]  v    vector
//  * @param[out] dest rotated vector
//  */
// CGLM_INLINE
// void
// glm_vec3_rotate_m3(mat3 m, vec3 v, vec3 dest) {
//   vec4 res, x, y, z;

//   glm_vec4((vec3)(m->m + 0), 0.0f, x);
//   glm_vec4((vec3)(m->m + 1), 0.0f, y);
//   glm_vec4((vec3)(m->m + 2), 0.0f, z);

//   glm_vec4_normalize(x);
//   glm_vec4_normalize(y);
//   glm_vec4_normalize(z);

//   glm_vec4_scale(x,   (*v).x, res);
//   glm_vec4_muladds(y, (*v).y, res);
//   glm_vec4_muladds(z, (*v).z, res);

//   glm_vec3(res, dest);
// }

// /*!
//  * @brief project a vector onto b vector
//  *
//  * @param[in]  a    vector1
//  * @param[in]  b    vector2
//  * @param[out] dest projected vector
//  */
// CGLM_INLINE
// void
// glm_vec3_proj(vec3 a, vec3 b, vec3 dest) {
//   glm_vec3_scale(b,
//                  glm_vec3_dot(a, b) / glm_vec3_norm2(b),
//                  dest);
// }

// /**
//  * @brief find center point of two vector
//  *
//  * @param[in]  a    vector1
//  * @param[in]  b    vector2
//  * @param[out] dest center point
//  */
// CGLM_INLINE
// void
// glm_vec3_center(vec3 a, vec3 b, vec3 dest) {
//   glm_vec3_add(a, b, dest);
//   glm_vec3_scale(dest, 0.5f, dest);
// }

// /**
//  * @brief squared distance between two vectors
//  *
//  * @param[in] a vector1
//  * @param[in] b vector2
//  * @return returns squared distance (distance * distance)
//  */
// CGLM_INLINE
// float
// glm_vec3_distance2(vec3 a, vec3 b) {
//   return glm_pow2((*b).x - (*a).x)
//        + glm_pow2((*b).y - (*a).y)
//        + glm_pow2((*b).z - (*a).z);
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
// glm_vec3_distance(vec3 a, vec3 b) {
//   return sqrtf(glm_vec3_distance2(a, b));
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
// glm_vec3_maxv(vec3 a, vec3 b, vec3 dest) {
//   (*dest).x = glm_max((*a).x, (*b).x);
//   (*dest).y = glm_max((*a).y, (*b).y);
//   (*dest).z = glm_max((*a).z, (*b).z);
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
// glm_vec3_minv(vec3 a, vec3 b, vec3 dest) {
//   (*dest).x = glm_min((*a).x, (*b).x);
//   (*dest).y = glm_min((*a).y, (*b).y);
//   (*dest).z = glm_min((*a).z, (*b).z);
// }

// /*!
//  * @brief possible orthogonal/perpendicular vector
//  *
//  * @param[in]  v    vector
//  * @param[out] dest orthogonal/perpendicular vector
//  */
// CGLM_INLINE
// void
// glm_vec3_ortho(vec3 v, vec3 dest) {
//   (*dest).x = (*v).y - (*v).z;
//   (*dest).y = (*v).z - (*v).x;
//   (*dest).z = (*v).x - (*v).y;
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
// glm_vec3_clamp(vec3 v, float minVal, float maxVal) {
//   (*v).x = glm_clamp((*v).x, minVal, maxVal);
//   (*v).y = glm_clamp((*v).y, minVal, maxVal);
//   (*v).z = glm_clamp((*v).z, minVal, maxVal);
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
// glm_vec3_lerp(vec3 from, vec3 to, float t, vec3 dest) {
//   vec3 s, v;

//   /* from + s * (to - from) */
//   glm_vec3_broadcast(glm_clamp_zo(t), s);
//   glm_vec3_sub(to, from, v);
//   glm_vec3_mul(s, v, v);
//   glm_vec3_add(from, v, dest);
// }

// /*!
//  * @brief vec3 cross product
//  *
//  * this is just convenient wrapper
//  *
//  * @param[in]  a source 1
//  * @param[in]  b source 2
//  * @param[out] d destination
//  */
// CGLM_INLINE
// void
// glm_cross(vec3 a, vec3 b, vec3 d) {
//   glm_vec3_cross(a, b, d);
// }

// /*!
//  * @brief vec3 dot product
//  *
//  * this is just convenient wrapper
//  *
//  * @param[in] a vector1
//  * @param[in] b vector2
//  *
//  * @return dot product
//  */
// CGLM_INLINE
// float
// glm_dot(vec3 a, vec3 b) {
//   return glm_vec3_dot(a, b);
// }

// /*!
//  * @brief normalize vec3 and store result in same vec
//  *
//  * this is just convenient wrapper
//  *
//  * @param[in, out] v vector
//  */
// CGLM_INLINE
// void
// glm_normalize(vec3 v) {
//   glm_vec3_normalize(v);
// }

// /*!
//  * @brief normalize vec3 to dest
//  *
//  * this is just convenient wrapper
//  *
//  * @param[in]  v    source
//  * @param[out] dest destination
//  */
// CGLM_INLINE
// void
// glm_normalize_to(vec3 v, vec3 dest) {
//   glm_vec3_normalize_to(v, dest);
// }

#endif /* cglm_vec3_h */
