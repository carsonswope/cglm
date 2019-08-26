/*
 * Copyright (c), Recep Aslantas.
 *
 * MIT License (MIT), http://opensource.org/licenses/MIT
 * Full license can be found in the LICENSE file
 */

/*!
 * @brief SIMD like functions
 */

/*
 Functions:
   CGLM_INLINE void  glm_vec3_broadcast(float val, vec3 d);
   CGLM_INLINE bool  glm_vec3_eq(vec3 v, float val);
   CGLM_INLINE bool  glm_vec3_eq_eps(vec3 v, float val);
   CGLM_INLINE bool  glm_vec3_eq_all(vec3 v);
   CGLM_INLINE bool  glm_vec3_eqv(vec3 a, vec3 b);
   CGLM_INLINE bool  glm_vec3_eqv_eps(vec3 a, vec3 b);
   CGLM_INLINE float glm_vec3_max(vec3 v);
   CGLM_INLINE float glm_vec3_min(vec3 v);
   CGLM_INLINE bool  glm_vec3_isnan(vec3 v);
   CGLM_INLINE bool  glm_vec3_isinf(vec3 v);
   CGLM_INLINE bool  glm_vec3_isvalid(vec3 v);
   CGLM_INLINE void  glm_vec3_sign(vec3 v, vec3 dest);
   CGLM_INLINE void  glm_vec3_sqrt(vec3 v, vec3 dest);
 */

#ifndef cglm_vec3_ext_h
#define cglm_vec3_ext_h

#include "common.h"
#include "util.h"

/*!
 * @brief fill a vector with specified value
 *
 * @param[in]  val value
 * @param[out] d   dest
 */
CGLM_INLINE
void
glm_vec3_broadcast(float val, vec3 d) {
  d->x = d->y = d->z = val;
}

/*!
 * @brief check if vector is equal to value (without epsilon)
 *
 * @param[in] v   vector
 * @param[in] val value
 */
CGLM_INLINE
bool
glm_vec3_eq(vec3 v, float val) {
  return v->x == val && v->x == v->y && v->x == v->z;
}

/*!
 * @brief check if vector is equal to value (with epsilon)
 *
 * @param[in] v   vector
 * @param[in] val value
 */
CGLM_INLINE
bool
glm_vec3_eq_eps(vec3 v, float val) {
  return fabsf(v->x - val) <= FLT_EPSILON
         && fabsf(v->y - val) <= FLT_EPSILON
         && fabsf(v->z - val) <= FLT_EPSILON;
}

/*!
 * @brief check if vectors members are equal (without epsilon)
 *
 * @param[in] v   vector
 */
CGLM_INLINE
bool
glm_vec3_eq_all(vec3 v) {
  return v->x == v->y && v->x == v->z;
}

/*!
 * @brief check if vector is equal to another (without epsilon)
 *
 * @param[in] a vector
 * @param[in] b vector
 */
CGLM_INLINE
bool
glm_vec3_eqv(vec3 a, vec3 b) {
  return a->x == b->x
         && a->y == b->y
         && a->z == b->z;
}

/*!
 * @brief check if vector is equal to another (with epsilon)
 *
 * @param[in] a vector
 * @param[in] b vector
 */
CGLM_INLINE
bool
glm_vec3_eqv_eps(vec3 a, vec3 b) {
  return fabsf(a->x - b->x) <= FLT_EPSILON
         && fabsf(a->y - b->y) <= FLT_EPSILON
         && fabsf(a->z - b->z) <= FLT_EPSILON;
}

/*!
 * @brief max value of vector
 *
 * @param[in] v vector
 */
CGLM_INLINE
float
glm_vec3_max(vec3 v) {
  float max;

  max = v->x;
  if (v->y > max)
    max = v->y;
  if (v->z > max)
    max = v->z;

  return max;
}

/*!
 * @brief min value of vector
 *
 * @param[in] v vector
 */
CGLM_INLINE
float
glm_vec3_min(vec3 v) {
  float min;

  min = v->x;
  if (v->y < min)
    min = v->y;
  if (v->z < min)
    min = v->z;

  return min;
}

/*!
 * @brief check if all items are NaN (not a number)
 *        you should only use this in DEBUG mode or very critical asserts
 *
 * @param[in] v vector
 */
CGLM_INLINE
bool
glm_vec3_isnan(vec3 v) {
  return isnan(v->x) || isnan(v->y) || isnan(v->z);
}

/*!
 * @brief check if all items are INFINITY
 *        you should only use this in DEBUG mode or very critical asserts
 *
 * @param[in] v vector
 */
CGLM_INLINE
bool
glm_vec3_isinf(vec3 v) {
  return isinf(v->x) || isinf(v->y) || isinf(v->z);
}

/*!
 * @brief check if all items are valid number
 *        you should only use this in DEBUG mode or very critical asserts
 *
 * @param[in] v vector
 */
CGLM_INLINE
bool
glm_vec3_isvalid(vec3 v) {
  return !glm_vec3_isnan(v) && !glm_vec3_isinf(v);
}

/*!
 * @brief get sign of 32 bit float as +1, -1, 0
 *
 * Important: It returns 0 for zero/NaN input
 *
 * @param v vector
 */
CGLM_INLINE
void
glm_vec3_sign(vec3 v, vec3 dest) {
  dest->x = glm_signf(v->x);
  dest->y = glm_signf(v->y);
  dest->z = glm_signf(v->z);
}

/*!
 * @brief square root of each vector item
 *
 * @param[in]  v    vector
 * @param[out] dest destination vector
 */
CGLM_INLINE
void
glm_vec3_sqrt(vec3 v, vec3 dest) {
  dest->x = sqrtf(v->x);
  dest->y = sqrtf(v->y);
  dest->z = sqrtf(v->z);
}

#endif /* cglm_vec3_ext_h */
