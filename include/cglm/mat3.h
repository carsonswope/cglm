/*
 * Copyright (c), Recep Aslantas.
 *
 * MIT License (MIT), http://opensource.org/licenses/MIT
 * Full license can be found in the LICENSE file
 */

/*
 Macros:
   GLM_MAT3_IDENTITY_INIT
   GLM_MAT3_ZERO_INIT
   GLM_MAT3_IDENTITY
   GLM_MAT3_ZERO
   glm_mat3_dup(mat, dest)

 Functions:
   CGLM_INLINE void  glm_mat3_copy(mat3 mat, mat3 dest);
   CGLM_INLINE void  glm_mat3_identity(mat3 mat);
   CGLM_INLINE void  glm_mat3_identity_array(mat3 * restrict mat, size_t count);
   CGLM_INLINE void  glm_mat3_zero(mat3 mat);
   CGLM_INLINE void  glm_mat3_mul(mat3 m1, mat3 m2, mat3 dest);
   CGLM_INLINE void  glm_mat3_transpose_to(mat3 m, mat3 dest);
   CGLM_INLINE void  glm_mat3_transpose(mat3 m);
   CGLM_INLINE void  glm_mat3_mulv(mat3 m, vec3 v, vec3 dest);
   CGLM_INLINE float glm_mat3_trace(mat3 m);
   CGLM_INLINE void  glm_mat3_quat(mat3 m, versor dest);
   CGLM_INLINE void  glm_mat3_scale(mat3 m, float s);
   CGLM_INLINE float glm_mat3_det(mat3 mat);
   CGLM_INLINE void  glm_mat3_inv(mat3 mat, mat3 dest);
   CGLM_INLINE void  glm_mat3_swap_col(mat3 mat, int col1, int col2);
   CGLM_INLINE void  glm_mat3_swap_row(mat3 mat, int row1, int row2);
   CGLM_INLINE float glm_mat3_rmc(vec3 r, mat3 m, vec3 c);
 */

#ifndef cglm_mat3_h
#define cglm_mat3_h

#include "common.h"
#include "vec3.h"

#ifdef CGLM_SSE_FP
#  include "simd/sse2/mat3.h"
#endif

#define GLM_MAT3_IDENTITY_INIT  {{{1.0f, 0.0f, 0.0f, 0.0f},                  \
                                  {0.0f, 1.0f, 0.0f, 0.0f},                  \
                                  {0.0f, 0.0f, 1.0f, 0.0f},                  \
                                  {0.0f, 0.0f, 0.0f, 0.0f}}}  // last row: 0s, because it isn't supposed to be 4x4!
#define GLM_MAT3_ZERO_INIT      {{{0.0f, 0.0f, 0.0f, 0.0f},                  \
                                  {0.0f, 0.0f, 0.0f, 0.0f},                  \
                                  {0.0f, 0.0f, 0.0f, 0.0f},                  \
                                  {0.0f, 0.0f, 0.0f, 0.0f}}}


/* for C only */
#define GLM_MAT3_IDENTITY ((_mat3)GLM_MAT3_IDENTITY_INIT)
#define GLM_MAT3_ZERO     ((_mat3)GLM_MAT3_ZERO_INIT)

/* DEPRECATED! use _copy, _ucopy versions */
#define glm_mat3_dup(mat, dest) glm_mat3_copy(mat, dest)

// /*!
//  * @brief copy all members of [mat] to [dest]
//  *
//  * @param[in]  mat  source
//  * @param[out] dest destination
//  */
// CGLM_INLINE
// void
// glm_mat3_copy(mat3 mat, mat3 dest) {
//   dest->m[0].x = mat->m[0].x;
//   dest->m[0].y = mat->m[0].y;
//   dest->m[0].z = mat->m[0].z;

//   dest->m[1].x = mat->m[1].x;
//   dest->m[1].y = mat->m[1].y;
//   dest->m[1].z = mat->m[1].z;

//   dest->m[2].x = mat->m[2].x;
//   dest->m[2].y = mat->m[2].y;
//   dest->m[2].z = mat->m[2].z;
// }

// /*!
//  * @brief make given matrix identity. It is identical with below,
//  *        but it is more easy to do that with this func especially for members
//  *        e.g. glm_mat3_identity(aStruct->aMatrix);
//  *
//  * @code
//  * glm_mat3_copy(GLM_MAT3_IDENTITY, mat); // C only
//  *
//  * // or
//  * mat3 mat = GLM_MAT3_IDENTITY_INIT;
//  * @endcode
//  *
//  * @param[in, out]  mat  destination
//  */
// CGLM_INLINE
// void
// glm_mat3_identity(mat3 mat) {
//   CGLM_ALIGN_MAT _mat3 t = GLM_MAT3_IDENTITY;
//   glm_mat3_copy(&t, mat);
// }

// /*!
//  * @brief make given matrix array's each element identity matrix
//  *
//  * @param[in, out]  mat   matrix array (must be aligned (16/32)
//  *                        if alignment is not disabled)
//  *
//  * @param[in]       count count of matrices
//  */
// CGLM_INLINE
// void
// glm_mat3_identity_array(mat3 * __restrict mat, size_t count) {
//   CGLM_ALIGN_MAT _mat3 t = GLM_MAT3_IDENTITY;
//   size_t i;

//   for (i = 0; i < count; i++) {
//     glm_mat3_copy(&t, mat[i]);
//   }
// }

// /*!
//  * @brief make given matrix zero.
//  *
//  * @param[in, out]  mat  matrix
//  */
// CGLM_INLINE
// void
// glm_mat3_zero(mat3 mat) {
//   CGLM_ALIGN_MAT _mat3 t = GLM_MAT3_ZERO;
//   glm_mat3_copy(&t, mat);
// }

// /*!
//  * @brief multiply m1 and m2 to dest
//  *
//  * m1, m2 and dest matrices can be same matrix, it is possible to write this:
//  *
//  * @code
//  * mat3 m = GLM_MAT3_IDENTITY_INIT;
//  * glm_mat3_mul(m, m, m);
//  * @endcode
//  *
//  * @param[in]  m1   left matrix
//  * @param[in]  m2   right matrix
//  * @param[out] dest destination matrix
//  */
// CGLM_INLINE
// void
// glm_mat3_mul(mat3 m1, mat3 m2, mat3 dest) {
// #if defined( __SSE__ ) || defined( __SSE2__ )
//   glm_mat3_mul_sse2(m1, m2, dest);
// #else
//   float a00 = m1->m[0].x, a01 = m1->m[0].y, a02 = m1->m[0].z,
//         a10 = m1->m[1].x, a11 = m1->m[1].y, a12 = m1->m[1].z,
//         a20 = m1->m[2].x, a21 = m1->m[2].y, a22 = m1->m[2].z,

//         b00 = m2->m[0].x, b01 = m2->m[0].y, b02 = m2->m[0].z,
//         b10 = m2->m[1].x, b11 = m2->m[1].y, b12 = m2->m[1].z,
//         b20 = m2->m[2].x, b21 = m2->m[2].y, b22 = m2->m[2].z;

//   dest->m[0].x = a00 * b00 + a10 * b01 + a20 * b02;
//   dest->m[0].y = a01 * b00 + a11 * b01 + a21 * b02;
//   dest->m[0].z = a02 * b00 + a12 * b01 + a22 * b02;
//   dest->m[1].x = a00 * b10 + a10 * b11 + a20 * b12;
//   dest->m[1].y = a01 * b10 + a11 * b11 + a21 * b12;
//   dest->m[1].z = a02 * b10 + a12 * b11 + a22 * b12;
//   dest->m[2].x = a00 * b20 + a10 * b21 + a20 * b22;
//   dest->m[2].y = a01 * b20 + a11 * b21 + a21 * b22;
//   dest->m[2].z = a02 * b20 + a12 * b21 + a22 * b22;
// #endif
// }

// /*!
//  * @brief transpose mat3 and store in dest
//  *
//  * source matrix will not be transposed unless dest is m
//  *
//  * @param[in]  m     matrix
//  * @param[out] dest  result
//  */
// CGLM_INLINE
// void
// glm_mat3_transpose_to(mat3 m, mat3 dest) {
//   dest->m[0].x = m->m[0].x;
//   dest->m[0].y = m->m[1].x;
//   dest->m[0].z = m->m[2].x;
//   dest->m[1].x = m->m[0].y;
//   dest->m[1].y = m->m[1].y;
//   dest->m[1].z = m->m[2].y;
//   dest->m[2].x = m->m[0].z;
//   dest->m[2].y = m->m[1].z;
//   dest->m[2].z = m->m[2].z;
// }

// /*!
//  * @brief tranpose mat3 and store result in same matrix
//  *
//  * @param[in, out] m source and dest
//  */
// CGLM_INLINE
// void
// glm_mat3_transpose(mat3 m) {
//   CGLM_ALIGN_MAT mat3 tmp;

//   tmp->m[0].y = m->m[1].x;
//   tmp->m[0].z = m->m[2].x;
//   tmp->m[1].x = m->m[0].y;
//   tmp->m[1].z = m->m[2].y;
//   tmp->m[2].x = m->m[0].z;
//   tmp->m[2].y = m->m[1].z;

//   m->m[0].y = tmp->m[0].y;
//   m->m[0].z = tmp->m[0].z;
//   m->m[1].x = tmp->m[1].x;
//   m->m[1].z = tmp->m[1].z;
//   m->m[2].x = tmp->m[2].x;
//   m->m[2].y = tmp->m[2].y;
// }

// /*!
//  * @brief multiply mat3 with vec3 (column vector) and store in dest vector
//  *
//  * @param[in]  m    mat3 (left)
//  * @param[in]  v    vec3 (right, column vector)
//  * @param[out] dest vec3 (result, column vector)
//  */
// CGLM_INLINE
// void
// glm_mat3_mulv(mat3 m, vec3 v, vec3 dest) {
//   dest[0] = m->m[0].x * v[0] + m->m[1].x * v[1] + m->m[2].x * v[2];
//   dest[1] = m->m[0].y * v[0] + m->m[1].y * v[1] + m->m[2].y * v[2];
//   dest[2] = m->m[0].z * v[0] + m->m[1].z * v[1] + m->m[2].z * v[2];
// }

// /*!
//  * @brief trace of matrix
//  *
//  * sum of the elements on the main diagonal from upper left to the lower right
//  *
//  * @param[in]  m matrix
//  */
// CGLM_INLINE
// float
// glm_mat3_trace(mat3 m) {
//   return m->m[0].x + m->m[1].y + m->m[2].z;
// }

// /*!
//  * @brief convert mat3 to quaternion
//  *
//  * @param[in]  m    rotation matrix
//  * @param[out] dest destination quaternion
//  */
// CGLM_INLINE
// void
// glm_mat3_quat(mat3 m, versor dest) {
//   float trace, r, rinv;

//   /* it seems using like m12 instead of m->m[1].z causes extra instructions */

//   trace = m->m[0].x + m->m[1].y + m->m[2].z;
//   if (trace >= 0.0f) {
//     r       = sqrtf(1.0f + trace);
//     rinv    = 0.5f / r;

//     dest[0] = rinv * (m->m[1].z - m->m[2].y);
//     dest[1] = rinv * (m->m[2].x - m->m[0].z);
//     dest[2] = rinv * (m->m[0].y - m->m[1].x);
//     dest[3] = r    * 0.5f;
//   } else if (m->m[0].x >= m->m[1].y && m->m[0].x >= m->m[2].z) {
//     r       = sqrtf(1.0f - m->m[1].y - m->m[2].z + m->m[0].x);
//     rinv    = 0.5f / r;

//     dest[0] = r    * 0.5f;
//     dest[1] = rinv * (m->m[0].y + m->m[1].x);
//     dest[2] = rinv * (m->m[0].z + m->m[2].x);
//     dest[3] = rinv * (m->m[1].z - m->m[2].y);
//   } else if (m->m[1].y >= m->m[2].z) {
//     r       = sqrtf(1.0f - m->m[0].x - m->m[2].z + m->m[1].y);
//     rinv    = 0.5f / r;

//     dest[0] = rinv * (m->m[0].y + m->m[1].x);
//     dest[1] = r    * 0.5f;
//     dest[2] = rinv * (m->m[1].z + m->m[2].y);
//     dest[3] = rinv * (m->m[2].x - m->m[0].z);
//   } else {
//     r       = sqrtf(1.0f - m->m[0].x - m->m[1].y + m->m[2].z);
//     rinv    = 0.5f / r;

//     dest[0] = rinv * (m->m[0].z + m->m[2].x);
//     dest[1] = rinv * (m->m[1].z + m->m[2].y);
//     dest[2] = r    * 0.5f;
//     dest[3] = rinv * (m->m[0].y - m->m[1].x);
//   }
// }

// /*!
//  * @brief scale (multiply with scalar) matrix
//  *
//  * multiply matrix with scalar
//  *
//  * @param[in, out] m matrix
//  * @param[in]      s scalar
//  */
// CGLM_INLINE
// void
// glm_mat3_scale(mat3 m, float s) {
//   m->m[0].x *= s; m->m[0].y *= s; m->m[0].z *= s;
//   m->m[1].x *= s; m->m[1].y *= s; m->m[1].z *= s;
//   m->m[2].x *= s; m->m[2].y *= s; m->m[2].z *= s;
// }

// /*!
//  * @brief mat3 determinant
//  *
//  * @param[in] mat matrix
//  *
//  * @return determinant
//  */
// CGLM_INLINE
// float
// glm_mat3_det(mat3 mat) {
//   float a = mat->m[0].x, b = mat->m[0].y, c = mat->m[0].z,
//         d = mat->m[1].x, e = mat->m[1].y, f = mat->m[1].z,
//         g = mat->m[2].x, h = mat->m[2].y, i = mat->m[2].z;

//   return a * (e * i - h * f) - d * (b * i - c * h) + g * (b * f - c * e);
// }

// /*!
//  * @brief inverse mat3 and store in dest
//  *
//  * @param[in]  mat  matrix
//  * @param[out] dest inverse matrix
//  */
// CGLM_INLINE
// void
// glm_mat3_inv(mat3 mat, mat3 dest) {
//   float det;
//   float a = mat->m[0].x, b = mat->m[0].y, c = mat->m[0].z,
//         d = mat->m[1].x, e = mat->m[1].y, f = mat->m[1].z,
//         g = mat->m[2].x, h = mat->m[2].y, i = mat->m[2].z;

//   dest->m[0].x =   e * i - f * h;
//   dest->m[0].y = -(b * i - h * c);
//   dest->m[0].z =   b * f - e * c;
//   dest->m[1].x = -(d * i - g * f);
//   dest->m[1].y =   a * i - c * g;
//   dest->m[1].z = -(a * f - d * c);
//   dest->m[2].x =   d * h - g * e;
//   dest->m[2].y = -(a * h - g * b);
//   dest->m[2].z =   a * e - b * d;

//   det = 1.0f / (a * dest->m[0].x + b * dest->m[1].x + c * dest->m[2].x);

//   glm_mat3_scale(dest, det);
// }

// /*!
//  * @brief swap two matrix columns
//  *
//  * @param[in,out] mat  matrix
//  * @param[in]     col1 col1
//  * @param[in]     col2 col2
//  */
// CGLM_INLINE
// void
// glm_mat3_swap_col(mat3 mat, int col1, int col2) {
//   float3 tmp;
//   //vec3 tmp;
//   glm_vec3_copy((vec3)mat->m + col1, &tmp);
//   glm_vec3_copy((vec3)mat->m + col2, (vec3)mat->m + col1);
//   glm_vec3_copy(&tmp, (vec3)mat-> m + col2);
// }

// /*!
//  * @brief swap two matrix rows
//  *
//  * @param[in,out] mat  matrix
//  * @param[in]     row1 row1
//  * @param[in]     row2 row2
//  */
// CGLM_INLINE
// void
// glm_mat3_swap_row(mat3 mat, int row1, int row2) {
//   float3 tmp;
//   tmp.x = ((float*)mat->m + 0)[row1];
//   tmp.y = ((float*)mat->m + 1)[row1];
//   tmp.z = ((float*)mat->m + 2)[row1];

//   ((float*)mat->m + 0)[row1] = ((float*)mat->m + 0)[row2];
//   ((float*)mat->m + 1)[row1] = ((float*)mat->m + 1)[row2];
//   ((float*)mat->m + 2)[row1] = ((float*)mat->m + 2)[row2];

//   ((float*)mat->m + 0)[row2] = tmp.x;
//   ((float*)mat->m + 1)[row2] = tmp.y;
//   ((float*)mat->m + 2)[row2] = tmp.z;
// }

// /*!
//  * @brief helper for  R (row vector) * M (matrix) * C (column vector)
//  *
//  * rmc stands for Row * Matrix * Column
//  *
//  * the result is scalar because R * M = Matrix1x3 (row vector),
//  * then Matrix1x3 * Vec3 (column vector) = Matrix1x1 (Scalar)
//  *
//  * @param[in]  r   row vector or matrix1x3
//  * @param[in]  m   matrix3x3
//  * @param[in]  c   column vector or matrix3x1
//  *
//  * @return scalar value e.g. Matrix1x1
//  */
// CGLM_INLINE
// float
// glm_mat3_rmc(vec3 r, mat3 m, vec3 c) {
//   vec3 tmp;
//   glm_mat3_mulv(m, c, tmp);
//   return glm_vec3_dot(r, tmp);
// }

#endif /* cglm_mat3_h */
