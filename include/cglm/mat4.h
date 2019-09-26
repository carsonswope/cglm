/*
 * Copyright (c), Recep Aslantas.
 *
 * MIT License (MIT), http://opensource.org/licenses/MIT
 * Full license can be found in the LICENSE file
 */

/*!
 * Most of functions in this header are optimized manually with SIMD
 * if available. You dont need to call/incude SIMD headers manually
 */

/*
 Macros:
   GLM_MAT4_IDENTITY_INIT
   GLM_MAT4_ZERO_INIT
   GLM_MAT4_IDENTITY
   GLM_MAT4_ZERO

 Functions:
   CGLM_INLINE void  glm_mat4_ucopy(mat4 mat, mat4 dest);
   CGLM_INLINE void  glm_mat4_copy(mat4 mat, mat4 dest);
   CGLM_INLINE void  glm_mat4_identity(mat4 mat);
   CGLM_INLINE void  glm_mat4_identity_array(mat4 * restrict mat, size_t count);
   CGLM_INLINE void  glm_mat4_zero(mat4 mat);
   CGLM_INLINE void  glm_mat4_pick3(mat4 mat, mat3 dest);
   CGLM_INLINE void  glm_mat4_pick3t(mat4 mat, mat3 dest);
   CGLM_INLINE void  glm_mat4_ins3(mat3 mat, mat4 dest);
   CGLM_INLINE void  glm_mat4_mul(mat4 m1, mat4 m2, mat4 dest);
   CGLM_INLINE void  glm_mat4_mulN(mat4 *matrices[], int len, mat4 dest);
   CGLM_INLINE void  glm_mat4_mulv(mat4 m, vec4 v, vec4 dest);
   CGLM_INLINE void  glm_mat4_mulv3(mat4 m, vec3 v, vec3 dest);
   CGLM_INLINE float glm_mat4_trace(mat4 m);
   CGLM_INLINE float glm_mat4_trace3(mat4 m);
   CGLM_INLINE void  glm_mat4_quat(mat4 m, versor dest) ;
   CGLM_INLINE void  glm_mat4_transpose_to(mat4 m, mat4 dest);
   CGLM_INLINE void  glm_mat4_transpose(mat4 m);
   CGLM_INLINE void  glm_mat4_scale_p(mat4 m, float s);
   CGLM_INLINE void  glm_mat4_scale(mat4 m, float s);
   CGLM_INLINE float glm_mat4_det(mat4 mat);
   CGLM_INLINE void  glm_mat4_inv(mat4 mat, mat4 dest);
   CGLM_INLINE void  glm_mat4_inv_fast(mat4 mat, mat4 dest);
   CGLM_INLINE void  glm_mat4_swap_col(mat4 mat, int col1, int col2);
   CGLM_INLINE void  glm_mat4_swap_row(mat4 mat, int row1, int row2);
   CGLM_INLINE float glm_mat4_rmc(vec4 r, mat4 m, vec4 c);
 */

#ifndef cglm_mat_h
#define cglm_mat_h

#include "common.h"
#include "vec4.h"
#include "vec3.h"

#ifdef CGLM_SSE_FP
#  include "simd/sse2/mat4.h"
#endif

#ifdef CGLM_AVX_FP
#  include "simd/avx/mat4.h"
#endif

#ifdef CGLM_NEON_FP
#  include "simd/neon/mat4.h"
#endif

#ifdef DEBUG
# include <assert.h>
#endif

#define GLM_MAT4_IDENTITY_INIT  {{(float4){1.0f, 0.0f, 0.0f, 0.0f},                    \
                                 (float4){0.0f, 1.0f, 0.0f, 0.0f},                    \
                                 (float4){0.0f, 0.0f, 1.0f, 0.0f},                    \
                                 (float4){0.0f, 0.0f, 0.0f, 1.0f}}}

#define GLM_MAT4_ZERO_INIT      {{(float4){0.0f, 0.0f, 0.0f, 0.0f},                    \
                                 (float4){0.0f, 0.0f, 0.0f, 0.0f},                    \
                                 (float4){0.0f, 0.0f, 0.0f, 0.0f},                    \
                                 (float4){0.0f, 0.0f, 0.0f, 0.0f}}}

/* for C only */
#define GLM_MAT4_IDENTITY ((_mat4)GLM_MAT4_IDENTITY_INIT)
#define GLM_MAT4_ZERO     ((_mat4)GLM_MAT4_ZERO_INIT)

/* DEPRECATED! use _copy, _ucopy versions */
#define glm_mat4_udup(mat, dest) glm_mat4_ucopy(mat, dest)
#define glm_mat4_dup(mat, dest)  glm_mat4_copy(mat, dest)

/* DEPRECATED! default is precise now. */
#define glm_mat4_inv_precise(mat, dest) glm_mat4_inv(mat, dest)

/*!
 * @brief copy all members of [mat] to [dest]
 *
 * matrix may not be aligned, u stands for unaligned, this may be useful when
 * copying a matrix from external source e.g. asset importer...
 *
 * @param[in]  mat  source
 * @param[out] dest destination
 */
CGLM_INLINE
void glm_mat4_ucopy(mat4 mat, mat4 dest) {
  dest->m[0].x = mat->m[0].x;  dest->m[1].x = mat->m[1].x;
  dest->m[0].y = mat->m[0].y;  dest->m[1].y = mat->m[1].y;
  dest->m[0].z = mat->m[0].z;  dest->m[1].z = mat->m[1].z;
  dest->m[0].w = mat->m[0].w;  dest->m[1].w = mat->m[1].w;

  dest->m[2].x = mat->m[2].x;  dest->m[3].x = mat->m[3].x;
  dest->m[2].y = mat->m[2].y;  dest->m[3].y = mat->m[3].y;
  dest->m[2].z = mat->m[2].z;  dest->m[3].z = mat->m[3].z;
  dest->m[2].w = mat->m[2].w;  dest->m[3].w = mat->m[3].w;
}

/*!
 * @brief copy all members of [mat] to [dest]
 *
 * @param[in]  mat  source
 * @param[out] dest destination
 */
CGLM_INLINE
void
glm_mat4_copy(mat4 mat, mat4 dest) {
  glm_mat4_ucopy(mat, dest);
}

/*!
 * @brief make given matrix identity. It is identical with below, 
 *        but it is more easy to do that with this func especially for members
 *        e.g. glm_mat4_identity(aStruct->aMatrix);
 *
 * @code
 * glm_mat4_copy(GLM_MAT4_IDENTITY, mat); // C only
 *
 * // or
 * mat4 mat = GLM_MAT4_IDENTITY_INIT;
 * @endcode
 *
 * @param[in, out]  mat  destination
 */
CGLM_INLINE
void
glm_mat4_identity(mat4 mat) {
  CGLM_ALIGN_MAT _mat4 t = GLM_MAT4_IDENTITY;
  glm_mat4_copy(&t, mat);
}

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
// glm_mat4_identity_array(mat4 * __restrict mat, size_t count) {
//   CGLM_ALIGN_MAT _mat4 t = GLM_MAT4_IDENTITY;
//   size_t i;

//   for (i = 0; i < count; i++) {
//     glm_mat4_copy(&t, mat[i]);
//   }
// }

// /*!
//  * @brief make given matrix zero.
//  *
//  * @param[in, out]  mat  matrix
//  */
// CGLM_INLINE
// void
// glm_mat4_zero(mat4 mat) {
//   CGLM_ALIGN_MAT _mat4 t = GLM_MAT4_ZERO;
//   glm_mat4_copy(&t, mat);
// }

// /*!
//  * @brief copy upper-left of mat4 to mat3
//  *
//  * @param[in]  mat  source
//  * @param[out] dest destination
//  */
// CGLM_INLINE
// void
// glm_mat4_pick3(mat4 mat, mat3 dest) {
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
//  * @brief copy upper-left of mat4 to mat3 (transposed)
//  *
//  * the postfix t stands for transpose
//  *
//  * @param[in]  mat  source
//  * @param[out] dest destination
//  */
// CGLM_INLINE
// void
// glm_mat4_pick3t(mat4 mat, mat3 dest) {
//   dest->m[0].x = mat->m[0].x;
//   dest->m[0].y = mat->m[1].x;
//   dest->m[0].z = mat->m[2].x;

//   dest->m[1].x = mat->m[0].y;
//   dest->m[1].y = mat->m[1].y;
//   dest->m[1].z = mat->m[2].y;

//   dest->m[2].x = mat->m[0].z;
//   dest->m[2].y = mat->m[1].z;
//   dest->m[2].z = mat->m[2].z;
// }

// /*!
//  * @brief copy mat3 to mat4's upper-left
//  *
//  * @param[in]  mat  source
//  * @param[out] dest destination
//  */
// CGLM_INLINE
// void
// glm_mat4_ins3(mat3 mat, mat4 dest) {
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

/*!
 * @brief multiply m1 and m2 to dest
 *
 * m1, m2 and dest matrices can be same matrix, it is possible to write this:
 *
 * @code
 * mat4 m = GLM_MAT4_IDENTITY_INIT;
 * glm_mat4_mul(m, m, m);
 * @endcode
 *
 * @param[in]  m1   left matrix
 * @param[in]  m2   right matrix
 * @param[out] dest destination matrix
 */
CGLM_INLINE
void glm_mat4_mul(mat4 m1, mat4 m2, mat4 dest) {
  float a00 = m1->m[0].x, a01 = m1->m[0].y, a02 = m1->m[0].z, a03 = m1->m[0].w,
        a10 = m1->m[1].x, a11 = m1->m[1].y, a12 = m1->m[1].z, a13 = m1->m[1].w,
        a20 = m1->m[2].x, a21 = m1->m[2].y, a22 = m1->m[2].z, a23 = m1->m[2].w,
        a30 = m1->m[3].x, a31 = m1->m[3].y, a32 = m1->m[3].z, a33 = m1->m[3].w,

        b00 = m2->m[0].x, b01 = m2->m[0].y, b02 = m2->m[0].z, b03 = m2->m[0].w,
        b10 = m2->m[1].x, b11 = m2->m[1].y, b12 = m2->m[1].z, b13 = m2->m[1].w,
        b20 = m2->m[2].x, b21 = m2->m[2].y, b22 = m2->m[2].z, b23 = m2->m[2].w,
        b30 = m2->m[3].x, b31 = m2->m[3].y, b32 = m2->m[3].z, b33 = m2->m[3].w;

  dest->m[0].x = a00 * b00 + a10 * b01 + a20 * b02 + a30 * b03;
  dest->m[0].y = a01 * b00 + a11 * b01 + a21 * b02 + a31 * b03;
  dest->m[0].z = a02 * b00 + a12 * b01 + a22 * b02 + a32 * b03;
  dest->m[0].w = a03 * b00 + a13 * b01 + a23 * b02 + a33 * b03;
  dest->m[1].x = a00 * b10 + a10 * b11 + a20 * b12 + a30 * b13;
  dest->m[1].y = a01 * b10 + a11 * b11 + a21 * b12 + a31 * b13;
  dest->m[1].z = a02 * b10 + a12 * b11 + a22 * b12 + a32 * b13;
  dest->m[1].w = a03 * b10 + a13 * b11 + a23 * b12 + a33 * b13;
  dest->m[2].x = a00 * b20 + a10 * b21 + a20 * b22 + a30 * b23;
  dest->m[2].y = a01 * b20 + a11 * b21 + a21 * b22 + a31 * b23;
  dest->m[2].z = a02 * b20 + a12 * b21 + a22 * b22 + a32 * b23;
  dest->m[2].w = a03 * b20 + a13 * b21 + a23 * b22 + a33 * b23;
  dest->m[3].x = a00 * b30 + a10 * b31 + a20 * b32 + a30 * b33;
  dest->m[3].y = a01 * b30 + a11 * b31 + a21 * b32 + a31 * b33;
  dest->m[3].z = a02 * b30 + a12 * b31 + a22 * b32 + a32 * b33;
  dest->m[3].w = a03 * b30 + a13 * b31 + a23 * b32 + a33 * b33;
}

// /*!
//  * @brief mupliply N mat4 matrices and store result in dest
//  *
//  * this function lets you multiply multiple (more than two or more...) matrices
//  * <br><br>multiplication will be done in loop, this may reduce instructions
//  * size but if <b>len</b> is too small then compiler may unroll whole loop,
//  * usage:
//  * @code
//  * mat m1, m2, m3, m4, res;
//  *
//  * glm_mat4_mulN((mat4 *[]){&m1, &m2, &m3, &m4}, 4, res);
//  * @endcode
//  *
//  * @warning matrices parameter is pointer array not mat4 array!
//  *
//  * @param[in]  matrices mat4 * array
//  * @param[in]  len      matrices count
//  * @param[out] dest     result
//  */
// CGLM_INLINE
// void
// glm_mat4_mulN(mat4 * __restrict matrices[], uint32_t len, mat4 dest) {
//   uint32_t i;

// #ifdef DEBUG
//   assert(len > 1 && "there must be least 2 matrices to go!");
// #endif

//   glm_mat4_mul(*matrices[0], *matrices[1], dest);

//   for (i = 2; i < len; i++)
//     glm_mat4_mul(dest, *matrices[i], dest);
// }

/*!
 * @brief multiply mat4 with vec4 (column vector) and store in dest vector
 *
 * @param[in]  m    mat4 (left)
 * @param[in]  v    vec4 (right, column vector)
 * @param[out] dest vec4 (result, column vector)
 */
CGLM_INLINE
void
glm_mat4_mulv(mat4 m, vec4 v, vec4 dest) {
  float4 res;
  res.x = m->m[0].x * (*v).x + m->m[1].x * (*v).y + m->m[2].x * (*v).z + m->m[3].x * (*v).w;
  res.y = m->m[0].y * (*v).x + m->m[1].y * (*v).y + m->m[2].y * (*v).z + m->m[3].y * (*v).w;
  res.z = m->m[0].z * (*v).x + m->m[1].z * (*v).y + m->m[2].z * (*v).z + m->m[3].z * (*v).w;
  res.w = m->m[0].w * (*v).x + m->m[1].w * (*v).y + m->m[2].w * (*v).z + m->m[3].w * (*v).w;
  glm_vec4_copy(&res, dest);
}

// /*!
//  * @brief trace of matrix
//  *
//  * sum of the elements on the main diagonal from upper left to the lower right
//  *
//  * @param[in]  m matrix
//  */
// CGLM_INLINE
// float
// glm_mat4_trace(mat4 m) {
//   return m->m[0].x + m->m[1].y + m->m[2].z + m->m[3].w;
// }

// /*!
//  * @brief trace of matrix (rotation part)
//  *
//  * sum of the elements on the main diagonal from upper left to the lower right
//  *
//  * @param[in]  m matrix
//  */
// CGLM_INLINE
// float
// glm_mat4_trace3(mat4 m) {
//   return m->m[0].x + m->m[1].y + m->m[2].z;
// }

// /*!
//  * @brief convert mat4's rotation part to quaternion
//  *
//  * @param[in]  m    affine matrix
//  * @param[out] dest destination quaternion
//  */
// CGLM_INLINE
// void
// glm_mat4_quat(mat4 m, versor dest) {
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
//  * @brief multiply vector with mat4
//  *
//  * @param[in]  m    mat4(affine transform)
//  * @param[in]  v    vec3
//  * @param[in]  last 4th item to make it vec4
//  * @param[out] dest result vector (vec3)
//  */
// CGLM_INLINE
// void
// glm_mat4_mulv3(mat4 m, vec3 v, float last, vec3 dest) {
//   vec4 res;
//   glm_vec4(v, last, res);
//   glm_mat4_mulv(m, res, res);
//   glm_vec3(res, dest);
// }

/*!
 * @brief transpose mat4 and store in dest
 *
 * source matrix will not be transposed unless dest is m
 *
 * @param[in]  m    matrix
 * @param[out] dest result
 */
CGLM_INLINE
void
glm_mat4_transpose_to(mat4 m, mat4 dest) {
  dest->m[0].x = m->m[0].x; dest->m[1].x = m->m[0].y;
  dest->m[0].y = m->m[1].x; dest->m[1].y = m->m[1].y;
  dest->m[0].z = m->m[2].x; dest->m[1].z = m->m[2].y;
  dest->m[0].w = m->m[3].x; dest->m[1].w = m->m[3].y;
  dest->m[2].x = m->m[0].z; dest->m[3].x = m->m[0].w;
  dest->m[2].y = m->m[1].z; dest->m[3].y = m->m[1].w;
  dest->m[2].z = m->m[2].z; dest->m[3].z = m->m[2].w;
  dest->m[2].w = m->m[3].z; dest->m[3].w = m->m[3].w;
}

/*!
 * @brief tranpose mat4 and store result in same matrix
 *
 * @param[in, out] m source and dest
 */
CGLM_INLINE
void glm_mat4_transpose(mat4 m) {
  _mat4 d;
  glm_mat4_transpose_to(m, &d);
  glm_mat4_ucopy(&d, m);
}

// /*!
//  * @brief scale (multiply with scalar) matrix without simd optimization
//  *
//  * multiply matrix with scalar
//  *
//  * @param[in, out] m matrix
//  * @param[in]      s scalar
//  */
// CGLM_INLINE
// void
// glm_mat4_scale_p(mat4 m, float s) {
//   m->m[0].x *= s; m->m[0].y *= s; m->m[0].z *= s; m->m[0].w *= s;
//   m->m[1].x *= s; m->m[1].y *= s; m->m[1].z *= s; m->m[1].w *= s;
//   m->m[2].x *= s; m->m[2].y *= s; m->m[2].z *= s; m->m[2].w *= s;
//   m->m[3].x *= s; m->m[3].y *= s; m->m[3].z *= s; m->m[3].w *= s;
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
// glm_mat4_scale(mat4 m, float s) {
// #if defined( __SSE__ ) || defined( __SSE2__ )
//   glm_mat4_scale_sse2(m, s);
// #elif defined(CGLM_NEON_FP)
//   float32x4_t v0;
//   v0 = vdupq_n_f32(s);
//   vst1q_f32(m[0], vmulq_f32(vld1q_f32(m[0]), v0));
//   vst1q_f32(m[1], vmulq_f32(vld1q_f32(m[1]), v0));
//   vst1q_f32(m[2], vmulq_f32(vld1q_f32(m[2]), v0));
//   vst1q_f32(m[3], vmulq_f32(vld1q_f32(m[3]), v0));
// #else
//   glm_mat4_scale_p(m, s);
// #endif
// }

// /*!
//  * @brief mat4 determinant
//  *
//  * @param[in] mat matrix
//  *
//  * @return determinant
//  */
// CGLM_INLINE
// float
// glm_mat4_det(mat4 mat) {
// #if defined( __SSE__ ) || defined( __SSE2__ )
//   return glm_mat4_det_sse2(mat);
// #else
//   /* [square] det(A) = det(At) */
//   float t[6];
//   float a = mat->m[0].x, b = mat->m[0].y, c = mat->m[0].z, d = mat->m[0].w,
//         e = mat->m[1].x, f = mat->m[1].y, g = mat->m[1].z, h = mat->m[1].w,
//         i = mat->m[2].x, j = mat->m[2].y, k = mat->m[2].z, l = mat->m[2].w,
//         m = mat->m[3].x, n = mat->m[3].y, o = mat->m[3].z, p = mat->m[3].w;

//   t[0] = k * p - o * l;
//   t[1] = j * p - n * l;
//   t[2] = j * o - n * k;
//   t[3] = i * p - m * l;
//   t[4] = i * o - m * k;
//   t[5] = i * n - m * j;

//   return a * (f * t[0] - g * t[1] + h * t[2])
//        - b * (e * t[0] - g * t[3] + h * t[4])
//        + c * (e * t[1] - f * t[3] + h * t[5])
//        - d * (e * t[2] - f * t[4] + g * t[5]);
// #endif
// }

// /*!
//  * @brief inverse mat4 and store in dest
//  *
//  * @param[in]  mat  matrix
//  * @param[out] dest inverse matrix
//  */
// CGLM_INLINE
// void
// glm_mat4_inv(mat4 mat, mat4 dest) {
// #if defined( __SSE__ ) || defined( __SSE2__ )
//   glm_mat4_inv_sse2(mat, dest);
// #else
//   float t[6];
//   float det;
//   float a = mat->m[0].x, b = mat->m[0].y, c = mat->m[0].z, d = mat->m[0].w,
//         e = mat->m[1].x, f = mat->m[1].y, g = mat->m[1].z, h = mat->m[1].w,
//         i = mat->m[2].x, j = mat->m[2].y, k = mat->m[2].z, l = mat->m[2].w,
//         m = mat->m[3].x, n = mat->m[3].y, o = mat->m[3].z, p = mat->m[3].w;

//   t[0] = k * p - o * l; t[1] = j * p - n * l; t[2] = j * o - n * k;
//   t[3] = i * p - m * l; t[4] = i * o - m * k; t[5] = i * n - m * j;

//   dest->m[0].x =  f * t[0] - g * t[1] + h * t[2];
//   dest->m[1].x =-(e * t[0] - g * t[3] + h * t[4]);
//   dest->m[2].x =  e * t[1] - f * t[3] + h * t[5];
//   dest->m[3].x =-(e * t[2] - f * t[4] + g * t[5]);

//   dest->m[0].y =-(b * t[0] - c * t[1] + d * t[2]);
//   dest->m[1].y =  a * t[0] - c * t[3] + d * t[4];
//   dest->m[2].y =-(a * t[1] - b * t[3] + d * t[5]);
//   dest->m[3].y =  a * t[2] - b * t[4] + c * t[5];

//   t[0] = g * p - o * h; t[1] = f * p - n * h; t[2] = f * o - n * g;
//   t[3] = e * p - m * h; t[4] = e * o - m * g; t[5] = e * n - m * f;

//   dest->m[0].z =  b * t[0] - c * t[1] + d * t[2];
//   dest->m[1].z =-(a * t[0] - c * t[3] + d * t[4]);
//   dest->m[2].z =  a * t[1] - b * t[3] + d * t[5];
//   dest->m[3].z =-(a * t[2] - b * t[4] + c * t[5]);

//   t[0] = g * l - k * h; t[1] = f * l - j * h; t[2] = f * k - j * g;
//   t[3] = e * l - i * h; t[4] = e * k - i * g; t[5] = e * j - i * f;

//   dest->m[0].w =-(b * t[0] - c * t[1] + d * t[2]);
//   dest->m[1].w =  a * t[0] - c * t[3] + d * t[4];
//   dest->m[2].w =-(a * t[1] - b * t[3] + d * t[5]);
//   dest->m[3].w =  a * t[2] - b * t[4] + c * t[5];

//   det = 1.0f / (a * dest->m[0].x + b * dest->m[1].x
//               + c * dest->m[2].x + d * dest->m[3].x);

//   glm_mat4_scale_p(dest, det);
// #endif
// }

// /*!
//  * @brief inverse mat4 and store in dest
//  *
//  * this func uses reciprocal approximation without extra corrections
//  * e.g Newton-Raphson. this should work faster than normal,
//  * to get more precise use glm_mat4_inv version.
//  *
//  * NOTE: You will lose precision, glm_mat4_inv is more accurate
//  *
//  * @param[in]  mat  matrix
//  * @param[out] dest inverse matrix
//  */
// CGLM_INLINE
// void
// glm_mat4_inv_fast(mat4 mat, mat4 dest) {
// #if defined( __SSE__ ) || defined( __SSE2__ )
//   glm_mat4_inv_fast_sse2(mat, dest);
// #else
//   glm_mat4_inv(mat, dest);
// #endif
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
// glm_mat4_swap_col(mat4 mat, int col1, int col2) {
//   CGLM_ALIGN(16) vec4 tmp;
//   glm_vec4_copy(mat->m + col1, tmp);
//   glm_vec4_copy(mat->m + col2, mat->m + col1);
//   glm_vec4_copy(tmp, mat->m + col2);
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
// glm_mat4_swap_row(mat4 mat, int row1, int row2) {
//   CGLM_ALIGN(16) vec4 tmp;
//   (*tmp).x = ((float*)mat->m + 0)[row1];
//   (*tmp).y = ((float*)mat->m + 1)[row1];
//   (*tmp).z = ((float*)mat->m + 2)[row1];
//   (*tmp).w = ((float*)mat->m + 3)[row1];

//   ((float*)mat->m + 0)[row1] = ((float*)mat->m + 0)[row2];
//   ((float*)mat->m + 1)[row1] = ((float*)mat->m + 1)[row2];
//   ((float*)mat->m + 2)[row1] = ((float*)mat->m + 2)[row2];
//   ((float*)mat->m + 3)[row1] = ((float*)mat->m + 3)[row2];

//   ((float*)mat->m + 0)[row2] = (*tmp).x;
//   ((float*)mat->m + 1)[row2] = (*tmp).y;
//   ((float*)mat->m + 2)[row2] = (*tmp).z;
//   ((float*)mat->m + 3)[row2] = (*tmp).w;
// }

// /*!
//  * @brief helper for  R (row vector) * M (matrix) * C (column vector)
//  *
//  * rmc stands for Row * Matrix * Column
//  *
//  * the result is scalar because R * M = Matrix1x4 (row vector),
//  * then Matrix1x4 * Vec4 (column vector) = Matrix1x1 (Scalar)
//  *
//  * @param[in]  r   row vector or matrix1x4
//  * @param[in]  m   matrix4x4
//  * @param[in]  c   column vector or matrix4x1
//  *
//  * @return scalar value e.g. B(s)
//  */
// CGLM_INLINE
// float
// glm_mat4_rmc(vec4 r, mat4 m, vec4 c) {
//   vec4 tmp;
//   glm_mat4_mulv(m, c, tmp);
//   return glm_vec4_dot(r, tmp);
// }

#endif /* cglm_mat_h */
