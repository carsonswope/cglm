/*
 * Copyright (c), Recep Aslantas.
 *
 * MIT License (MIT), http://opensource.org/licenses/MIT
 * Full license can be found in the LICENSE file
 */

/*
 Functions:
   CGLM_INLINE void glm_mul(mat4 m1, mat4 m2, mat4 dest);
   CGLM_INLINE void glm_inv_tr(mat4 mat);
 */

#ifndef cglm_affine_mat_h
#define cglm_affine_mat_h

#include "common.h"
#include "mat4.h"
#include "mat3.h"

#ifdef CGLM_SSE_FP
#  include "simd/sse2/affine.h"
#endif

#ifdef CGLM_AVX_FP
#  include "simd/avx/affine.h"
#endif

/*!
 * @brief this is similar to glm_mat4_mul but specialized to affine transform
 *
 * Matrix format should be:
 *   R  R  R  X
 *   R  R  R  Y
 *   R  R  R  Z
 *   0  0  0  W
 *
 * this reduces some multiplications. It should be faster than mat4_mul.
 * if you are not sure about matrix format then DON'T use this! use mat4_mul
 *
 * @param[in]   m1    affine matrix 1
 * @param[in]   m2    affine matrix 2
 * @param[out]  dest  result matrix
 */
// CGLM_INLINE
// void
// glm_mul(mat4 m1, mat4 m2, mat4 dest) {
// #ifdef __AVX__
//   glm_mul_avx(m1, m2, dest);
// #elif defined( __SSE__ ) || defined( __SSE2__ )
//   glm_mul_sse2(m1, m2, dest);
// #else
//   float a00 = m1->m[0].x, a01 = m1->m[0].y, a02 = m1->m[0].z, a03 = m1->m[0].w,
//         a10 = m1->m[1].x, a11 = m1->m[1].y, a12 = m1->m[1].z, a13 = m1->m[1].w,
//         a20 = m1->m[2].x, a21 = m1->m[2].y, a22 = m1->m[2].z, a23 = m1->m[2].w,
//         a30 = m1->m[3].x, a31 = m1->m[3].y, a32 = m1->m[3].z, a33 = m1->m[3].w,

//         b00 = m2->m[0].x, b01 = m2->m[0].y, b02 = m2->m[0].z,
//         b10 = m2->m[1].x, b11 = m2->m[1].y, b12 = m2->m[1].z,
//         b20 = m2->m[2].x, b21 = m2->m[2].y, b22 = m2->m[2].z,
//         b30 = m2->m[3].x, b31 = m2->m[3].y, b32 = m2->m[3].z, b33 = m2->m[3].w;

//   dest->m[0].x = a00 * b00 + a10 * b01 + a20 * b02;
//   dest->m[0].y = a01 * b00 + a11 * b01 + a21 * b02;
//   dest->m[0].z = a02 * b00 + a12 * b01 + a22 * b02;
//   dest->m[0].w = a03 * b00 + a13 * b01 + a23 * b02;

//   dest->m[1].x = a00 * b10 + a10 * b11 + a20 * b12;
//   dest->m[1].y = a01 * b10 + a11 * b11 + a21 * b12;
//   dest->m[1].z = a02 * b10 + a12 * b11 + a22 * b12;
//   dest->m[1].w = a03 * b10 + a13 * b11 + a23 * b12;

//   dest->m[2].x = a00 * b20 + a10 * b21 + a20 * b22;
//   dest->m[2].y = a01 * b20 + a11 * b21 + a21 * b22;
//   dest->m[2].z = a02 * b20 + a12 * b21 + a22 * b22;
//   dest->m[2].w = a03 * b20 + a13 * b21 + a23 * b22;

//   dest->m[3].x = a00 * b30 + a10 * b31 + a20 * b32 + a30 * b33;
//   dest->m[3].y = a01 * b30 + a11 * b31 + a21 * b32 + a31 * b33;
//   dest->m[3].z = a02 * b30 + a12 * b31 + a22 * b32 + a32 * b33;
//   dest->m[3].w = a03 * b30 + a13 * b31 + a23 * b32 + a33 * b33;
// #endif
// }

// !
//  * @brief this is similar to glm_mat4_mul but specialized to affine transform
//  *
//  * Right Matrix format should be:
//  *   R  R  R  0
//  *   R  R  R  0
//  *   R  R  R  0
//  *   0  0  0  1
//  *
//  * this reduces some multiplications. It should be faster than mat4_mul.
//  * if you are not sure about matrix format then DON'T use this! use mat4_mul
//  *
//  * @param[in]   m1    affine matrix 1
//  * @param[in]   m2    affine matrix 2
//  * @param[out]  dest  result matrix
 
// CGLM_INLINE
// void
// glm_mul_rot(mat4 m1, mat4 m2, mat4 dest) {
// #if defined( __SSE__ ) || defined( __SSE2__ )
//   glm_mul_rot_sse2(m1, m2, dest);
// #else
//   float a00 = m1->m[0].x, a01 = m1->m[0].y, a02 = m1->m[0].z, a03 = m1->m[0].w,
//         a10 = m1->m[1].x, a11 = m1->m[1].y, a12 = m1->m[1].z, a13 = m1->m[1].w,
//         a20 = m1->m[2].x, a21 = m1->m[2].y, a22 = m1->m[2].z, a23 = m1->m[2].w,
//         a30 = m1->m[3].x, a31 = m1->m[3].y, a32 = m1->m[3].z, a33 = m1->m[3].w,

//         b00 = m2->m[0].x, b01 = m2->m[0].y, b02 = m2->m[0].z,
//         b10 = m2->m[1].x, b11 = m2->m[1].y, b12 = m2->m[1].z,
//         b20 = m2->m[2].x, b21 = m2->m[2].y, b22 = m2->m[2].z;

//   dest->m[0].x = a00 * b00 + a10 * b01 + a20 * b02;
//   dest->m[0].y = a01 * b00 + a11 * b01 + a21 * b02;
//   dest->m[0].z = a02 * b00 + a12 * b01 + a22 * b02;
//   dest->m[0].w = a03 * b00 + a13 * b01 + a23 * b02;

//   dest->m[1].x = a00 * b10 + a10 * b11 + a20 * b12;
//   dest->m[1].y = a01 * b10 + a11 * b11 + a21 * b12;
//   dest->m[1].z = a02 * b10 + a12 * b11 + a22 * b12;
//   dest->m[1].w = a03 * b10 + a13 * b11 + a23 * b12;

//   dest->m[2].x = a00 * b20 + a10 * b21 + a20 * b22;
//   dest->m[2].y = a01 * b20 + a11 * b21 + a21 * b22;
//   dest->m[2].z = a02 * b20 + a12 * b21 + a22 * b22;
//   dest->m[2].w = a03 * b20 + a13 * b21 + a23 * b22;

//   dest->m[3].x = a30;
//   dest->m[3].y = a31;
//   dest->m[3].z = a32;
//   dest->m[3].w = a33;
// #endif
// }

// /*!
//  * @brief inverse orthonormal rotation + translation matrix (ridig-body)
//  *
//  * @code
//  * X = | R  T |   X' = | R' -R'T |
//  *     | 0  1 |        | 0     1 |
//  * @endcode
//  *
//  * @param[in,out]  mat  matrix
//  */
// CGLM_INLINE
// void
// glm_inv_tr(mat4 mat) {
// #if defined( __SSE__ ) || defined( __SSE2__ )
//   glm_inv_tr_sse2(mat);
// #else
//   CGLM_ALIGN_MAT _mat3 r;
//   CGLM_ALIGN(8)  float3 t;

//   /* rotate */
//   glm_mat4_pick3t(mat, &r);
//   glm_mat4_ins3(&r, mat);

//   /* translate */
//   glm_mat3_mulv(&r, (vec3)(mat->m + 3), &t);
//   glm_vec3_negate(&t);
//   glm_vec3_copy(&t, (vec3)(mat->m + 3));
// #endif
// }

#endif /* cglm_affine_mat_h */
