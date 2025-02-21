/*
 * Copyright (c), Recep Aslantas.
 *
 * MIT License (MIT), http://opensource.org/licenses/MIT
 * Full license can be found in the LICENSE file
 */

#ifndef cglm_types_h
#define cglm_types_h

#define CGLM_ALL_UNALIGNED
#define CGLM_ALIGN(X) /* no alignment */
#define CGLM_ALIGN_IF(X) /* no alignment */
#define CGLM_ALIGN_MAT

typedef float3*					   vec3;
typedef int3                      ivec3;
typedef float4*					   vec4;
typedef vec4                    versor;

// mat3 is essentially a mat4 under the hood!
typedef struct _mat3 { float4 m[4]; } _mat3;
typedef _mat3*						  mat3;

// _mat4 to be used as parameters to kernel functions.
// cglm types are essentially pointer types ( float[4] is equivalent to float*, and thus can't be
// passed as kernel argument )
typedef struct _mat4 { float4 m[4]; } _mat4;

// mat4 pointer type to be used as parameters to cglm functions
typedef _mat4*						mat4;

#define GLM_E         2.71828182845904523536028747135266250   /* e           */
#define GLM_LOG2E     1.44269504088896340735992468100189214   /* log2(e)     */
#define GLM_LOG10E    0.434294481903251827651128918916605082  /* log10(e)    */
#define GLM_LN2       0.693147180559945309417232121458176568  /* loge(2)     */
#define GLM_LN10      2.30258509299404568401799145468436421   /* loge(10)    */
#define GLM_PI        3.14159265358979323846264338327950288   /* pi          */
#define GLM_PI_2      1.57079632679489661923132169163975144   /* pi/2        */
#define GLM_PI_4      0.785398163397448309615660845819875721  /* pi/4        */
#define GLM_1_PI      0.318309886183790671537767526745028724  /* 1/pi        */
#define GLM_2_PI      0.636619772367581343075535053490057448  /* 2/pi        */
#define GLM_2_SQRTPI  1.12837916709551257389615890312154517   /* 2/sqrt(pi)  */
#define GLM_SQRT2     1.41421356237309504880168872420969808   /* sqrt(2)     */
#define GLM_SQRT1_2   0.707106781186547524400844362104849039  /* 1/sqrt(2)   */

#define GLM_Ef        ((float)GLM_E)
#define GLM_LOG2Ef    ((float)GLM_LOG2E)
#define GLM_LOG10Ef   ((float)GLM_LOG10E)
#define GLM_LN2f      ((float)GLM_LN2)
#define GLM_LN10f     ((float)GLM_LN10)
#define GLM_PIf       ((float)GLM_PI)
#define GLM_PI_2f     ((float)GLM_PI_2)
#define GLM_PI_4f     ((float)GLM_PI_4)
#define GLM_1_PIf     ((float)GLM_1_PI)
#define GLM_2_PIf     ((float)GLM_2_PI)
#define GLM_2_SQRTPIf ((float)GLM_2_SQRTPI)
#define GLM_SQRT2f    ((float)GLM_SQRT2)
#define GLM_SQRT1_2f  ((float)GLM_SQRT1_2)

/* DEPRECATED! use GLM_PI and friends */
#define CGLM_PI       GLM_PIf
#define CGLM_PI_2     GLM_PI_2f
#define CGLM_PI_4     GLM_PI_4f

#endif /* cglm_types_h */
