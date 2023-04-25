#ifndef _BATCHED_BLAS_FP16_H_
#define _BATCHED_BLAS_FP16_H_

#ifdef __aarch64__
#  if !defined(__FUJITSU) && !defined(__CLANG_FUJITSU)
//#    define FP16_NATIVE_SUPPORT
#    define FP16_AUTO_PROMOTION
#  else
#    define FP16_FUJITSU_TRAD_MODE
#  endif
#elif defined(BF_NMANT)
#  define FP16_BFLIKE_FLOAT
#  if BF_NMANT>7 || BF_NMANT<=1
#    error "too large or small mantissa for BFLIKE_FLOAT"
#  endif
#elif defined(__AVX2__)
#  define FP16_AVX2_EMULATION
#elif defined(__clang__) && __clang_major__ >= 8
#  define FP16_AUTO_PROMOTION
#else
#  define FP16_IS_NOT_SUPPORTED
#endif


#ifdef FP16_NATIVE_SUPPORT
typedef _Float16 fp16;
#endif

#ifdef FP16_AUTO_PROMOTION
typedef __fp16 fp16;
#endif

#ifdef FP16_FUJITSU_TRAD_MODE
typedef __fp16 fp16;
#endif

#ifdef FP16_NATIVE_SUPPORT
typedef _Float16 fp16;
#endif

#ifdef FP16_IS_NOT_SUPPORTED
//#warning "FP16 IS NOT SUPPORTED"
typedef unsigned short fp16;
#endif

#endif

