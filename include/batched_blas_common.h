#ifndef _BATCHED_BLAS_COMMON_H_
#define _BATCHED_BLAS_COMMON_H_

#include <bblas_types.h>
#include <bblas_error.h>

#ifdef _CBLAS_
#include <cblas.h>
#ifdef _ARMPL_
#include <armpl.h>
#endif
#else
#include <mkl.h>
#endif

int gl2g(const int group_no, const int local_no, const int *group_head); 

#include <batched_blas_fp16.h>

#define __CBLAS_trick__

#ifdef _CBLAS_
  #if defined(__FUJITSU) || defined(__CLANG_FUJITSU)
    // Fujitsu SSL-II
  #elif defined(_ARMPL_) || defined(CBLAS_ENUM_DEFINED_H)
    // ARM Performance Library or ATLAS
    #undef __CBLAS_trick__
    #define __CBLAS_trick__	enum
    typedef __CBLAS_trick__ CBLAS_ORDER CBLAS_LAYOUT;
  #elif defined(OPENBLAS_CONST)
    // for compatibility cblas.h in OpenBLAS
    typedef CBLAS_ORDER CBLAS_LAYOUT;
  #endif
#endif

// change from bblas parameter to cblas parameter
 CBLAS_LAYOUT layout_cblas(bblas_enum_t layout);
__CBLAS_trick__
 CBLAS_TRANSPOSE transpose_cblas(bblas_enum_t transpose);
__CBLAS_trick__
 CBLAS_UPLO uplo_cblas(bblas_enum_t uplo);
__CBLAS_trick__
 CBLAS_DIAG diag_cblas(bblas_enum_t uplo);
__CBLAS_trick__
 CBLAS_SIDE side_cblas(bblas_enum_t uplo);

#endif
