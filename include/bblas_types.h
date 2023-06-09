#ifndef BBLAS_TYPES_H
#define BBLAS_TYPES_H

#include <complex.h>
#ifndef CBLAS_SADDR
#define CBLAS_SADDR(var) &(var)
#endif

enum {
  BblasInvalid       = -1,
  
  BblasRowMajor      = 101,
  BblasColMajor      = 102,

  BblasNoTrans       = 111,
  BblasTrans         = 112,
  BblasConjTrans     = 113,
  Bblas_ConjTrans    = BblasConjTrans,

  BblasUpper         = 121,
  BblasLower         = 122,
  BblasGeneral       = 123,
  BblasGeneralBand   = 124,

  BblasNonUnit       = 131,
  BblasUnit          = 132,

  BblasLeft          = 141,
  BblasRight         = 142,

  BblasOneNorm       = 171,
  BblasRealOneNorm   = 172,
  BblasTwoNorm       = 173,
  BblasFrobeniusNorm = 174,
  BblasInfNorm       = 175,
  BblasRealInfNorm   = 176,
  BblasMaxNorm       = 177,
  BblasRealMaxNorm   = 178,

  BblasForward       = 391,
  BblasBackward      = 392,

  BblasColumnwise    = 401,
  BblasRowwise       = 402,

  BblasW             = 501,
  BblasA2            = 502

};


enum {
  BblasSuccess = 0,
  BblasFail
};

enum {
  BblasErrorsReportAll = 0,
  BblasErrorsReportGroup,
  BblasErrorsReportAny,
  BblasErrorsReportNone 
};

/******************************************************************************/
typedef int bblas_enum_t;
typedef float  _Complex bblas_complex32_t;
typedef double _Complex bblas_complex64_t;

/******************************************************************************/
bblas_enum_t bblas_diag_const(char lapack_char);
bblas_enum_t bblas_direct_const(char lapack_char);
bblas_enum_t bblas_norm_const(char lapack_char);
bblas_enum_t bblas_side_const(char lapack_char);
bblas_enum_t bblas_storev_const(char lapack_char);
bblas_enum_t bblas_trans_const(char lapack_char);
bblas_enum_t bblas_uplo_const(char lapack_char);
bblas_enum_t bblas_info_const(char lapack_char);

/******************************************************************************/
static inline int imin(int a, int b)
{
  if (a < b){
    return a;
  }else{
    return b;
  }
}

/******************************************************************************/
static inline int imax(int a, int b)
{
  if (a > b){
    return a;
  }else{
    return b;
  }
}

#endif  // BBLAS_TYPES_H 
