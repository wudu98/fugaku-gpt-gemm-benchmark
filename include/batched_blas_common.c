#include "batched_blas_common.h"

// Group no. and local no. --> global no. 
int gl2g(const int group_no, const int local_no, const int *group_head){ 
  return group_head[group_no] + local_no; 
} 

CBLAS_LAYOUT layout_cblas(bblas_enum_t layout)
{
  if(layout == BblasRowMajor)  return CblasRowMajor;
  if(layout == BblasColMajor)  return CblasColMajor;
  return CblasRowMajor;
}
__CBLAS_trick__
CBLAS_TRANSPOSE transpose_cblas(bblas_enum_t transpose)
{
    if( transpose == BblasTrans )     return CblasTrans;
    if( transpose == BblasNoTrans )   return CblasNoTrans;
    if( transpose == BblasConjTrans ) return CblasConjTrans;
    return CblasNoTrans;
}
__CBLAS_trick__
CBLAS_UPLO uplo_cblas(bblas_enum_t uplo)
{
    if( uplo == BblasLower )  return CblasLower;
    if( uplo == BblasUpper )  return CblasUpper;
    return CblasLower;
}
__CBLAS_trick__
CBLAS_DIAG diag_cblas(bblas_enum_t diag)
{
    if( diag == BblasUnit )    return CblasUnit;
    if( diag == BblasNonUnit ) return CblasNonUnit;
    return CblasUnit;
}
__CBLAS_trick__
CBLAS_SIDE side_cblas(bblas_enum_t side)
{
    if( side == BblasRight ) return CblasRight;
    if( side == BblasLeft )  return CblasLeft;
    return CblasRight;
}
