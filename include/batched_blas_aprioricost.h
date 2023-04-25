#ifndef _BATCHED_BLAS_APRIORICOST_H_
#define _BATCHED_BLAS_APRIORICOST_H_

#include "batched_blas_common.h"

void blas_sgemm_batchf_aprioricost(const int group_size,const bblas_enum_t layout,const bblas_enum_t transa,const bblas_enum_t transb,const int m,const int n,const int k,const float alpha,const float ** a,const int lda,const float ** b,const int ldb,const float beta,float ** c,const int ldc,int *info);
void blas_my_sgemm_batchf_aprioricost(const int group_size,const bblas_enum_t layout,const bblas_enum_t transa,const bblas_enum_t transb,const int m,const int n,const int k,const float alpha,const float ** a,const int lda,const float ** b,const int ldb,const float beta,float ** c,const int ldc,int *info);
void blas_sgemm_batch_aprioricost(const int group_count,const int *group_size,const bblas_enum_t layout,const bblas_enum_t* transa,const bblas_enum_t* transb,const int* m,const int* n,const int* k,const float* alpha,const float ** a,const int* lda,const float ** b,const int* ldb,const float* beta,float ** c,const int* ldc,int *info);
void blas_my_sgemm_batch_aprioricost(const int group_count,const int *group_size,const bblas_enum_t layout,const bblas_enum_t* transa,const bblas_enum_t* transb,const int* m,const int* n,const int* k,const float* alpha,const float ** a,const int* lda,const float ** b,const int* ldb,const float* beta,float ** c,const int* ldc,int *info);

#endif
