#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#include "batched_blas_aprioricost.h"

#define typ float

static double getTime(){
  struct timeval t;
  gettimeofday(&t,NULL);
  return t.tv_sec + (double) t.tv_usec*1e-6;
}

int main(int argc, char *argv[])
{
  int i,j,l;
  int M, N, K,gsize,group_count,*group_size,*group_head;
  bblas_enum_t layout;
  bblas_enum_t *transa,*transb;
  
  int *m,*n,*k,*lda,*ldb,*ldc;
  typ *alpha,*beta;
  typ **a,**b;
  typ **c,**c_temp;
  double ts,te,gflops;
  double error,max_error;
  int total_batch_count;
  int a_size,b_size,c_size;
  
  if(argc < 6){
    printf("usage : ./(load module) <group_count> <group_size> <M> <N> <K>\n");
    exit(1);
  }

  group_count = atoi(argv[1]);
  gsize = atoi(argv[2]);
  M = atoi(argv[3]);
  N = atoi(argv[4]);
  K = atoi(argv[5]);
  
  // printf("name= %s, group_count= %d, group_size= %d, nsize= %d\n", argv[0], group_count, gsize, nsize);

  srand48(123L);

  group_size = (int *)malloc(sizeof(int) * group_count);
  group_head = (int *)malloc(sizeof(int) * group_count);
  transa = malloc(sizeof(bblas_enum_t) * group_count);
  transb = malloc(sizeof(bblas_enum_t) * group_count);
  m = (int *)malloc(sizeof(int) * group_count);
  n = (int *)malloc(sizeof(int) * group_count);
  k = (int *)malloc(sizeof(int) * group_count);
  lda = (int *)malloc(sizeof(int) * group_count);
  ldb = (int *)malloc(sizeof(int) * group_count);
  ldc = (int *)malloc(sizeof(int) * group_count);
  alpha = (typ *)malloc(sizeof(typ) * group_count);
  beta = (typ *)malloc(sizeof(typ) * group_count);

  // set group size
  total_batch_count = 0;
  for(i = 0; i < group_count; i++){
    group_size[i] = gsize;
    total_batch_count += group_size[i];
  }
  // set group head
  group_head[0] = 0;
  for(i = 1; i < group_count; i++){
    group_head[i] = group_head[i - 1] + group_size[i - 1];
  }

  a = (typ **)malloc(sizeof(typ *) * total_batch_count);
  b = (typ **)malloc(sizeof(typ *) * total_batch_count);
  c = (typ **)malloc(sizeof(typ *) * total_batch_count);
  c_temp = (typ **)malloc(sizeof(typ *) * total_batch_count);
  
  // set parameter and value
  // layout = BblasColMajor;
  layout = BblasRowMajor;
  for(i = 0; i < group_count; i++){
    transa[i] = BblasNoTrans;
    transb[i] = BblasNoTrans;
    m[i] = M;
    n[i] = N;
    k[i] = K;
    lda[i] = K;
    ldb[i] = N;
    ldc[i] = N;
    alpha[i] = 1;
    beta[i] = 0;
    for(j = 0; j < group_size[i]; j++){
      if((layout == BblasColMajor && transa[i] == BblasNoTrans) ||
	 (layout == BblasRowMajor && (transa[i] == BblasTrans || transa[i] == BblasConjTrans)))
	a_size = lda[i] * k[i];
      else
	a_size = lda[i] * m[i];
      if((layout == BblasColMajor && transb[i] == BblasNoTrans) ||
	 (layout == BblasRowMajor && (transb[i] == BblasTrans || transb[i] == BblasConjTrans)))
	b_size = ldb[i] * n[i];
      else
	b_size = ldb[i] * k[i];
      if(layout == BblasColMajor)
	c_size = ldc[i] * n[i];
      else
	c_size = ldc[i] * m[i];
      a[group_head[i]+j] = (typ *)malloc(sizeof(typ) * a_size);
      b[group_head[i]+j] = (typ *)malloc(sizeof(typ) * b_size);
      c[group_head[i]+j] = (typ *)malloc(sizeof(typ) * c_size);
      c_temp[group_head[i]+j] = (typ *)malloc(sizeof(typ) * c_size);
      for(l = 0; l < a_size; l++){
	a[group_head[i]+j][l] = drand48();
      }
      for(l = 0; l < b_size; l++){
	b[group_head[i]+j][l] = drand48();
      }
      for(l = 0; l < c_size; l++){
	c_temp[group_head[i]+j][l] = c[group_head[i]+j][l] = drand48();
      }
    }
  }
  
  typ *mem_refresh = (typ *)malloc(sizeof(typ) * 4096 * 4096);

  int info_size   = 1;
  int info_option = BblasErrorsReportNone;
  int *info = (int*) malloc((size_t)info_size*sizeof(int))  ;
  info[0] = info_option;
  
  // for(l = 0; l < 4096 * 4096; l++)
	//   mem_refresh[l] = drand48();

  // // batched
  // ts = getTime();
  // blas_sgemm_batch_aprioricost(group_count, group_size, layout, transa, transb, m, n, k, alpha, (const typ **)a, lda, (const typ **)b, ldb, beta, (typ **)c, ldc, info);
  // te = getTime();
  // gflops = 2.0 * total_batch_count * M * N * K / (te - ts) / 1000000000;
  // printf("batched= %.6f s, %.2f gflops, ", te - ts, gflops);

  // for(l = 0; l < 4096 * 4096; l++)
	//   mem_refresh[l] = drand48();

  // // my_batched
  // ts = getTime();
  // blas_my_sgemm_batch_aprioricost(group_count, group_size, layout, transa, transb, m, n, k, alpha, (const typ **)a, lda, (const typ **)b, ldb, beta, (typ **)c, ldc, info);
  // te = getTime();
  // gflops = 2.0 * total_batch_count * M * N * K / (te - ts) / 1000000000;
  // printf("my_batched= %.6f s, %.2f gflops, ", te - ts, gflops);

  for(l = 0; l < 4096 * 4096; l++)
	  mem_refresh[l] = drand48();

  // non-batched
  ts = getTime();
  for(i = 0; i < group_count; i++){
    for(j = 0; j < group_size[i]; j++){
      cblas_sgemm(layout_cblas( layout ) , transpose_cblas( transa[i] ) , transpose_cblas( transb[i] ) ,
		  m[i], n[i], k[i], alpha[i], a[group_head[i]+j], lda[i],
		  b[group_head[i]+j], ldb[i], beta[i], c_temp[group_head[i]+j], ldc[i]);
    }
  }
  te = getTime();
  gflops = 2.0 * total_batch_count * M * N * K / (te - ts) / 1000000000;
  printf("no_batched= %.6f s, %.2f gflops, ", te - ts, gflops);

  // // check
  // max_error = 0.0;

  // for(i = 0; i < group_count; i++){
  //   if(layout == BblasColMajor)
  //     c_size = ldc[i] * n[i];
  //   else
  //     c_size = ldc[i] * m[i];
  //   for(j = 0; j < group_size[i]; j++){
  //     for(l = 0; l < c_size; l++){
	// error = fabs(c_temp[group_head[i]+j][l]-c[group_head[i]+j][l]);
  // // printf("%f,  %f\n", c[group_head[i]+j][l], c_temp[group_head[i]+j][l]);
	// if(error > max_error)
	//   max_error = error;
  //     }
  //   }
  // }

  // printf("max_error= %e", max_error);
  free(info);
  free(mem_refresh);
  printf("\n");

  for(i = 0; i < group_count; i++){
    for(j = 0; j < group_size[i]; j++){
      free(a[group_head[i]+j]);
      free(b[group_head[i]+j]);
      free(c[group_head[i]+j]);
      free(c_temp[group_head[i]+j]);
    }
  }
  free(a);
  free(b);
  free(c);
  free(c_temp);
  
  free(group_size);
  free(group_head);
  free(transa);
  free(transb);
  free(m);
  free(n);
  free(k);
  free(lda);
  free(ldb);
  free(ldc);
  free(alpha);
  free(beta);
  
}

