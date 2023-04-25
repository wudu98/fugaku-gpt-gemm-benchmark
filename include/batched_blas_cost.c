#include "batched_blas_cost.h" 

double get_cost_const(struct _cost_param *cost_param, int group_no)
{
  return 1.0;
}

double get_cost_n3(struct _cost_param *cost_param, int group_no)
{
  int n1,n2,n3;
  
  n1 = ((int *)cost_param->ptr[0])[group_no];
  n2 = ((int *)cost_param->ptr[1])[group_no];
  n3 = ((int *)cost_param->ptr[2])[group_no];
  
  return (double)n1 * (double)n2 * (double)n3;
}

double get_cost_side_n3(struct _cost_param *cost_param, int group_no)
{
  int n1,n2;
__CBLAS_trick__
  CBLAS_SIDE side;
  
  side = ((__CBLAS_trick__ CBLAS_SIDE *)cost_param->ptr[0])[group_no];
  n1 = ((int *)cost_param->ptr[1])[group_no];
  n2 = ((int *)cost_param->ptr[2])[group_no];

  if(side == CblasLeft)
    return (double)n1 * (double)n1 * (double)n2;
  else
    return (double)n1 * (double)n2 * (double)n2;
    
}

double get_cost_trans_n3(struct _cost_param *cost_param, int group_no)
{
  int n1,n2;
__CBLAS_trick__
  CBLAS_TRANSPOSE trans;
  
  trans = (( __CBLAS_trick__ CBLAS_TRANSPOSE *)cost_param->ptr[0])[group_no];
  n1 = ((int *)cost_param->ptr[1])[group_no];
  n2 = ((int *)cost_param->ptr[2])[group_no];

  if(trans == CblasNoTrans)
    return (double)n1 * (double)n1 * (double)n2;
  else
    return (double)n1 * (double)n2 * (double)n2;
    
}

double get_cost_n2(struct _cost_param *cost_param, int group_no)
{
  int n1,n2;
  
  n1 = ((int *)cost_param->ptr[0])[group_no];
  n2 = ((int *)cost_param->ptr[1])[group_no];
  
  return (double)n1 * (double)n2;
}

double get_cost_band_trans_n2(struct _cost_param *cost_param, int group_no)
{
  int n1,n2;
  int kl,ku;
__CBLAS_trick__
  CBLAS_TRANSPOSE trans;
  
  trans = (( __CBLAS_trick__ CBLAS_TRANSPOSE *)cost_param->ptr[0])[group_no];
  n1 = ((int *)cost_param->ptr[1])[group_no];
  n2 = ((int *)cost_param->ptr[2])[group_no];
  kl = ((int *)cost_param->ptr[3])[group_no];
  ku = ((int *)cost_param->ptr[4])[group_no];

  if(trans == CblasNoTrans)
    return ((double)kl + (double)ku + 1.0) * (double)n2;
  else
    return (double)n1 * ((double)kl + (double)ku + 1.0);
    
}

double get_cost_band_n2(struct _cost_param *cost_param, int group_no)
{
  int n1;
  int kl,ku;
  
  n1 = ((int *)cost_param->ptr[0])[group_no];
  kl = ((int *)cost_param->ptr[1])[group_no];
  ku = ((int *)cost_param->ptr[2])[group_no];

  return (double)n1 * ((double)kl + (double)ku + 1.0);
    
}

double get_cost_n1(struct _cost_param *cost_param, int group_no)
{
  int n1;
  
  n1 = ((int *)cost_param->ptr[0])[group_no];
  
  return (double)n1;
}

