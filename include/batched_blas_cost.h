#ifndef _BATCHED_BLAS_COST_H_
#define _BATCHED_BLAS_COST_H_

#include "batched_blas_common.h"

#define MAX_COSTPARAM 20

struct _cost_param{
	void *ptr[MAX_COSTPARAM];
};

double get_cost_const(struct _cost_param *cost_param, int group_no);
double get_cost_n3(struct _cost_param *cost_param, int group_no);
double get_cost_side_n3(struct _cost_param *cost_param, int group_no);
double get_cost_trans_n3(struct _cost_param *cost_param, int group_no);
double get_cost_n2(struct _cost_param *cost_param, int group_no);
double get_cost_n1(struct _cost_param *cost_param, int group_no);
double get_cost_band_trans_n2(struct _cost_param *cost_param, int group_no);
double get_cost_band_n2(struct _cost_param *cost_param, int group_no);

#endif
