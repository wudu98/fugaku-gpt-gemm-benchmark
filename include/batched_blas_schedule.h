#ifndef _BATCHED_BLAS_SCHEDULE_H_
#define _BATCHED_BLAS_SCHEDULE_H_

#include "batched_blas_common.h"
#include "batched_blas_cost.h"

void schedule_batch(const int group_count, const int *group_size,
		    double (*cost_func)(struct _cost_param *cost_param, int group_no), struct _cost_param *cost_param, int *which_thread);
void schedule_batch_inorder (const int group_count, const int *group_size, int *which_thread);
int use_batch(void);

#endif
