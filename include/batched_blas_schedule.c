#include <stdlib.h>
#include <omp.h>
#include "batched_blas_common.h"
#include "batched_blas_cost.h"
#include "batched_blas_schedule.h"

static int compare_cost(const void *a, const void *b, void *cost_by_group) { 
  return ((double *)cost_by_group)[*(int *)b] - ((double *)cost_by_group)[*(int *)a]; 
} 

static int pick_min_thread(const double *cost, const int n)
{ 
  int i; 
  int current_min_i = 0; 
  for (i = 1; i < n; ++i) { 
    if (cost[i] < cost[current_min_i]) { 
      current_min_i = i; 
    } 
  } 
  return current_min_i; 
} 

void schedule_batch_inorder (const int group_count, const int *group_size, int *which_thread)
{
  const int num_threads = omp_get_max_threads();
  int local_no, group_no, tid;
  int group_head[group_count]; // The global no. of i-th group's first task is group_head[i]. 

  // set head index of the group
  group_head[0] = 0; 
  for (group_no = 1; group_no < group_count; group_no++) { 
    group_head[group_no] = group_head[group_no - 1] + group_size[group_no-1]; 
  } 

  // Task allocation 
  int threadno = 0;
  for (group_no = 0; group_no < group_count; group_no++) { 
    for (local_no = 0; local_no < group_size[group_no]; local_no++) { 
      which_thread[gl2g(group_no, local_no, group_head)] = threadno % num_threads;
	  threadno++;
    } 
  } 
}
 
void schedule_batch(const int group_count, const int *group_size,
   double (*cost_func)(struct _cost_param *cost_param, int group_no), struct _cost_param *cost_param, int *which_thread)
{
  int total_batch_count = 0; 
  const int num_threads = omp_get_max_threads();
  int local_no, group_no, tid;
  double cost_by_group[group_count]; 
  int cost_rank[group_count];
  int group_head[group_count]; // The global no. of i-th group's first task is group_head[i]. 
  double current_cost[num_threads]; // current_cost[i] = j : Current cost of the i-th thread is j. 

  // set total batch count
  for (group_no = 0; group_no < group_count; group_no++) { 
    total_batch_count += group_size[group_no];
  }

  // set head index of the group
  group_head[0] = 0; 
  for (group_no = 1; group_no < group_count; group_no++) { 
    group_head[group_no] = group_head[group_no - 1] + group_size[group_no-1]; 
  } 

  // Greedy task allocation 
  for (tid = 0; tid < num_threads; tid++) { 
    current_cost[tid] = 0.0; 
  } 

  // compute cost of each group 
  for (group_no = 0; group_no < group_count; group_no++) { 
    cost_by_group[group_no] = cost_func(cost_param,group_no); 
  }

  // Sort by cost 
  for (group_no = 0; group_no < group_count; group_no++) { 
    cost_rank[group_no] = group_no; 
  }
  qsort_r(cost_rank, group_no, sizeof(int), compare_cost, cost_by_group);

  // Task allocation 
  for (group_no = 0; group_no < group_count; group_no++) { 
    for (local_no = 0; local_no < group_size[cost_rank[group_no]]; local_no++) { 
      const int min_thread = pick_min_thread(current_cost, num_threads); 
      which_thread[gl2g(cost_rank[group_no], local_no, group_head)] = min_thread;
      current_cost[which_thread[gl2g(cost_rank[group_no], local_no, group_head)]] += cost_by_group[cost_rank[group_no]]; 
    } 
  } 

}

// condition for using batch
//  return 1:use , 0:not use
int use_batch()
{
  const int num_threads = omp_get_max_threads();

  if(num_threads > 1)
    return 1;
  else
    return 0;
}

