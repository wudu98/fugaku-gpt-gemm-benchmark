# fugaku-gpt-gemm-benchmark

Batched BLAS Generator 1.2
Copyright (C) 2020-2022 RIKEN.

This is a micro benchmark with the representative GEMM sizes in Fugaku GPT Project.

The interface of this project is BMM in Batched BLAS Generator (https://www.r-ccs.riken.jp/labs/lpnctrt/projects/batchedblas/index.html).

**Build dependencies:**
* Fujitsu cross-compiler 
* SSL-II

**Compilation**

```
 # run on Login node
 bash ./make_ssl2.sh 
```

**Benchmark**

```
 # run on Computing node
 bash ./run_benchmark.sh [THREADS]
```

Results in (https://docs.google.com/spreadsheets/d/1JenD_QRAMcWaumgb4V4j8zzijcSz8wgZw2MmX_3T4Sg/edit#gid=0)
