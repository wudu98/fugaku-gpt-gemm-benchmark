ifeq (x$(CC),xcc)
ifdef MKLROOT
CC=icc
else
ifdef FJSVXTCLANGA
CC=fccpx
else
CC=gcc
endif
endif
endif

ifeq (x$(CC),xicc)
CCFLAG=-O3 -xHost -qopenmp -I./
BLAS=-lmkl_rt -lm
else
ifeq (x$(CC),xfccpx)
CCFLAG=-std=gnu11 -Kfast,ocl,openmp -Nclang -D_CBLAS_ -I./
BLAS=-SSL2BLAMP -lm
else
ifeq (x$(CC),xarmclang)
CCFLAG=-O3 -fopenmp -Wall -armpl=sve -mcpu=a64fx -D_CBLAS_ -I./
BLAS=-armpl -lm
else
CCFLAG=-O3 -fopenmp -Wall -D_CBLAS_ -I./ -I/usr/include/openblas
BLAS=-lopenblas -lm
endif
endif
endif

DEF = -I../bblas_src_aprioricost
LIB = -Wl,-rpath,../bblas_src_aprioricost -L../bblas_src_aprioricost -lbblas_aprioricost -lm -lfjomp

SRC=batched_blas_sgemm_sample.c
EXE=batched_blas_sgemm_sample

benchmark : 
	$(CC) $(SRC) $(DEF) $(CCFLAG) $(LIB) $(BLAS) -o $(EXE)

clean :
	rm -f *.o *~ $(OUT)


