#!/bin/bash
set -e

tmp=`dirname $0`
PROJECT_ROOT=`cd $tmp; pwd`
cd ${PROJECT_ROOT}

export CC=fccpx
export BASEBLAS=SSL2

which python >& /dev/null
if [ $? -eq 0 ]; then
	PYTHON=python
fi
which python2 >& /dev/null
if [ $? -eq 0 ]; then
	PYTHON=python2
fi
which python3 >& /dev/null
if [ $? -eq 0 ]; then
	PYTHON=python3
fi
if [ x$PYTHON == x ]; then
	echo 'Can not find python command'
	exit 1
fi

make clean

$PYTHON data/bblas_data_gen.py data/bblas_data_ssl2.csv aprioricost
$PYTHON bblas_aprioricost.py data/bblas_data_ssl2_aprioricost.csv
make aprioricost 
cp bblas_src_aprioricost/*.a lib
cp bblas_src_aprioricost/*.so lib
cp bblas_src_aprioricost/batched_blas_aprioricost.h include

cd benchmark
make
