/*
Copyright 2014 NAKATA Maho. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.
Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY NAKATA Maho ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
EVENT SHALL NAKATA Maho OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of NAKATA Maho.
*/

/* vanilla version of dsyev benchmark (using OpenBLAS) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#if defined USE_MKL
#include <mkl_lapacke.h>
#else 
#include <lapacke.h>
#endif
#include "mt19937ar.h"

double gettime(void)
{
    struct timeval mytime;
    gettimeofday(&mytime, (struct timezone *) 0);
    return (double) (mytime.tv_sec) + (double) (mytime.tv_usec) * 1.e-6;
}

static void bench_syev(int size, double *A, double *w)
{
    int i, ret, iter = 1;
    double start, stop;
    double wall, peak;

    if ((double) size * (double) size * (double) size <= 1.e12)	iter = 2;
    if ((double) size * (double) size * (double) size <= 1.e10)	iter = 3;
    if ((double) size * (double) size * (double) size <= 1.e8)	iter = 5;

    peak = 0.;
    for (i = 0; i < iter; i++) {
	start = gettime();
        ret = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', size, A, size, w);
	stop = gettime();
	wall = stop - start;
	if (peak < wall) peak = wall;
    }
    printf("%4d, %10.3f\n", size, peak);
    fflush(stdout);
    return;
}

int main(int argc, char *argv[])
{
    int size, j, p;
    init_genrand(0UL);
    for (p = 512; p >= 500; p--) {
	size = p;
	double *A = new double[size * size];
	double *w = new double[size];
	if (A == NULL || w == NULL) {
	    printf("Out of Memory!!\n");
	    exit(1);
	}
	for (j = 0; j < size * size; j++) A[j] = genrand_real1();
	for (j = 0; j < size ; j++) w[j] = genrand_real1();
	bench_syev(size, A, w);
	delete[]w;
	delete[]A;
    }
    return 0;
}
