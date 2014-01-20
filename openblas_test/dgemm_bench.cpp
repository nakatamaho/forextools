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

/* vanilla version of dgemm benchmark (using OpenBLAS) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#if defined USE_MACPORTS
#include <cblas_openblas.h>
#else
#include <cblas.h>
#endif
#include "mt19937ar.h"

double gettime(void)
{
    struct timeval mytime;
    gettimeofday(&mytime, (struct timezone *) 0);
    return (double) (mytime.tv_sec) + (double) (mytime.tv_usec) * 1.e-6;
}

static double getgflops(int k, double secs)
{
    double gflops;
    if (secs == 0.) return 0.;
    gflops = 1.e-9 * 2. * (double) k *(double) k *(double) k / secs;
    return gflops;
}

static void bench_gemm(int size, double *A, double *B, double *C)
{
    int i, iter = 1;
    double start, stop;
    double alpha = 1.;
    double beta = 1.;
    double gflops, peak;

    if ((double) size * (double) size * (double) size <= 1.e16)	iter = 3;
    if ((double) size * (double) size * (double) size <= 1.e12)	iter = 5;
    if ((double) size * (double) size * (double) size <= 1.e10)	iter = 10;
    if ((double) size * (double) size * (double) size <= 1.e9)	iter = 20;
    if ((double) size * (double) size * (double) size <= 1.e8)	iter = 40;

    peak = 0.;
    for (i = 0; i < iter; i++) {
	start = gettime();
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, size, size, size, alpha, A, size, B, size, beta, C, size);
	stop = gettime();
	gflops = getgflops(size, stop - start);
	if (peak < gflops) peak = gflops;
    }
    printf("%4d, %10.3f\n", size, peak);
    fflush(stdout);
    return;
}

int main(int argc, char *argv[])
{
    int size, j, p;
    init_genrand(0UL);
    for (p = 1024; p >= 1; p--) {
	size = p;
	double *A = new double[size * size];
	double *B = new double[size * size];
	double *C = new double[size * size];
	if (A == NULL || B == NULL || C == NULL) {
	    printf("Out of Memory!!\n");
	    exit(1);
	}

	for (j = 0; j < size * size; j++) A[j] = genrand_real1();
	for (j = 0; j < size * size; j++) B[j] = genrand_real1();
	for (j = 0; j < size * size; j++) C[j] = genrand_real1();

	bench_gemm(p, A, B, C);

	delete[]C; delete[]B; delete[]A;
    }
    return 0;
}
