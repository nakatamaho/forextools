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

#include <stdio.h>
#if defined USE_MACPORTS
#include <cblas_openblas.h>
#elif defined USE_MKL
#include <mkl_cblas.h>
#define blasint MKL_INT
#else 
#include <cblas.h>
#endif

//Matlab/Octave format 
void printmat(blasint N, blasint M, double *A, blasint LDA)
{
    double mtmp;
    printf("[ ");
    for (blasint i = 0; i < N; i++) {
	printf("[ ");
	for (blasint j = 0; j < M; j++) {
	    mtmp = A[i + j * LDA];
	    printf("%5.2e", mtmp);
	    if (j < M - 1)
		printf(", ");
	}
	if (i < N - 1)    printf("]; ");
	else    printf("] ");
    }
    printf("]");
}

int main()
{
    blasint n = 3;
    double alpha, beta;
    double *A = new double[n * n];
    double *B = new double[n * n];
    double *C = new double[n * n];

    A[0 + 0 * n] = 1;    A[0 + 1 * n] = 8;     A[0 + 2 * n] = 3;
    A[1 + 0 * n] = 2;    A[1 + 1 * n] = 10;    A[1 + 2 * n] = 8;
    A[2 + 0 * n] = 9;    A[2 + 1 * n] = -5;    A[2 + 2 * n] = -1;

    B[0 + 0 * n] = 9;    B[0 + 1 * n] = 8;     B[0 + 2 * n] = 3;
    B[1 + 0 * n] = 3;    B[1 + 1 * n] = 11;    B[1 + 2 * n] = 2.3;
    B[2 + 0 * n] = -8;   B[2 + 1 * n] = 6;     B[2 + 2 * n] = 1;

    C[0 + 0 * n] = 3;    C[0 + 1 * n] = 3;     C[0 + 2 * n] = 1.2;
    C[1 + 0 * n] = 8;    C[1 + 1 * n] = 4;     C[1 + 2 * n] = 8;
    C[2 + 0 * n] = 6;    C[2 + 1 * n] = 1;     C[2 + 2 * n] = -2;

    printf("# dgemm demo...\n");
    printf("A =");    printmat(n, n, A, n);    printf("\n");
    printf("B =");    printmat(n, n, B, n);    printf("\n");
    printf("C =");    printmat(n, n, C, n);    printf("\n");
    alpha = 3.0;
    beta = -2.0;
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, alpha, A, n, B, n, beta, C, n);
    printf("alpha = %5.3e\n", alpha);
    printf("beta  = %5.3e\n", beta);
    printf("ans=");    printmat(n, n, C, n);    printf("\n");
    printf("#check by Matlab/Octave by:\n");
    printf("alpha * A * B + beta * C\n");
    delete[]C;    delete[]B;    delete[]A;
}
