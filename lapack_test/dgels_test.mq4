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

#property copyright   "2014- NAKATA Maho"
#property link        "https://github.com/nakatamaho/forextools"

#property indicator_chart_window

#include "mql_lapack.mqh"

//Matlab/Octave format
void printmat(int N, int M, double &A[], int LDA)
{
    double mtmp;
    string STR;
    STR = "[ ";
    for (int i = 0; i < N; i++) {
	STR = STR + ("[ ");
	for (int j = 0; j < M; j++) {
	    mtmp = A[i + j * LDA];
	    STR = STR + StringFormat("%5.2e", mtmp);
	    if (j < M - 1)
		STR = STR + (", ");
	}
	if (i < N - 1)
	    STR = STR + ("]; ");
	else
	    STR = STR + ("] ");
    }
    STR = STR + ("]");
    Print(STR);
}

int init()
{
    int m = 5, n = 3, nrhs = 2, lda = 5, ldb = 5;
    int ret;
    double A[], B[];

    ArrayResize(A, m * n);
    ArrayResize(B, m * nrhs);

    //setting A matrix
    A[0 + 0 * lda] = 1;    A[0 + 1 * lda] = 1;    A[0 + 2 * lda] = 1;
    A[1 + 0 * lda] = 2;    A[1 + 1 * lda] = 3;    A[1 + 2 * lda] = 4;
    A[2 + 0 * lda] = 3;    A[2 + 1 * lda] = 5;    A[2 + 2 * lda] =  2;
    A[3 + 0 * lda] = 4;    A[3 + 1 * lda] = 2;    A[3 + 2 * lda] =  5;
    A[4 + 0 * lda] = 5;    A[4 + 1 * lda] = 4;    A[4 + 2 * lda] =  3;

    //setting A matrix
    B[0 + 0 * ldb] = -10;    B[0 + 1 * ldb] = -3;
    B[1 + 0 * ldb] =  12;    B[1 + 1 * ldb] = 14;
    B[2 + 0 * ldb] =  14;    B[2 + 1 * ldb] = 12;
    B[3 + 0 * ldb] =  16;    B[3 + 1 * ldb] = 16;
    B[4 + 0 * ldb] =  18;    B[4 + 1 * ldb] = 16;

    printf("A =");    printmat(m, n, A, lda);    printf("\n");
    printf("B =");    printmat(m, nrhs, B, ldb); printf("\n");

    mql_dgels(LAPACK_COL_MAJOR,'N', m, n, nrhs, A, lda, B, ldb);

    printf("B =");    printmat(n, nrhs, B, ldb); printf("\n");
    return (0);
}

int deinit()
{
    return (0);
}

int start()
{
    return (0);
}
