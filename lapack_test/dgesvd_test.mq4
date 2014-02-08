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
#define min(X,Y) ((X) < (Y) ? (X) : (Y))

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
    int m = 3, n = 4, ret;
    double A[], U[], VT[], S[], superb[];

    ArrayResize(A, m * n);
    ArrayResize(U, m * m);
    ArrayResize(VT, n * n);
    ArrayResize(S, min(m, n));
    ArrayResize(superb, min(m, n));

    //setting A matrix
    A[0 + 0 * m] = 1;    A[0 + 1 * m] = 2;    A[0 + 2 * m] = 3;    A[0 + 3 * m] = 4;
    A[1 + 0 * m] = 5;    A[1 + 1 * m] = 6;    A[1 + 2 * m] = 7;    A[1 + 3 * m] = 8;
    A[2 + 0 * m] = 9;    A[2 + 1 * m] = 10;   A[2 + 2 * m] = 11;   A[2 + 3 * m] = 12;

    printf("A =");    printmat(m, n, A, m);    printf("\n");
    ret = mql_dgesvd(LAPACK_COL_MAJOR, 'A', 'A', m, n, A, m, S, U, m, VT, n, superb);

    //print out some results.
    printf("#singular values\n");
    printf("S =");    printmat(min(m, n), 1, S, 1);    printf("\n");
    printf("#U matrix\n");
    printf("U =");    printmat(m, m, U, m);    printf("\n");
    printf("#V^t matrix\n");
    printf("Vt =");   printmat(n, n, VT, n);    printf("\n");
    printf("#You can check result by octave/matlab by\n");
    printf("[u, s, v] = svd (A) \n");
    printf("#u*s*v' should be equal to A\n");
    printf("u*s*v' \n");
    printf("U*diag(S,%d,%d)*Vt \n",m,n);

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
