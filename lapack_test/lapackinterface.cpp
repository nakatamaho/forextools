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
//#property copyright "NAKATA Maho"
//#property link      "https://github.com/nakatamaho"

#define WIN32_LEAN_AND_MEAN	// Exclude rarely-used stuff from Windows headers
#include <windows.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "blas_subset.h"
#include <lapacke.h>

#define _DLLAPI extern "C" __declspec(dllexport)

_DLLAPI void __stdcall mql_dgemm(char TransA, char TransB, int M, int N, int K, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc)
{
 dgemm_f77(&TransA, &TransB, &M, &N, &K, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}

_DLLAPI double __stdcall mql_ddot(int n, double *x, int incx, double *y, int incy)
{
  return ddot_f77(&n, x, &incx, y, &incy);
}

_DLLAPI double __stdcall mql_dasm(int n, double *x, int incx)
{
  return dasum_f77(&n, x, &incx);
}

_DLLAPI int __stdcall mql_dsyev(char jobz, char uplo, lapack_int n, double* a, lapack_int lda, double* w)
{
  int matrix_order=LAPACK_COL_MAJOR;
  return LAPACKE_dsyev(matrix_order, jobz, uplo, (lapack_int)n, a, (lapack_int)lda, w);
}

_DLLAPI int __stdcall mql_dgesvd(char jobu, char jobvt, lapack_int m, lapack_int n, double* a, lapack_int lda, double* s, double* u, lapack_int ldu, double* vt, lapack_int ldvt, double* superb)
{
  int matrix_order=LAPACK_COL_MAJOR;
  return LAPACKE_dgesvd(matrix_order, jobu, jobvt, (lapack_int)m, (lapack_int)n, a, (lapack_int)lda, s, u, (lapack_int)ldu, vt, (lapack_int)ldvt, superb);
}

_DLLAPI int __stdcall mql_dgesdd(char jobz, int m, int n, double* a, int lda, double* s, double* u, int ldu, double* vt, int ldvt)
{
  int matrix_order=LAPACK_COL_MAJOR;
  return LAPACKE_dgesdd(matrix_order, jobz, (lapack_int)m, (lapack_int)n, a, (lapack_int)lda, s, u, (lapack_int)ldu, vt, (lapack_int)ldvt);
}

_DLLAPI int __stdcall mql_dgels(char trans, int m, int n, int nrhs, double *a, int lda, double *b, int ldb)
{
  int matrix_order=LAPACK_COL_MAJOR;
  return LAPACKE_dgels(matrix_order, trans, (lapack_int)m, (lapack_int)n, (lapack_int) nrhs, a, (lapack_int) lda, b, (lapack_int)ldb );
}
