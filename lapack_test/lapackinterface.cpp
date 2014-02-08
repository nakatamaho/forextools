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

#include <cblas.h>
#include <lapacke.h>

//----
#define _DLLAPI extern "C" __declspec(dllexport)

_DLLAPI void __stdcall mql_dgemm(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, int M, int N, int K, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc)
{
  cblas_dgemm(Order, TransA, TransB, (blasint)M, (blasint)N, (blasint)K, alpha, A, (blasint)lda, B, (blasint)ldb, beta, C, (blasint)ldc);
}

_DLLAPI double __stdcall mql_ddot(int n, double *x, int incx, double *y, int incy)
{
  return cblas_ddot(n, x, incx, y, incy);
}

_DLLAPI double __stdcall mql_dasm(int n, double *x, int incx)
{
  return cblas_dasum(n, x, incx);
}

_DLLAPI int __stdcall mql_dsyev(int matrix_order, char jobz, char uplo, lapack_int n, double* a, lapack_int lda, double* w)
{
  return LAPACKE_dsyev(matrix_order, jobz, uplo, n, a, lda, w);
}

_DLLAPI int __stdcall mql_dgesvd(int matrix_order, char jobu, char jobvt, lapack_int m, lapack_int n, double* a, lapack_int lda, double* s, double* u, lapack_int ldu, double* vt, lapack_int ldvt, double* superb)
{
  return LAPACKE_dgesvd(matrix_order, jobu, jobvt, (lapack_int)m, (lapack_int)n, a, (lapack_int)lda, s, u, (lapack_int)ldu, vt, (lapack_int)ldvt, superb);
}

_DLLAPI int __stdcall mql_dgesdd(int matrix_order, char jobz, int m, int n, double* a, int lda, double* s, double* u, int ldu, double* vt, int ldvt)
{
  return LAPACKE_dgesdd(matrix_order, jobz, (lapack_int)m, (lapack_int)n, a, (lapack_int)lda, s, u, (lapack_int)ldu, vt, (lapack_int)ldvt);
}
