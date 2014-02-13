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

_DLLAPI double __stdcall mql_dasm(int n, double *x, int incx)
{
    return dasum_f77(&n, x, &incx);
}

_DLLAPI void __stdcall mql_daxpy(int n, double da, double *dx, int incx, double *dy, int incy)
{
    daxpy_f77(&n, &da, dx, &incx, dy, &incy);
    return;
}

_DLLAPI void __stdcall mql_dcopy(int n, double *dx, int incx, double *dy, int incy)
{
    dcopy_f77(&n, dx, &incx, dy, &incy);
    return;
}

_DLLAPI double __stdcall mql_ddot(int n, double *x, int incx, double *y, int incy)
{
    return ddot_f77(&n, x, &incx, y, &incy);
}

_DLLAPI void __stdcall mql_dgbmv(char trans, int m, int n, int kl, int ku, double alpha, double *a, int lda, double *x, int incx, double beta, double *y, int incy)
{
    dgbmv_f77(&trans, &m, &n, &kl, &ku, &alpha, a, &lda, x, &incx, &beta, y, &incy);
    return;
}

_DLLAPI void __stdcall mql_dgemm(char TransA, char TransB, int M, int N, int K, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc)
{
    dgemm_f77(&TransA, &TransB, &M, &N, &K, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}

_DLLAPI void __stdcall mql_dgemv(char trans, int m, int n, double alpha, double *a, int lda, double *x, int incx, double beta, double *y, int incy)
{
    dgemv_f77(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

_DLLAPI void __stdcall mql_dger(int m, int n, double alpha, double *x, int incx, double *y, int incy, double *a, int lda)
{
    dger_f77(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
}

_DLLAPI double __stdcall mql_dnrm2(int n, double *x, int incx)
{
    return dnrm2_f77(&n, x, &incx);
}

_DLLAPI void __stdcall mql_drot(int n, double *dx, int incx, double *dy, int incy, double c, double s)
{
   drot_f77(&n, dx, &incx, dy, &incy, &c, &s);
}

_DLLAPI void __stdcall mql_drotg(double *da, double *db, double *c, double *s)
{
    drotg_f77(da, db, c, s);
}

_DLLAPI void __stdcall mql_drotm(int n, double *dx, int incx, double *dy, int incy, double *dparam)
{
    drotm_f77(&n, dx, &incx, dy, &incy, dparam);
}

_DLLAPI void __stdcall mql_drotmg(double *dd1, double *dd2, double *dx1, double dy1, double *dparam)
{
    drotmg_f77(dd1, dd2, dx1, &dy1, dparam);
}

_DLLAPI void __stdcall mql_dsbmv(char uplo, int n, int k, double alpha, double *a, int lda, double *x, int incx, double beta, double *y, int incy)
{
    dsbmv_f77(&uplo, &n, &k, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

_DLLAPI void __stdcall mql_dscal(int n, double *da, double *dx, int incx)
{
    dscal_f77(&n, da, dx, &incx);
}

_DLLAPI void __stdcall mql_dspmv(char uplo, int n, double alpha, double *ap, double *x, int incx, double beta, double *y, int incy)
{
    dspmv_f77(&uplo, &n, &alpha, ap, x, &incx, &beta, y, &incy);
}

_DLLAPI void __stdcall mql_dspr(char uplo, int n, double alpha, double *x, int incx, double *ap)
{
    dspr_f77(&uplo, &n, &alpha, x, &incx, ap);
}

_DLLAPI void __stdcall mql_dspr2(char uplo, int n, double alpha, double *x, int incx, double *y, int incy, double *ap)
{
    dspr2_f77(&uplo, &n, &alpha, x, &incx, y, &incy, ap);
}

_DLLAPI void __stdcall mql_dswap(int n, double *dx, int incx, double *dy, int incy)
{
    dswap_f77(&n, dx, &incx, dy, &incy);
}

_DLLAPI void __stdcall mql_dsymm(char side, char uplo, int m, int n, double alpha, double *a, int lda, double *b, int ldb, double beta, double *c, int ldc)
{
    dsymm_f77(&side, &uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

_DLLAPI void __stdcall mql_dsymv(char uplo, int n, double alpha, double *a, int lda, double *x, int incx, double beta, double *y, int incy)
{
    dsymv_f77(&uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

_DLLAPI void __stdcall mql_dsyr(char uplo, int n, double alpha, double *x, int incx, double *a, int lda)
{
    dsyr_f77(&uplo, &n, &alpha, x, &incx, a, &lda);
}

_DLLAPI void __stdcall mql_dsyr2(char uplo, int n, double alpha, double *x, int incx, double *y, int incy, double *a, int lda)
{
    dsyr2_f77(&uplo, &n, &alpha, x, &incx, y, &incy, a, &lda);
}

_DLLAPI void __stdcall mql_dsyr2k(char uplo, char trans, int n, int k, double alpha, double *a, int lda, double *b, int ldb, double beta, double *c, int ldc)
{
    dsyr2k_f77(&uplo, &trans, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

_DLLAPI void __stdcall mql_dsyrk(char uplo, char trans, int n, int k, double alpha, double *a, int lda, double beta, double *c, int ldc)
{
    dsyrk_f77(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
}

_DLLAPI void __stdcall mql_dtbmv(char uplo, char trans, char diag, int n, int k, double *a, int lda, double *x, int incx)
{
    dtbmv_f77(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx);
}

_DLLAPI void __stdcall mql_dtbsv(char uplo, char trans, char diag, int n, int k, double *a, int lda, double *x, int incx)
{
    dtbsv_f77(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx);
}

_DLLAPI void __stdcall mql_dtpmv(char uplo, char trans, char diag, int n, double *ap, double *x, int incx)
{
    dtpmv_f77(&uplo, &trans, &diag, &n, ap, x, &incx);
}

_DLLAPI void __stdcall mql_dtpsv(char uplo, char trans, char diag, int n, double *ap, double *x, int incx)
{
    dtpsv_f77(&uplo, &trans, &diag, &n, ap, x, &incx);
}

_DLLAPI void __stdcall mql_dtrmm(char side, char uplo, char transa, char diag, int m, int n, double alpha, double *a, int lda, double *b, int ldb)
{
    dtrmm_f77(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
}

_DLLAPI void __stdcall mql_dtrmv(char uplo, char trans, char diag, int n, double *a, int lda, double *x, int incx)
{
    dtrmv_f77(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
}

_DLLAPI void __stdcall mql_dtrsm(char side, char uplo, char transa, char diag, int m, int n, double alpha, double *a, int lda, double *b, int ldb)
{
    dtrsm_f77(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb);
}

_DLLAPI void __stdcall mql_dtrsv(char uplo, char trans, char diag, int n, double *a, int lda, double *x, int incx)
{
    dtrsv_f77(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
}

_DLLAPI void __stdcall mql_idamax(int n, double *dx, int incx)
{
    idamax_f77(&n, dx, &incx);
}

_DLLAPI int __stdcall mql_dsyev(char jobz, char uplo, lapack_int n, double *a, lapack_int lda, double *w)
{
    int matrix_order = LAPACK_COL_MAJOR;
    return LAPACKE_dsyev(matrix_order, jobz, uplo, (lapack_int) n, a, (lapack_int) lda, w);
}

_DLLAPI int __stdcall mql_dgesvd(char jobu, char jobvt, lapack_int m, lapack_int n, double *a, lapack_int lda, double *s, double *u, lapack_int ldu, double *vt, lapack_int ldvt, double *superb)
{
    int matrix_order = LAPACK_COL_MAJOR;
    return LAPACKE_dgesvd(matrix_order, jobu, jobvt, (lapack_int) m, (lapack_int) n, a, (lapack_int) lda, s, u, (lapack_int) ldu, vt, (lapack_int) ldvt, superb);
}

_DLLAPI int __stdcall mql_dgesdd(char jobz, int m, int n, double *a, int lda, double *s, double *u, int ldu, double *vt, int ldvt)
{
    int matrix_order = LAPACK_COL_MAJOR;
    return LAPACKE_dgesdd(matrix_order, jobz, (lapack_int) m, (lapack_int) n, a, (lapack_int) lda, s, u, (lapack_int) ldu, vt, (lapack_int) ldvt);
}

_DLLAPI int __stdcall mql_dgels(char trans, int m, int n, int nrhs, double *a, int lda, double *b, int ldb)
{
    int matrix_order = LAPACK_COL_MAJOR;
    return LAPACKE_dgels(matrix_order, trans, (lapack_int) m, (lapack_int) n, (lapack_int) nrhs, a, (lapack_int) lda, b, (lapack_int) ldb);
}
