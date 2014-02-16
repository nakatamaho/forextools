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

_DLLAPI double __stdcall mql_dasum(int n, double *x, int incx)
{
    return dasum_f77(&n, x, &incx);
}

_DLLAPI int __stdcall mql_idamax(int n, double *x, int incx)
{
    return idamax_f77(&n, x, &incx);
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

_DLLAPI int __stdcall mql_dbdsdc(char uplo, char compq, int n, double *d, double *e, double *u, int ldu, double *vt, int ldvt, double *q, int *iq)
{
    return LAPACKE_dbdsdc(LAPACK_COL_MAJOR, uplo, compq, (lapack_int) n, d, e, u, (lapack_int) ldu, vt, (lapack_int) ldvt, q, iq);
}

_DLLAPI int __stdcall mql_dbdsqr(char uplo, int n, int ncvt, int nru, int ncc, double *d, double *e, double *vt, int ldvt, double *u, int ldu, double *c, int ldc)
{
    return LAPACKE_dbdsqr(LAPACK_COL_MAJOR, uplo, (lapack_int) n, (lapack_int) ncvt, (lapack_int) nru, (lapack_int) ncc, d, e, vt, (lapack_int) ldvt, u, (lapack_int) ldu, c, (lapack_int) ldc);
}

_DLLAPI int __stdcall mql_ddisna(char job, int m, int n, double *d, double *sep)
{
    return LAPACKE_ddisna(job, (lapack_int) m, (lapack_int) n, d, sep);
}

_DLLAPI int __stdcall mql_dgbbrd(char vect, int m, int n, int ncc, int kl, int ku, double *ab, int ldab, double *d, double *e, double *q, int ldq, double *pt, int ldpt, double *c, int ldc)
{
    return LAPACKE_dgbbrd(LAPACK_COL_MAJOR, vect, (lapack_int) m, (lapack_int) n, (lapack_int) ncc, (lapack_int) kl, (lapack_int) ku, ab, (lapack_int) ldab, d, e, q, (lapack_int) ldq, pt, (lapack_int) ldpt, c, (lapack_int) ldc);
}

_DLLAPI int __stdcall mql_dgbcon(char norm, int n, int kl, int ku, double *ab, int ldab, int *ipiv, double anorm, double *rcond)
{
    return LAPACKE_dgbcon(LAPACK_COL_MAJOR, norm, (lapack_int) n, (lapack_int) kl, (lapack_int) ku, ab, (lapack_int) ldab, ipiv, anorm, rcond);
}

_DLLAPI int __stdcall mql_dgbequ(int m, int n, int kl, int ku, double *ab, int ldab, double *r, double *c, double *rowcnd, double *colcnd, double *amax)
{
    return LAPACKE_dgbequ(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) n, (lapack_int) kl, (lapack_int) ku, ab, (lapack_int) ldab, r, c, rowcnd, colcnd, amax);
}

_DLLAPI int __stdcall mql_dgbequb(int m, int n, int kl, int ku, double *ab, int ldab, double *r, double *c, double *rowcnd, double *colcnd, double *amax)
{
    return LAPACKE_dgbequb(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) n, (lapack_int) kl, (lapack_int) ku, ab, (lapack_int) ldab, r, c, rowcnd, colcnd, amax);
}

_DLLAPI int __stdcall mql_dgbrfs(char trans, int n, int kl, int ku, int nrhs, double *ab, int ldab, double *afb, int ldafb, int *ipiv, double *b, int ldb, double *x, int ldx, double *ferr, double *berr)
{
    return LAPACKE_dgbrfs(LAPACK_COL_MAJOR, trans, (lapack_int) n, (lapack_int) kl, (lapack_int) ku, (lapack_int) nrhs, ab, (lapack_int) ldab, afb, (lapack_int) ldafb, ipiv, b, (lapack_int) ldb, x, (lapack_int) ldx, ferr, berr);
}

#if defined USEXBLAS
_DLLAPI int __stdcall mql_dgbrfsx(char trans, char equed, int n, int kl, int ku, int nrhs, double *ab, int ldab, double *afb, int ldafb, int *ipiv, double *r, double *c, double *b, int ldb, double *x, int ldx, double *rcond, double *berr, int n_err_bnds, double *err_bnds_norm, double *err_bnds_comp, int nparams, double *params)
{
    return LAPACKE_dgbrfsx(LAPACK_COL_MAJOR, trans, equed, (lapack_int) n, (lapack_int) kl, (lapack_int) ku, (lapack_int) nrhs, ab, (lapack_int) ldab, afb, (lapack_int) ldafb, ipiv, r, c, b, (lapack_int) ldb, x, (lapack_int) ldx, rcond, berr, (lapack_int) n_err_bnds, err_bnds_norm, err_bnds_comp, (lapack_int) nparams, params);
}
#endif

_DLLAPI int __stdcall mql_dgbsv(int n, int kl, int ku, int nrhs, double *ab, int ldab, int *ipiv, double *b, int ldb)
{
    return LAPACKE_dgbsv(LAPACK_COL_MAJOR, (lapack_int) n, (lapack_int) kl, (lapack_int) ku, (lapack_int) nrhs, ab, (lapack_int) ldab, ipiv, b, (lapack_int) ldb);
}

_DLLAPI int __stdcall mql_dgbsvx(char fact, char trans, int n, int kl, int ku, int nrhs, double *ab, int ldab, double *afb, int ldafb, int *ipiv, char *equed, double *r, double *c, double *b, int ldb, double *x, int ldx, double *rcond, double *ferr, double *berr, double *rpivot)
{
    return LAPACKE_dgbsvx(LAPACK_COL_MAJOR, fact, trans, (lapack_int) n, (lapack_int) kl, (lapack_int) ku, (lapack_int) nrhs, ab, (lapack_int) ldab, afb, (lapack_int) ldafb, ipiv, equed, r, c, b, (lapack_int) ldb, x, (lapack_int) ldx, rcond, ferr, berr, rpivot);
}

#if defined USEXBLAS
_DLLAPI int __stdcall mql_dgbsvxx(char fact, char trans, int n, int kl, int ku, int nrhs, double *ab, int ldab, double *afb, int ldafb, int *ipiv, char *equed, double *r, double *c, double *b, int ldb, double *x, int ldx, double *rcond, double *rpvgrw, double *berr, int n_err_bnds, double *err_bnds_norm, double *err_bnds_comp, int nparams, double *params)
{
    return LAPACKE_dgbsvxx(LAPACK_COL_MAJOR, fact, trans, (lapack_int) n, (lapack_int) kl, (lapack_int) ku, (lapack_int) nrhs, ab, (lapack_int) ldab, afb, (lapack_int) ldafb, ipiv, equed, r, c, b, (lapack_int) ldb, x, (lapack_int) ldx, rcond, rpvgrw, berr, (lapack_int) n_err_bnds, err_bnds_norm, err_bnds_comp, (lapack_int) nparams, params);
}
#endif

_DLLAPI int __stdcall mql_dgbtrf(int m, int n, int kl, int ku, double *ab, int ldab, int *ipiv)
{
    return LAPACKE_dgbtrf(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) n, (lapack_int) kl, (lapack_int) ku, ab, (lapack_int) ldab, ipiv);
}

_DLLAPI int __stdcall mql_dgbtrs(char trans, int n, int kl, int ku, int nrhs, double *ab, int ldab, int *ipiv, double *b, int ldb)
{
    return LAPACKE_dgbtrs(LAPACK_COL_MAJOR, trans, (lapack_int) n, (lapack_int) kl, (lapack_int) ku, (lapack_int) nrhs, ab, (lapack_int) ldab, ipiv, b, (lapack_int) ldb);
}

_DLLAPI int __stdcall mql_dgebak(char job, char side, int n, int ilo, int ihi, double *scale, int m, double *v, int ldv)
{
    return LAPACKE_dgebak(LAPACK_COL_MAJOR, job, side, (lapack_int) n, (lapack_int) ilo, (lapack_int) ihi, scale, (lapack_int) m, v, (lapack_int) ldv);
}

_DLLAPI int __stdcall mql_dgebal(char job, int n, double *a, int lda, int *ilo, int *ihi, double *scale)
{
    return LAPACKE_dgebal(LAPACK_COL_MAJOR, job, (lapack_int) n, a, (lapack_int) lda, ilo, ihi, scale);
}

_DLLAPI int __stdcall mql_dgebrd(int m, int n, double *a, int lda, double *d, double *e, double *tauq, double *taup)
{
    return LAPACKE_dgebrd(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) n, a, (lapack_int) lda, d, e, tauq, taup);
}

_DLLAPI int __stdcall mql_dgecon(char norm, int n, double *a, int lda, double anorm, double *rcond)
{
    return LAPACKE_dgecon(LAPACK_COL_MAJOR, norm, (lapack_int) n, a, (lapack_int) lda, anorm, rcond);
}

_DLLAPI int __stdcall mql_dgeequ(int m, int n, double *a, int lda, double *r, double *c, double *rowcnd, double *colcnd, double *amax)
{
    return LAPACKE_dgeequ(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) n, a, (lapack_int) lda, r, c, rowcnd, colcnd, amax);
}

_DLLAPI int __stdcall mql_dgeequb(int m, int n, double *a, int lda, double *r, double *c, double *rowcnd, double *colcnd, double *amax)
{
    return LAPACKE_dgeequb(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) n, a, (lapack_int) lda, r, c, rowcnd, colcnd, amax);
}

_DLLAPI int __stdcall mql_dgeev(char jobvl, char jobvr, int n, double *a, int lda, double *wr, double *wi, double *vl, int ldvl, double *vr, int ldvr)
{
    return LAPACKE_dgeev(LAPACK_COL_MAJOR, jobvl, jobvr, (lapack_int) n, a, (lapack_int) lda, wr, wi, vl, (lapack_int) ldvl, vr, (lapack_int) ldvr);
}

_DLLAPI int __stdcall mql_dgeevx(char balanc, char jobvl, char jobvr, char sense, int n, double *a, int lda, double *wr, double *wi, double *vl, int ldvl, double *vr, int ldvr, int *ilo, int *ihi, double *scale, double *abnrm, double *rconde, double *rcondv)
{
    return LAPACKE_dgeevx(LAPACK_COL_MAJOR, balanc, jobvl, jobvr, sense, (lapack_int) n, a, (lapack_int) lda, wr, wi, vl, (lapack_int) ldvl, vr, (lapack_int) ldvr, ilo, ihi, scale, abnrm, rconde, rcondv);
}

_DLLAPI int __stdcall mql_dgehrd(int n, int ilo, int ihi, double *a, int lda, double *tau)
{
    return LAPACKE_dgehrd(LAPACK_COL_MAJOR, (lapack_int) n, (lapack_int) ilo, (lapack_int) ihi, a, (lapack_int) lda, tau);
}

_DLLAPI int __stdcall mql_dgejsv(char joba, char jobu, char jobv, char jobr, char jobt, char jobp, int m, int n, double *a, int lda, double *sva, double *u, int ldu, double *v, int ldv, double *stat, int *istat)
{
    return LAPACKE_dgejsv(LAPACK_COL_MAJOR, joba, jobu, jobv, jobr, jobt, jobp, (lapack_int) m, (lapack_int) n, a, (lapack_int) lda, sva, u, (lapack_int) ldu, v, (lapack_int) ldv, stat, istat);
}

_DLLAPI int __stdcall mql_dgelq2(int m, int n, double *a, int lda, double *tau)
{
    return LAPACKE_dgelq2(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) n, a, (lapack_int) lda, tau);
}

_DLLAPI int __stdcall mql_dgelqf(int m, int n, double *a, int lda, double *tau)
{
    return LAPACKE_dgelqf(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) n, a, (lapack_int) lda, tau);
}

_DLLAPI int __stdcall mql_dgels(char trans, int m, int n, int nrhs, double *a, int lda, double *b, int ldb)
{
    return LAPACKE_dgels(LAPACK_COL_MAJOR, trans, (lapack_int) m, (lapack_int) n, (lapack_int) nrhs, a, (lapack_int) lda, b, (lapack_int) ldb);
}

_DLLAPI int __stdcall mql_dgelsd(int m, int n, int nrhs, double *a, int lda, double *b, int ldb, double *s, double rcond, int *rank)
{
    return LAPACKE_dgelsd(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) n, (lapack_int) nrhs, a, (lapack_int) lda, b, (lapack_int) ldb, s, rcond, rank);
}

_DLLAPI int __stdcall mql_dgelss(int m, int n, int nrhs, double *a, int lda, double *b, int ldb, double *s, double rcond, int *rank)
{
    return LAPACKE_dgelss(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) n, (lapack_int) nrhs, a, (lapack_int) lda, b, (lapack_int) ldb, s, rcond, rank);
}

_DLLAPI int __stdcall mql_dgelsy(int m, int n, int nrhs, double *a, int lda, double *b, int ldb, int *jpvt, double rcond, int *rank)
{
    return LAPACKE_dgelsy(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) n, (lapack_int) nrhs, a, (lapack_int) lda, b, (lapack_int) ldb, jpvt, rcond, rank);
}

_DLLAPI int __stdcall mql_dgeqlf(int m, int n, double *a, int lda, double *tau)
{
    return LAPACKE_dgeqlf(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) n, a, (lapack_int) lda, tau);
}

_DLLAPI int __stdcall mql_dgeqp3(int m, int n, double *a, int lda, int *jpvt, double *tau)
{
    return LAPACKE_dgeqp3(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) n, a, (lapack_int) lda, jpvt, tau);
}

_DLLAPI int __stdcall mql_dgeqpf(int m, int n, double *a, int lda, int *jpvt, double *tau)
{
    return LAPACKE_dgeqpf(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) n, a, (lapack_int) lda, jpvt, tau);
}

_DLLAPI int __stdcall mql_dgeqr2(int m, int n, double *a, int lda, double *tau)
{
    return LAPACKE_dgeqr2(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) n, a, (lapack_int) lda, tau);
}

_DLLAPI int __stdcall mql_dgeqrf(int m, int n, double *a, int lda, double *tau)
{
    return LAPACKE_dgeqrf(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) n, a, (lapack_int) lda, tau);
}

_DLLAPI int __stdcall mql_dgeqrfp(int m, int n, double *a, int lda, double *tau)
{
    return LAPACKE_dgeqrfp(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) n, a, (lapack_int) lda, tau);
}

_DLLAPI int __stdcall mql_dgerfs(char trans, int n, int nrhs, double *a, int lda, double *af, int ldaf, int *ipiv, double *b, int ldb, double *x, int ldx, double *ferr, double *berr)
{
    return LAPACKE_dgerfs(LAPACK_COL_MAJOR, trans, (lapack_int) n, (lapack_int) nrhs, a, (lapack_int) lda, af, (lapack_int) ldaf, ipiv, b, (lapack_int) ldb, x, (lapack_int) ldx, ferr, berr);
}

#if defined USEXBLAS
_DLLAPI int __stdcall mql_dgerfsx(char trans, char equed, int n, int nrhs, double *a, int lda, double *af, int ldaf, int *ipiv, double *r, double *c, double *b, int ldb, double *x, int ldx, double *rcond, double *berr, int n_err_bnds, double *err_bnds_norm, double *err_bnds_comp, int nparams, double *params)
{
    return LAPACKE_dgerfsx(LAPACK_COL_MAJOR, trans, equed, (lapack_int) n, (lapack_int) nrhs, a, (lapack_int) lda, af, (lapack_int) ldaf, ipiv, r, c, b, (lapack_int) ldb, x, (lapack_int) ldx, rcond, berr, (lapack_int) n_err_bnds, err_bnds_norm, err_bnds_comp, (lapack_int) nparams, params);
}
#endif

_DLLAPI int __stdcall mql_dgerqf(int m, int n, double *a, int lda, double *tau)
{
    return LAPACKE_dgerqf(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) n, a, (lapack_int) lda, tau);
}

_DLLAPI int __stdcall mql_dgesdd(char jobz, int m, int n, double *a, int lda, double *s, double *u, int ldu, double *vt, int ldvt)
{
    return LAPACKE_dgesdd(LAPACK_COL_MAJOR, jobz, (lapack_int) m, (lapack_int) n, a, (lapack_int) lda, s, u, (lapack_int) ldu, vt, (lapack_int) ldvt);
}

_DLLAPI int __stdcall mql_dgesv(int n, int nrhs, double *a, int lda, int *ipiv, double *b, int ldb)
{
    return LAPACKE_dgesv(LAPACK_COL_MAJOR, (lapack_int) n, (lapack_int) nrhs, a, (lapack_int) lda, ipiv, b, (lapack_int) ldb);
}

_DLLAPI int __stdcall mql_dsgesv(int n, int nrhs, double *a, int lda, int *ipiv, double *b, int ldb, double *x, int ldx, int *iter)
{
    return LAPACKE_dsgesv(LAPACK_COL_MAJOR, (lapack_int) n, (lapack_int) nrhs, a, (lapack_int) lda, ipiv, b, (lapack_int) ldb, x, (lapack_int) ldx, iter);
}

_DLLAPI int __stdcall mql_dgesvd(char jobu, char jobvt, int m, int n, double *a, int lda, double *s, double *u, int ldu, double *vt, int ldvt, double *superb)
{
    return LAPACKE_dgesvd(LAPACK_COL_MAJOR, jobu, jobvt, (lapack_int) m, (lapack_int) n, a, (lapack_int) lda, s, u, (lapack_int) ldu, vt, (lapack_int) ldvt, superb);
}

_DLLAPI int __stdcall mql_dgesvj(char joba, char jobu, char jobv, int m, int n, double *a, int lda, double *sva, int mv, double *v, int ldv, double *stat)
{
    return LAPACKE_dgesvj(LAPACK_COL_MAJOR, joba, jobu, jobv, (lapack_int) m, (lapack_int) n, a, (lapack_int) lda, sva, (lapack_int) mv, v, (lapack_int) ldv, stat);
}

_DLLAPI int __stdcall mql_dgesvx(char fact, char trans, int n, int nrhs, double *a, int lda, double *af, int ldaf, int *ipiv, char *equed, double *r, double *c, double *b, int ldb, double *x, int ldx, double *rcond, double *ferr, double *berr, double *rpivot)
{
    return LAPACKE_dgesvx(LAPACK_COL_MAJOR, fact, trans, (lapack_int) n, (lapack_int) nrhs, a, (lapack_int) lda, af, (lapack_int) ldaf, ipiv, equed, r, c, b, (lapack_int) ldb, x, (lapack_int) ldx, rcond, ferr, berr, rpivot);
}

#if defined USEXBLAS
_DLLAPI int __stdcall mql_dgesvxx(char fact, char trans, int n, int nrhs, double *a, int lda, double *af, int ldaf, int *ipiv, char *equed, double *r, double *c, double *b, int ldb, double *x, int ldx, double *rcond, double *rpvgrw, double *berr, int n_err_bnds, double *err_bnds_norm, double *err_bnds_comp, int nparams, double *params)
{
    return LAPACKE_dgesvxx(LAPACK_COL_MAJOR, fact, trans, (lapack_int) n, (lapack_int) nrhs, a, (lapack_int) lda, af, (lapack_int) ldaf, ipiv, equed, r, c, b, (lapack_int) ldb, x, (lapack_int) ldx, rcond, rpvgrw, berr, (lapack_int) n_err_bnds, err_bnds_norm, err_bnds_comp, (lapack_int) nparams, params);
}
#endif

_DLLAPI int __stdcall mql_dgetf2(int m, int n, double *a, int lda, int *ipiv)
{
    return LAPACKE_dgetf2(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) n, a, (lapack_int) lda, ipiv);
}

_DLLAPI int __stdcall mql_dgetrf(int m, int n, double *a, int lda, int *ipiv)
{
    return LAPACKE_dgetrf(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) n, a, (lapack_int) lda, ipiv);
}

_DLLAPI int __stdcall mql_dgetri(int n, double *a, int lda, int *ipiv)
{
    return LAPACKE_dgetri(LAPACK_COL_MAJOR, (lapack_int) n, a, (lapack_int) lda, ipiv);
}

_DLLAPI int __stdcall mql_dgetrs(char trans, int n, int nrhs, double *a, int lda, int *ipiv, double *b, int ldb)
{
    return LAPACKE_dgetrs(LAPACK_COL_MAJOR, trans, (lapack_int) n, (lapack_int) nrhs, a, (lapack_int) lda, ipiv, b, (lapack_int) ldb);
}

_DLLAPI int __stdcall mql_dggbak(char job, char side, int n, int ilo, int ihi, double *lscale, double *rscale, int m, double *v, int ldv)
{
    return LAPACKE_dggbak(LAPACK_COL_MAJOR, job, side, (lapack_int) n, (lapack_int) ilo, (lapack_int) ihi, lscale, rscale, (lapack_int) m, v, (lapack_int) ldv);
}

_DLLAPI int __stdcall mql_dggbal(char job, int n, double *a, int lda, double *b, int ldb, int *ilo, int *ihi, double *lscale, double *rscale)
{
    return LAPACKE_dggbal(LAPACK_COL_MAJOR, job, (lapack_int) n, a, (lapack_int) lda, b, (lapack_int) ldb, ilo, ihi, lscale, rscale);
}

_DLLAPI int __stdcall mql_dggev(char jobvl, char jobvr, int n, double *a, int lda, double *b, int ldb, double *alphar, double *alphai, double *beta, double *vl, int ldvl, double *vr, int ldvr)
{
    return LAPACKE_dggev(LAPACK_COL_MAJOR, jobvl, jobvr, (lapack_int) n, a, (lapack_int) lda, b, (lapack_int) ldb, alphar, alphai, beta, vl, (lapack_int) ldvl, vr, (lapack_int) ldvr);
}

_DLLAPI int __stdcall mql_dggevx(char balanc, char jobvl, char jobvr, char sense, int n, double *a, int lda, double *b, int ldb, double *alphar, double *alphai, double *beta, double *vl, int ldvl, double *vr, int ldvr, int *ilo, int *ihi, double *lscale, double *rscale, double *abnrm, double *bbnrm, double *rconde, double *rcondv)
{
    return LAPACKE_dggevx(LAPACK_COL_MAJOR, balanc, jobvl, jobvr, sense, (lapack_int) n, a, (lapack_int) lda, b, (lapack_int) ldb, alphar, alphai, beta, vl, (lapack_int) ldvl, vr, (lapack_int) ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv);
}

_DLLAPI int __stdcall mql_dggglm(int n, int m, int p, double *a, int lda, double *b, int ldb, double *d, double *x, double *y)
{
    return LAPACKE_dggglm(LAPACK_COL_MAJOR, (lapack_int) n, (lapack_int) m, (lapack_int) p, a, (lapack_int) lda, b, (lapack_int) ldb, d, x, y);
}

_DLLAPI int __stdcall mql_dgghrd(char compq, char compz, int n, int ilo, int ihi, double *a, int lda, double *b, int ldb, double *q, int ldq, double *z, int ldz)
{
    return LAPACKE_dgghrd(LAPACK_COL_MAJOR, compq, compz, (lapack_int) n, (lapack_int) ilo, (lapack_int) ihi, a, (lapack_int) lda, b, (lapack_int) ldb, q, (lapack_int) ldq, z, (lapack_int) ldz);
}

_DLLAPI int __stdcall mql_dgglse(int m, int n, int p, double *a, int lda, double *b, int ldb, double *c, double *d, double *x)
{
    return LAPACKE_dgglse(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) n, (lapack_int) p, a, (lapack_int) lda, b, (lapack_int) ldb, c, d, x);
}

_DLLAPI int __stdcall mql_dggqrf(int n, int m, int p, double *a, int lda, double *taua, double *b, int ldb, double *taub)
{
    return LAPACKE_dggqrf(LAPACK_COL_MAJOR, (lapack_int) n, (lapack_int) m, (lapack_int) p, a, (lapack_int) lda, taua, b, (lapack_int) ldb, taub);
}

_DLLAPI int __stdcall mql_dggrqf(int m, int p, int n, double *a, int lda, double *taua, double *b, int ldb, double *taub)
{
    return LAPACKE_dggrqf(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) p, (lapack_int) n, a, (lapack_int) lda, taua, b, (lapack_int) ldb, taub);
}

_DLLAPI int __stdcall mql_dggsvd(char jobu, char jobv, char jobq, int m, int n, int p, int *k, int *l, double *a, int lda, double *b, int ldb, double *alpha, double *beta, double *u, int ldu, double *v, int ldv, double *q, int ldq, int *iwork)
{
    return LAPACKE_dggsvd(LAPACK_COL_MAJOR, jobu, jobv, jobq, (lapack_int) m, (lapack_int) n, (lapack_int) p, k, l, a, (lapack_int) lda, b, (lapack_int) ldb, alpha, beta, u, (lapack_int) ldu, v, (lapack_int) ldv, q, (lapack_int) ldq, iwork);
}

_DLLAPI int __stdcall mql_dggsvp(char jobu, char jobv, char jobq, int m, int p, int n, double *a, int lda, double *b, int ldb, double tola, double tolb, int *k, int *l, double *u, int ldu, double *v, int ldv, double *q, int ldq)
{
    return LAPACKE_dggsvp(LAPACK_COL_MAJOR, jobu, jobv, jobq, (lapack_int) m, (lapack_int) p, (lapack_int) n, a, (lapack_int) lda, b, (lapack_int) ldb, tola, tolb, k, l, u, (lapack_int) ldu, v, (lapack_int) ldv, q, (lapack_int) ldq);
}

_DLLAPI int __stdcall mql_dgtcon(char norm, int n, double *dl, double *d, double *du, double *du2, int *ipiv, double anorm, double *rcond)
{
    return LAPACKE_dgtcon(norm, (lapack_int) n, dl, d, du, du2, ipiv, anorm, rcond);
}

_DLLAPI int __stdcall mql_dgtrfs(char trans, int n, int nrhs, double *dl, double *d, double *du, double *dlf, double *df, double *duf, double *du2, int *ipiv, double *b, int ldb, double *x, int ldx, double *ferr, double *berr)
{
    return LAPACKE_dgtrfs(LAPACK_COL_MAJOR, trans, (lapack_int) n, (lapack_int) nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, (lapack_int) ldb, x, (lapack_int) ldx, ferr, berr);
}

_DLLAPI int __stdcall mql_dgtsv(int n, int nrhs, double *dl, double *d, double *du, double *b, int ldb)
{
    return LAPACKE_dgtsv(LAPACK_COL_MAJOR, (lapack_int) n, (lapack_int) nrhs, dl, d, du, b, (lapack_int) ldb);
}

_DLLAPI int __stdcall mql_dgtsvx(char fact, char trans, int n, int nrhs, double *dl, double *d, double *du, double *dlf, double *df, double *duf, double *du2, int *ipiv, double *b, int ldb, double *x, int ldx, double *rcond, double *ferr, double *berr)
{
    return LAPACKE_dgtsvx(LAPACK_COL_MAJOR, fact, trans, (lapack_int) n, (lapack_int) nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, (lapack_int) ldb, x, (lapack_int) ldx, rcond, ferr, berr);
}

_DLLAPI int __stdcall mql_dgttrf(int n, double *dl, double *d, double *du, double *du2, int *ipiv)
{
    return LAPACKE_dgttrf((lapack_int) n, dl, d, du, du2, ipiv);
}

_DLLAPI int __stdcall mql_dgttrs(char trans, int n, int nrhs, double *dl, double *d, double *du, double *du2, int *ipiv, double *b, int ldb)
{
    return LAPACKE_dgttrs(LAPACK_COL_MAJOR, trans, (lapack_int) n, (lapack_int) nrhs, dl, d, du, du2, ipiv, b, (lapack_int) ldb);
}

_DLLAPI int __stdcall mql_dhgeqz(char job, char compq, char compz, int n, int ilo, int ihi, double *h, int ldh, double *t, int ldt, double *alphar, double *alphai, double *beta, double *q, int ldq, double *z, int ldz)
{
    return LAPACKE_dhgeqz(LAPACK_COL_MAJOR, job, compq, compz, (lapack_int) n, (lapack_int) ilo, (lapack_int) ihi, h, (lapack_int) ldh, t, (lapack_int) ldt, alphar, alphai, beta, q, (lapack_int) ldq, z, (lapack_int) ldz);
}

_DLLAPI int __stdcall mql_dhseqr(char job, char compz, int n, int ilo, int ihi, double *h, int ldh, double *wr, double *wi, double *z, int ldz)
{
    return LAPACKE_dhseqr(LAPACK_COL_MAJOR, job, compz, (lapack_int) n, (lapack_int) ilo, (lapack_int) ihi, h, (lapack_int) ldh, wr, wi, z, (lapack_int) ldz);
}

/*
_DLLAPI int __stdcall mql_dlacn2(int n, double *v, double *x, int *isgn, double *est, int *kase, int *isave)
{
    return LAPACKE_dlacn2((lapack_int) n, v, x, isgn, est, kase, isave);
}
*/

_DLLAPI int __stdcall mql_dlacpy(char uplo, int m, int n, double *a, int lda, double *b, int ldb)
{
    return LAPACKE_dlacpy(LAPACK_COL_MAJOR, uplo, (lapack_int) m, (lapack_int) n, a, (lapack_int) lda, b, (lapack_int) ldb);
}

#if defined TESTING
_DLLAPI int __stdcall mql_dlagge(int m, int n, int kl, int ku, double *d, double *a, int lda, int *iseed)
{
    return LAPACKE_dlagge(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) n, (lapack_int) kl, (lapack_int) ku, d, a, (lapack_int) lda, iseed);
}
#endif

_DLLAPI double __stdcall mql_dlamch(char cmach)
{
    return LAPACKE_dlamch(cmach);
}

_DLLAPI double __stdcall mql_dlange(char norm, int m, int n, double *a, int lda)
{
    return LAPACKE_dlange(LAPACK_COL_MAJOR, norm, (lapack_int) m, (lapack_int) n, a, (lapack_int) lda);
}

_DLLAPI double __stdcall mql_dlansy(char norm, char uplo, int n, double *a, int lda)
{
    return LAPACKE_dlansy(LAPACK_COL_MAJOR, norm, uplo, (lapack_int) n, a, (lapack_int) lda);
}

_DLLAPI double __stdcall mql_dlantr(char norm, char uplo, char diag, int m, int n, double *a, int lda)
{
    return LAPACKE_dlantr(LAPACK_COL_MAJOR, norm, uplo, diag, (lapack_int) m, (lapack_int) n, a, (lapack_int) lda);
}

_DLLAPI int __stdcall mql_dlarfb(char side, char trans, char direct, char storev, int m, int n, int k, double *v, int ldv, double *t, int ldt, double *c, int ldc)
{
    return LAPACKE_dlarfb(LAPACK_COL_MAJOR, side, trans, direct, storev, (lapack_int) m, (lapack_int) n, (lapack_int) k, v, (lapack_int) ldv, t, (lapack_int) ldt, c, (lapack_int) ldc);
}

_DLLAPI int __stdcall mql_dlarfg(int n, double *alpha, double *x, int incx, double *tau)
{
    return LAPACKE_dlarfg((lapack_int) n, alpha, x, (lapack_int) incx, tau);
}

_DLLAPI int __stdcall mql_dlarft(char direct, char storev, int n, int k, double *v, int ldv, double *tau, double *t, int ldt)
{
    return LAPACKE_dlarft(LAPACK_COL_MAJOR, direct, storev, (lapack_int) n, (lapack_int) k, v, (lapack_int) ldv, tau, t, (lapack_int) ldt);
}

_DLLAPI int __stdcall mql_dlarfx(char side, int m, int n, double *v, double tau, double *c, int ldc, double *work)
{
    return LAPACKE_dlarfx(LAPACK_COL_MAJOR, side, (lapack_int) m, (lapack_int) n, v, tau, c, (lapack_int) ldc, work);
}

_DLLAPI int __stdcall mql_dlarnv(int idist, int *iseed, int n, double *x)
{
    return LAPACKE_dlarnv((lapack_int) idist, iseed, (lapack_int) n, x);
}

_DLLAPI int __stdcall mql_dlaset(char uplo, int m, int n, double alpha, double beta, double *a, int lda)
{
    return LAPACKE_dlaset(LAPACK_COL_MAJOR, uplo, (lapack_int) m, (lapack_int) n, alpha, beta, a, (lapack_int) lda);
}

_DLLAPI int __stdcall mql_dlasrt(char id, int n, double *d)
{
    return LAPACKE_dlasrt(id, (lapack_int) n, d);
}

_DLLAPI int __stdcall mql_dlaswp(int n, double *a, int lda, int k1, int k2, int *ipiv, int incx)
{
    return LAPACKE_dlaswp(LAPACK_COL_MAJOR, (lapack_int) n, a, (lapack_int) lda, (lapack_int) k1, (lapack_int) k2, ipiv, (lapack_int) incx);
}

#if defined TESTING
_DLLAPI int __stdcall mql_dlatms(int m, int n, char dist, int *iseed, char sym, double *d, int mode, double cond, double dmax, int kl, int ku, char pack, double *a, int lda)
{
    return LAPACKE_dlatms(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) n, dist, iseed, sym, d, (lapack_int) mode, cond, dmax, (lapack_int) kl, (lapack_int) ku, pack, a, (lapack_int) lda);
}
#endif

_DLLAPI int __stdcall mql_dlauum(char uplo, int n, double *a, int lda)
{
    return LAPACKE_dlauum(LAPACK_COL_MAJOR, uplo, (lapack_int) n, a, (lapack_int) lda);
}

_DLLAPI int __stdcall mql_dopgtr(char uplo, int n, double *ap, double *tau, double *q, int ldq)
{
    return LAPACKE_dopgtr(LAPACK_COL_MAJOR, uplo, (lapack_int) n, ap, tau, q, (lapack_int) ldq);
}

_DLLAPI int __stdcall mql_dopmtr(char side, char uplo, char trans, int m, int n, double *ap, double *tau, double *c, int ldc)
{
    return LAPACKE_dopmtr(LAPACK_COL_MAJOR, side, uplo, trans, (lapack_int) m, (lapack_int) n, ap, tau, c, (lapack_int) ldc);
}

_DLLAPI int __stdcall mql_dorgbr(char vect, int m, int n, int k, double *a, int lda, double *tau)
{
    return LAPACKE_dorgbr(LAPACK_COL_MAJOR, vect, (lapack_int) m, (lapack_int) n, (lapack_int) k, a, (lapack_int) lda, tau);
}

_DLLAPI int __stdcall mql_dorghr(int n, int ilo, int ihi, double *a, int lda, double *tau)
{
    return LAPACKE_dorghr(LAPACK_COL_MAJOR, (lapack_int) n, (lapack_int) ilo, (lapack_int) ihi, a, (lapack_int) lda, tau);
}

_DLLAPI int __stdcall mql_dorglq(int m, int n, int k, double *a, int lda, double *tau)
{
    return LAPACKE_dorglq(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) n, (lapack_int) k, a, (lapack_int) lda, tau);
}

_DLLAPI int __stdcall mql_dorgql(int m, int n, int k, double *a, int lda, double *tau)
{
    return LAPACKE_dorgql(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) n, (lapack_int) k, a, (lapack_int) lda, tau);
}

_DLLAPI int __stdcall mql_dorgqr(int m, int n, int k, double *a, int lda, double *tau)
{
    return LAPACKE_dorgqr(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) n, (lapack_int) k, a, (lapack_int) lda, tau);
}

_DLLAPI int __stdcall mql_dorgrq(int m, int n, int k, double *a, int lda, double *tau)
{
    return LAPACKE_dorgrq(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) n, (lapack_int) k, a, (lapack_int) lda, tau);
}

_DLLAPI int __stdcall mql_dorgtr(char uplo, int n, double *a, int lda, double *tau)
{
    return LAPACKE_dorgtr(LAPACK_COL_MAJOR, uplo, (lapack_int) n, a, (lapack_int) lda, tau);
}

_DLLAPI int __stdcall mql_dormbr(char vect, char side, char trans, int m, int n, int k, double *a, int lda, double *tau, double *c, int ldc)
{
    return LAPACKE_dormbr(LAPACK_COL_MAJOR, vect, side, trans, (lapack_int) m, (lapack_int) n, (lapack_int) k, a, (lapack_int) lda, tau, c, (lapack_int) ldc);
}

_DLLAPI int __stdcall mql_dormhr(char side, char trans, int m, int n, int ilo, int ihi, double *a, int lda, double *tau, double *c, int ldc)
{
    return LAPACKE_dormhr(LAPACK_COL_MAJOR, side, trans, (lapack_int) m, (lapack_int) n, (lapack_int) ilo, (lapack_int) ihi, a, (lapack_int) lda, tau, c, (lapack_int) ldc);
}

_DLLAPI int __stdcall mql_dormlq(char side, char trans, int m, int n, int k, double *a, int lda, double *tau, double *c, int ldc)
{
    return LAPACKE_dormlq(LAPACK_COL_MAJOR, side, trans, (lapack_int) m, (lapack_int) n, (lapack_int) k, a, (lapack_int) lda, tau, c, (lapack_int) ldc);
}

_DLLAPI int __stdcall mql_dormql(char side, char trans, int m, int n, int k, double *a, int lda, double *tau, double *c, int ldc)
{
    return LAPACKE_dormql(LAPACK_COL_MAJOR, side, trans, (lapack_int) m, (lapack_int) n, (lapack_int) k, a, (lapack_int) lda, tau, c, (lapack_int) ldc);
}

_DLLAPI int __stdcall mql_dormqr(char side, char trans, int m, int n, int k, double *a, int lda, double *tau, double *c, int ldc)
{
    return LAPACKE_dormqr(LAPACK_COL_MAJOR, side, trans, (lapack_int) m, (lapack_int) n, (lapack_int) k, a, (lapack_int) lda, tau, c, (lapack_int) ldc);
}

_DLLAPI int __stdcall mql_dormrq(char side, char trans, int m, int n, int k, double *a, int lda, double *tau, double *c, int ldc)
{
    return LAPACKE_dormrq(LAPACK_COL_MAJOR, side, trans, (lapack_int) m, (lapack_int) n, (lapack_int) k, a, (lapack_int) lda, tau, c, (lapack_int) ldc);
}

_DLLAPI int __stdcall mql_dormrz(char side, char trans, int m, int n, int k, int l, double *a, int lda, double *tau, double *c, int ldc)
{
    return LAPACKE_dormrz(LAPACK_COL_MAJOR, side, trans, (lapack_int) m, (lapack_int) n, (lapack_int) k, (lapack_int) l, a, (lapack_int) lda, tau, c, (lapack_int) ldc);
}

_DLLAPI int __stdcall mql_dormtr(char side, char uplo, char trans, int m, int n, double *a, int lda, double *tau, double *c, int ldc)
{
    return LAPACKE_dormtr(LAPACK_COL_MAJOR, side, uplo, trans, (lapack_int) m, (lapack_int) n, a, (lapack_int) lda, tau, c, (lapack_int) ldc);
}

_DLLAPI int __stdcall mql_dpbcon(char uplo, int n, int kd, double *ab, int ldab, double anorm, double *rcond)
{
    return LAPACKE_dpbcon(LAPACK_COL_MAJOR, uplo, (lapack_int) n, (lapack_int) kd, ab, (lapack_int) ldab, anorm, rcond);
}

_DLLAPI int __stdcall mql_dpbequ(char uplo, int n, int kd, double *ab, int ldab, double *s, double *scond, double *amax)
{
    return LAPACKE_dpbequ(LAPACK_COL_MAJOR, uplo, (lapack_int) n, (lapack_int) kd, ab, (lapack_int) ldab, s, scond, amax);
}

_DLLAPI int __stdcall mql_dpbrfs(char uplo, int n, int kd, int nrhs, double *ab, int ldab, double *afb, int ldafb, double *b, int ldb, double *x, int ldx, double *ferr, double *berr)
{
    return LAPACKE_dpbrfs(LAPACK_COL_MAJOR, uplo, (lapack_int) n, (lapack_int) kd, (lapack_int) nrhs, ab, (lapack_int) ldab, afb, (lapack_int) ldafb, b, (lapack_int) ldb, x, (lapack_int) ldx, ferr, berr);
}

_DLLAPI int __stdcall mql_dpbstf(char uplo, int n, int kb, double *bb, int ldbb)
{
    return LAPACKE_dpbstf(LAPACK_COL_MAJOR, uplo, (lapack_int) n, (lapack_int) kb, bb, (lapack_int) ldbb);
}

_DLLAPI int __stdcall mql_dpbsv(char uplo, int n, int kd, int nrhs, double *ab, int ldab, double *b, int ldb)
{
    return LAPACKE_dpbsv(LAPACK_COL_MAJOR, uplo, (lapack_int) n, (lapack_int) kd, (lapack_int) nrhs, ab, (lapack_int) ldab, b, (lapack_int) ldb);
}

_DLLAPI int __stdcall mql_dpbsvx(char fact, char uplo, int n, int kd, int nrhs, double *ab, int ldab, double *afb, int ldafb, char *equed, double *s, double *b, int ldb, double *x, int ldx, double *rcond, double *ferr, double *berr)
{
    return LAPACKE_dpbsvx(LAPACK_COL_MAJOR, fact, uplo, (lapack_int) n, (lapack_int) kd, (lapack_int) nrhs, ab, (lapack_int) ldab, afb, (lapack_int) ldafb, equed, s, b, (lapack_int) ldb, x, (lapack_int) ldx, rcond, ferr, berr);
}

_DLLAPI int __stdcall mql_dpbtrf(char uplo, int n, int kd, double *ab, int ldab)
{
    return LAPACKE_dpbtrf(LAPACK_COL_MAJOR, uplo, (lapack_int) n, (lapack_int) kd, ab, (lapack_int) ldab);
}

_DLLAPI int __stdcall mql_dpbtrs(char uplo, int n, int kd, int nrhs, double *ab, int ldab, double *b, int ldb)
{
    return LAPACKE_dpbtrs(LAPACK_COL_MAJOR, uplo, (lapack_int) n, (lapack_int) kd, (lapack_int) nrhs, ab, (lapack_int) ldab, b, (lapack_int) ldb);
}

_DLLAPI int __stdcall mql_dpftrf(char transr, char uplo, int n, double *a)
{
    return LAPACKE_dpftrf(LAPACK_COL_MAJOR, transr, uplo, (lapack_int) n, a);
}

_DLLAPI int __stdcall mql_dpftri(char transr, char uplo, int n, double *a)
{
    return LAPACKE_dpftri(LAPACK_COL_MAJOR, transr, uplo, (lapack_int) n, a);
}

_DLLAPI int __stdcall mql_dpftrs(char transr, char uplo, int n, int nrhs, double *a, double *b, int ldb)
{
    return LAPACKE_dpftrs(LAPACK_COL_MAJOR, transr, uplo, (lapack_int) n, (lapack_int) nrhs, a, b, (lapack_int) ldb);
}

_DLLAPI int __stdcall mql_dpocon(char uplo, int n, double *a, int lda, double anorm, double *rcond)
{
    return LAPACKE_dpocon(LAPACK_COL_MAJOR, uplo, (lapack_int) n, a, (lapack_int) lda, anorm, rcond);
}

_DLLAPI int __stdcall mql_dpoequ(int n, double *a, int lda, double *s, double *scond, double *amax)
{
    return LAPACKE_dpoequ(LAPACK_COL_MAJOR, (lapack_int) n, a, (lapack_int) lda, s, scond, amax);
}

_DLLAPI int __stdcall mql_dpoequb(int n, double *a, int lda, double *s, double *scond, double *amax)
{
    return LAPACKE_dpoequb(LAPACK_COL_MAJOR, (lapack_int) n, a, (lapack_int) lda, s, scond, amax);
}

_DLLAPI int __stdcall mql_dporfs(char uplo, int n, int nrhs, double *a, int lda, double *af, int ldaf, double *b, int ldb, double *x, int ldx, double *ferr, double *berr)
{
    return LAPACKE_dporfs(LAPACK_COL_MAJOR, uplo, (lapack_int) n, (lapack_int) nrhs, a, (lapack_int) lda, af, (lapack_int) ldaf, b, (lapack_int) ldb, x, (lapack_int) ldx, ferr, berr);
}

#if defined USEXBLAS
_DLLAPI int __stdcall mql_dporfsx(char uplo, char equed, int n, int nrhs, double *a, int lda, double *af, int ldaf, double *s, double *b, int ldb, double *x, int ldx, double *rcond, double *berr, int n_err_bnds, double *err_bnds_norm, double *err_bnds_comp, int nparams, double *params)
{
    return LAPACKE_dporfsx(LAPACK_COL_MAJOR, uplo, equed, (lapack_int) n, (lapack_int) nrhs, a, (lapack_int) lda, af, (lapack_int) ldaf, s, b, (lapack_int) ldb, x, (lapack_int) ldx, rcond, berr, (lapack_int) n_err_bnds, err_bnds_norm, err_bnds_comp, (lapack_int) nparams, params);
}
#endif

_DLLAPI int __stdcall mql_dposv(char uplo, int n, int nrhs, double *a, int lda, double *b, int ldb)
{
    return LAPACKE_dposv(LAPACK_COL_MAJOR, uplo, (lapack_int) n, (lapack_int) nrhs, a, (lapack_int) lda, b, (lapack_int) ldb);
}

_DLLAPI int __stdcall mql_dsposv(char uplo, int n, int nrhs, double *a, int lda, double *b, int ldb, double *x, int ldx, int *iter)
{
    return LAPACKE_dsposv(LAPACK_COL_MAJOR, uplo, (lapack_int) n, (lapack_int) nrhs, a, (lapack_int) lda, b, (lapack_int) ldb, x, (lapack_int) ldx, iter);
}

_DLLAPI int __stdcall mql_dposvx(char fact, char uplo, int n, int nrhs, double *a, int lda, double *af, int ldaf, char *equed, double *s, double *b, int ldb, double *x, int ldx, double *rcond, double *ferr, double *berr)
{
    return LAPACKE_dposvx(LAPACK_COL_MAJOR, fact, uplo, (lapack_int) n, (lapack_int) nrhs, a, (lapack_int) lda, af, (lapack_int) ldaf, equed, s, b, (lapack_int) ldb, x, (lapack_int) ldx, rcond, ferr, berr);
}

#if defined USEXBLAS
_DLLAPI int __stdcall mql_dposvxx(char fact, char uplo, int n, int nrhs, double *a, int lda, double *af, int ldaf, char *equed, double *s, double *b, int ldb, double *x, int ldx, double *rcond, double *rpvgrw, double *berr, int n_err_bnds, double *err_bnds_norm, double *err_bnds_comp, int nparams, double *params)
{
    return LAPACKE_dposvxx(LAPACK_COL_MAJOR, fact, uplo, (lapack_int) n, (lapack_int) nrhs, a, (lapack_int) lda, af, (lapack_int) ldaf, equed, s, b, (lapack_int) ldb, x, (lapack_int) ldx, rcond, rpvgrw, berr, (lapack_int) n_err_bnds, err_bnds_norm, err_bnds_comp, (lapack_int) nparams, params);
}
#endif

_DLLAPI int __stdcall mql_dpotrf(char uplo, int n, double *a, int lda)
{
    return LAPACKE_dpotrf(LAPACK_COL_MAJOR, uplo, (lapack_int) n, a, (lapack_int) lda);
}

_DLLAPI int __stdcall mql_dpotri(char uplo, int n, double *a, int lda)
{
    return LAPACKE_dpotri(LAPACK_COL_MAJOR, uplo, (lapack_int) n, a, (lapack_int) lda);
}

_DLLAPI int __stdcall mql_dpotrs(char uplo, int n, int nrhs, double *a, int lda, double *b, int ldb)
{
    return LAPACKE_dpotrs(LAPACK_COL_MAJOR, uplo, (lapack_int) n, (lapack_int) nrhs, a, (lapack_int) lda, b, (lapack_int) ldb);
}

_DLLAPI int __stdcall mql_dppcon(char uplo, int n, double *ap, double anorm, double *rcond)
{
    return LAPACKE_dppcon(LAPACK_COL_MAJOR, uplo, (lapack_int) n, ap, anorm, rcond);
}

_DLLAPI int __stdcall mql_dppequ(char uplo, int n, double *ap, double *s, double *scond, double *amax)
{
    return LAPACKE_dppequ(LAPACK_COL_MAJOR, uplo, (lapack_int) n, ap, s, scond, amax);
}

_DLLAPI int __stdcall mql_dpprfs(char uplo, int n, int nrhs, double *ap, double *afp, double *b, int ldb, double *x, int ldx, double *ferr, double *berr)
{
    return LAPACKE_dpprfs(LAPACK_COL_MAJOR, uplo, (lapack_int) n, (lapack_int) nrhs, ap, afp, b, (lapack_int) ldb, x, (lapack_int) ldx, ferr, berr);
}

_DLLAPI int __stdcall mql_dppsv(char uplo, int n, int nrhs, double *ap, double *b, int ldb)
{
    return LAPACKE_dppsv(LAPACK_COL_MAJOR, uplo, (lapack_int) n, (lapack_int) nrhs, ap, b, (lapack_int) ldb);
}

_DLLAPI int __stdcall mql_dppsvx(char fact, char uplo, int n, int nrhs, double *ap, double *afp, char *equed, double *s, double *b, int ldb, double *x, int ldx, double *rcond, double *ferr, double *berr)
{
    return LAPACKE_dppsvx(LAPACK_COL_MAJOR, fact, uplo, (lapack_int) n, (lapack_int) nrhs, ap, afp, equed, s, b, (lapack_int) ldb, x, (lapack_int) ldx, rcond, ferr, berr);
}

_DLLAPI int __stdcall mql_dpptrf(char uplo, int n, double *ap)
{
    return LAPACKE_dpptrf(LAPACK_COL_MAJOR, uplo, (lapack_int) n, ap);
}

_DLLAPI int __stdcall mql_dpptri(char uplo, int n, double *ap)
{
    return LAPACKE_dpptri(LAPACK_COL_MAJOR, uplo, (lapack_int) n, ap);
}

_DLLAPI int __stdcall mql_dpptrs(char uplo, int n, int nrhs, double *ap, double *b, int ldb)
{
    return LAPACKE_dpptrs(LAPACK_COL_MAJOR, uplo, (lapack_int) n, (lapack_int) nrhs, ap, b, (lapack_int) ldb);
}

_DLLAPI int __stdcall mql_dpstrf(char uplo, int n, double *a, int lda, int *piv, int *rank, double tol)
{
    return LAPACKE_dpstrf(LAPACK_COL_MAJOR, uplo, (lapack_int) n, a, (lapack_int) lda, piv, rank, tol);
}

_DLLAPI int __stdcall mql_dptcon(int n, double *d, double *e, double anorm, double *rcond)
{
    return LAPACKE_dptcon((lapack_int) n, d, e, anorm, rcond);
}

_DLLAPI int __stdcall mql_dpteqr(char compz, int n, double *d, double *e, double *z, int ldz)
{
    return LAPACKE_dpteqr(LAPACK_COL_MAJOR, compz, (lapack_int) n, d, e, z, (lapack_int) ldz);
}

_DLLAPI int __stdcall mql_dptrfs(int n, int nrhs, double *d, double *e, double *df, double *ef, double *b, int ldb, double *x, int ldx, double *ferr, double *berr)
{
    return LAPACKE_dptrfs(LAPACK_COL_MAJOR, (lapack_int) n, (lapack_int) nrhs, d, e, df, ef, b, (lapack_int) ldb, x, (lapack_int) ldx, ferr, berr);
}

_DLLAPI int __stdcall mql_dptsv(int n, int nrhs, double *d, double *e, double *b, int ldb)
{
    return LAPACKE_dptsv(LAPACK_COL_MAJOR, (lapack_int) n, (lapack_int) nrhs, d, e, b, (lapack_int) ldb);
}

_DLLAPI int __stdcall mql_dptsvx(char fact, int n, int nrhs, double *d, double *e, double *df, double *ef, double *b, int ldb, double *x, int ldx, double *rcond, double *ferr, double *berr)
{
    return LAPACKE_dptsvx(LAPACK_COL_MAJOR, fact, (lapack_int) n, (lapack_int) nrhs, d, e, df, ef, b, (lapack_int) ldb, x, (lapack_int) ldx, rcond, ferr, berr);
}

_DLLAPI int __stdcall mql_dpttrf(int n, double *d, double *e)
{
    return LAPACKE_dpttrf((lapack_int) n, d, e);
}

_DLLAPI int __stdcall mql_dpttrs(int n, int nrhs, double *d, double *e, double *b, int ldb)
{
    return LAPACKE_dpttrs(LAPACK_COL_MAJOR, (lapack_int) n, (lapack_int) nrhs, d, e, b, (lapack_int) ldb);
}

_DLLAPI int __stdcall mql_dsbev(char jobz, char uplo, int n, int kd, double *ab, int ldab, double *w, double *z, int ldz)
{
    return LAPACKE_dsbev(LAPACK_COL_MAJOR, jobz, uplo, (lapack_int) n, (lapack_int) kd, ab, (lapack_int) ldab, w, z, (lapack_int) ldz);
}

_DLLAPI int __stdcall mql_dsbevd(char jobz, char uplo, int n, int kd, double *ab, int ldab, double *w, double *z, int ldz)
{
    return LAPACKE_dsbevd(LAPACK_COL_MAJOR, jobz, uplo, (lapack_int) n, (lapack_int) kd, ab, (lapack_int) ldab, w, z, (lapack_int) ldz);
}

_DLLAPI int __stdcall mql_dsbevx(char jobz, char range, char uplo, int n, int kd, double *ab, int ldab, double *q, int ldq, double vl, double vu, int il, int iu, double abstol, int *m, double *w, double *z, int ldz, int *ifail)
{
    return LAPACKE_dsbevx(LAPACK_COL_MAJOR, jobz, range, uplo, (lapack_int) n, (lapack_int) kd, ab, (lapack_int) ldab, q, (lapack_int) ldq, vl, vu, (lapack_int) il, (lapack_int) iu, abstol, m, w, z, (lapack_int) ldz, ifail);
}

_DLLAPI int __stdcall mql_dsbgst(char vect, char uplo, int n, int ka, int kb, double *ab, int ldab, double *bb, int ldbb, double *x, int ldx)
{
    return LAPACKE_dsbgst(LAPACK_COL_MAJOR, vect, uplo, (lapack_int) n, (lapack_int) ka, (lapack_int) kb, ab, (lapack_int) ldab, bb, (lapack_int) ldbb, x, (lapack_int) ldx);
}

_DLLAPI int __stdcall mql_dsbgv(char jobz, char uplo, int n, int ka, int kb, double *ab, int ldab, double *bb, int ldbb, double *w, double *z, int ldz)
{
    return LAPACKE_dsbgv(LAPACK_COL_MAJOR, jobz, uplo, (lapack_int) n, (lapack_int) ka, (lapack_int) kb, ab, (lapack_int) ldab, bb, (lapack_int) ldbb, w, z, (lapack_int) ldz);
}

_DLLAPI int __stdcall mql_dsbgvd(char jobz, char uplo, int n, int ka, int kb, double *ab, int ldab, double *bb, int ldbb, double *w, double *z, int ldz)
{
    return LAPACKE_dsbgvd(LAPACK_COL_MAJOR, jobz, uplo, (lapack_int) n, (lapack_int) ka, (lapack_int) kb, ab, (lapack_int) ldab, bb, (lapack_int) ldbb, w, z, (lapack_int) ldz);
}

_DLLAPI int __stdcall mql_dsbgvx(char jobz, char range, char uplo, int n, int ka, int kb, double *ab, int ldab, double *bb, int ldbb, double *q, int ldq, double vl, double vu, int il, int iu, double abstol, int *m, double *w, double *z, int ldz, int *ifail)
{
    return LAPACKE_dsbgvx(LAPACK_COL_MAJOR, jobz, range, uplo, (lapack_int) n, (lapack_int) ka, (lapack_int) kb, ab, (lapack_int) ldab, bb, (lapack_int) ldbb, q, (lapack_int) ldq, vl, vu, (lapack_int) il, (lapack_int) iu, abstol, m, w, z, (lapack_int) ldz, ifail);
}

_DLLAPI int __stdcall mql_dsbtrd(char vect, char uplo, int n, int kd, double *ab, int ldab, double *d, double *e, double *q, int ldq)
{
    return LAPACKE_dsbtrd(LAPACK_COL_MAJOR, vect, uplo, (lapack_int) n, (lapack_int) kd, ab, (lapack_int) ldab, d, e, q, (lapack_int) ldq);
}

_DLLAPI int __stdcall mql_dsfrk(char transr, char uplo, char trans, int n, int k, double alpha, double *a, int lda, double beta, double *c)
{
    return LAPACKE_dsfrk(LAPACK_COL_MAJOR, transr, uplo, trans, (lapack_int) n, (lapack_int) k, alpha, a, (lapack_int) lda, beta, c);
}

_DLLAPI int __stdcall mql_dspcon(char uplo, int n, double *ap, int *ipiv, double anorm, double *rcond)
{
    return LAPACKE_dspcon(LAPACK_COL_MAJOR, uplo, (lapack_int) n, ap, ipiv, anorm, rcond);
}

_DLLAPI int __stdcall mql_dspev(char jobz, char uplo, int n, double *ap, double *w, double *z, int ldz)
{
    return LAPACKE_dspev(LAPACK_COL_MAJOR, jobz, uplo, (lapack_int) n, ap, w, z, (lapack_int) ldz);
}

_DLLAPI int __stdcall mql_dspevd(char jobz, char uplo, int n, double *ap, double *w, double *z, int ldz)
{
    return LAPACKE_dspevd(LAPACK_COL_MAJOR, jobz, uplo, (lapack_int) n, ap, w, z, (lapack_int) ldz);
}

_DLLAPI int __stdcall mql_dspevx(char jobz, char range, char uplo, int n, double *ap, double vl, double vu, int il, int iu, double abstol, int *m, double *w, double *z, int ldz, int *ifail)
{
    return LAPACKE_dspevx(LAPACK_COL_MAJOR, jobz, range, uplo, (lapack_int) n, ap, vl, vu, (lapack_int) il, (lapack_int) iu, abstol, m, w, z, (lapack_int) ldz, ifail);
}

_DLLAPI int __stdcall mql_dspgst(int itype, char uplo, int n, double *ap, double *bp)
{
    return LAPACKE_dspgst(LAPACK_COL_MAJOR, (lapack_int) itype, uplo, (lapack_int) n, ap, bp);
}

_DLLAPI int __stdcall mql_dspgv(int itype, char jobz, char uplo, int n, double *ap, double *bp, double *w, double *z, int ldz)
{
    return LAPACKE_dspgv(LAPACK_COL_MAJOR, (lapack_int) itype, jobz, uplo, (lapack_int) n, ap, bp, w, z, (lapack_int) ldz);
}

_DLLAPI int __stdcall mql_dspgvd(int itype, char jobz, char uplo, int n, double *ap, double *bp, double *w, double *z, int ldz)
{
    return LAPACKE_dspgvd(LAPACK_COL_MAJOR, (lapack_int) itype, jobz, uplo, (lapack_int) n, ap, bp, w, z, (lapack_int) ldz);
}

_DLLAPI int __stdcall mql_dspgvx(int itype, char jobz, char range, char uplo, int n, double *ap, double *bp, double vl, double vu, int il, int iu, double abstol, int *m, double *w, double *z, int ldz, int *ifail)
{
    return LAPACKE_dspgvx(LAPACK_COL_MAJOR, (lapack_int) itype, jobz, range, uplo, (lapack_int) n, ap, bp, vl, vu, (lapack_int) il, (lapack_int) iu, abstol, m, w, z, (lapack_int) ldz, ifail);
}

_DLLAPI int __stdcall mql_dsprfs(char uplo, int n, int nrhs, double *ap, double *afp, int *ipiv, double *b, int ldb, double *x, int ldx, double *ferr, double *berr)
{
    return LAPACKE_dsprfs(LAPACK_COL_MAJOR, uplo, (lapack_int) n, (lapack_int) nrhs, ap, afp, ipiv, b, (lapack_int) ldb, x, (lapack_int) ldx, ferr, berr);
}

_DLLAPI int __stdcall mql_dspsv(char uplo, int n, int nrhs, double *ap, int *ipiv, double *b, int ldb)
{
    return LAPACKE_dspsv(LAPACK_COL_MAJOR, uplo, (lapack_int) n, (lapack_int) nrhs, ap, ipiv, b, (lapack_int) ldb);
}

_DLLAPI int __stdcall mql_dspsvx(char fact, char uplo, int n, int nrhs, double *ap, double *afp, int *ipiv, double *b, int ldb, double *x, int ldx, double *rcond, double *ferr, double *berr)
{
    return LAPACKE_dspsvx(LAPACK_COL_MAJOR, fact, uplo, (lapack_int) n, (lapack_int) nrhs, ap, afp, ipiv, b, (lapack_int) ldb, x, (lapack_int) ldx, rcond, ferr, berr);
}

_DLLAPI int __stdcall mql_dsptrd(char uplo, int n, double *ap, double *d, double *e, double *tau)
{
    return LAPACKE_dsptrd(LAPACK_COL_MAJOR, uplo, (lapack_int) n, ap, d, e, tau);
}

_DLLAPI int __stdcall mql_dsptrf(char uplo, int n, double *ap, int *ipiv)
{
    return LAPACKE_dsptrf(LAPACK_COL_MAJOR, uplo, (lapack_int) n, ap, ipiv);
}

_DLLAPI int __stdcall mql_dsptri(char uplo, int n, double *ap, int *ipiv)
{
    return LAPACKE_dsptri(LAPACK_COL_MAJOR, uplo, (lapack_int) n, ap, ipiv);
}

_DLLAPI int __stdcall mql_dsptrs(char uplo, int n, int nrhs, double *ap, int *ipiv, double *b, int ldb)
{
    return LAPACKE_dsptrs(LAPACK_COL_MAJOR, uplo, (lapack_int) n, (lapack_int) nrhs, ap, ipiv, b, (lapack_int) ldb);
}

_DLLAPI int __stdcall mql_dstebz(char range, char order, int n, double vl, double vu, int il, int iu, double abstol, double *d, double *e, int *m, int *nsplit, double *w, int *iblock, int *isplit)
{
    return LAPACKE_dstebz(range, order, (lapack_int) n, vl, vu, (lapack_int) il, (lapack_int) iu, abstol, d, e, m, nsplit, w, iblock, isplit);
}

_DLLAPI int __stdcall mql_dstedc(char compz, int n, double *d, double *e, double *z, int ldz)
{
    return LAPACKE_dstedc(LAPACK_COL_MAJOR, compz, (lapack_int) n, d, e, z, (lapack_int) ldz);
}

_DLLAPI int __stdcall mql_dstegr(char jobz, char range, int n, double *d, double *e, double vl, double vu, int il, int iu, double abstol, int *m, double *w, double *z, int ldz, int *isuppz)
{
    return LAPACKE_dstegr(LAPACK_COL_MAJOR, jobz, range, (lapack_int) n, d, e, vl, vu, (lapack_int) il, (lapack_int) iu, abstol, m, w, z, (lapack_int) ldz, isuppz);
}

_DLLAPI int __stdcall mql_dstein(int n, double *d, double *e, int m, double *w, int *iblock, int *isplit, double *z, int ldz, int *ifailv)
{
    return LAPACKE_dstein(LAPACK_COL_MAJOR, (lapack_int) n, d, e, (lapack_int) m, w, iblock, isplit, z, (lapack_int) ldz, ifailv);
}

_DLLAPI int __stdcall mql_dsteqr(char compz, int n, double *d, double *e, double *z, int ldz)
{
    return LAPACKE_dsteqr(LAPACK_COL_MAJOR, compz, (lapack_int) n, d, e, z, (lapack_int) ldz);
}

_DLLAPI int __stdcall mql_dsterf(int n, double *d, double *e)
{
    return LAPACKE_dsterf((lapack_int) n, d, e);
}

_DLLAPI int __stdcall mql_dstev(char jobz, int n, double *d, double *e, double *z, int ldz)
{
    return LAPACKE_dstev(LAPACK_COL_MAJOR, jobz, (lapack_int) n, d, e, z, (lapack_int) ldz);
}

_DLLAPI int __stdcall mql_dstevd(char jobz, int n, double *d, double *e, double *z, int ldz)
{
    return LAPACKE_dstevd(LAPACK_COL_MAJOR, jobz, (lapack_int) n, d, e, z, (lapack_int) ldz);
}

_DLLAPI int __stdcall mql_dstevr(char jobz, char range, int n, double *d, double *e, double vl, double vu, int il, int iu, double abstol, int *m, double *w, double *z, int ldz, int *isuppz)
{
    return LAPACKE_dstevr(LAPACK_COL_MAJOR, jobz, range, (lapack_int) n, d, e, vl, vu, (lapack_int) il, (lapack_int) iu, abstol, m, w, z, (lapack_int) ldz, isuppz);
}

_DLLAPI int __stdcall mql_dstevx(char jobz, char range, int n, double *d, double *e, double vl, double vu, int il, int iu, double abstol, int *m, double *w, double *z, int ldz, int *ifail)
{
    return LAPACKE_dstevx(LAPACK_COL_MAJOR, jobz, range, (lapack_int) n, d, e, vl, vu, (lapack_int) il, (lapack_int) iu, abstol, m, w, z, (lapack_int) ldz, ifail);
}

_DLLAPI int __stdcall mql_dsycon(char uplo, int n, double *a, int lda, int *ipiv, double anorm, double *rcond)
{
    return LAPACKE_dsycon(LAPACK_COL_MAJOR, uplo, (lapack_int) n, a, (lapack_int) lda, ipiv, anorm, rcond);
}

_DLLAPI int __stdcall mql_dsyequb(char uplo, int n, double *a, int lda, double *s, double *scond, double *amax)
{
    return LAPACKE_dsyequb(LAPACK_COL_MAJOR, uplo, (lapack_int) n, a, (lapack_int) lda, s, scond, amax);
}

_DLLAPI int __stdcall mql_dsyev(char jobz, char uplo, int n, double *a, int lda, double *w)
{
    return LAPACKE_dsyev(LAPACK_COL_MAJOR, jobz, uplo, (lapack_int) n, a, (lapack_int) lda, w);
}

_DLLAPI int __stdcall mql_dsyevd(char jobz, char uplo, int n, double *a, int lda, double *w)
{
    return LAPACKE_dsyevd(LAPACK_COL_MAJOR, jobz, uplo, (lapack_int) n, a, (lapack_int) lda, w);
}

_DLLAPI int __stdcall mql_dsyevr(char jobz, char range, char uplo, int n, double *a, int lda, double vl, double vu, int il, int iu, double abstol, int *m, double *w, double *z, int ldz, int *isuppz)
{
    return LAPACKE_dsyevr(LAPACK_COL_MAJOR, jobz, range, uplo, (lapack_int) n, a, (lapack_int) lda, vl, vu, (lapack_int) il, (lapack_int) iu, abstol, m, w, z, (lapack_int) ldz, isuppz);
}

_DLLAPI int __stdcall mql_dsyevx(char jobz, char range, char uplo, int n, double *a, int lda, double vl, double vu, int il, int iu, double abstol, int *m, double *w, double *z, int ldz, int *ifail)
{
    return LAPACKE_dsyevx(LAPACK_COL_MAJOR, jobz, range, uplo, (lapack_int) n, a, (lapack_int) lda, vl, vu, (lapack_int) il, (lapack_int) iu, abstol, m, w, z, (lapack_int) ldz, ifail);
}

_DLLAPI int __stdcall mql_dsygst(int itype, char uplo, int n, double *a, int lda, double *b, int ldb)
{
    return LAPACKE_dsygst(LAPACK_COL_MAJOR, (lapack_int) itype, uplo, (lapack_int) n, a, (lapack_int) lda, b, (lapack_int) ldb);
}

_DLLAPI int __stdcall mql_dsygv(int itype, char jobz, char uplo, int n, double *a, int lda, double *b, int ldb, double *w)
{
    return LAPACKE_dsygv(LAPACK_COL_MAJOR, (lapack_int) itype, jobz, uplo, (lapack_int) n, a, (lapack_int) lda, b, (lapack_int) ldb, w);
}

_DLLAPI int __stdcall mql_dsygvd(int itype, char jobz, char uplo, int n, double *a, int lda, double *b, int ldb, double *w)
{
    return LAPACKE_dsygvd(LAPACK_COL_MAJOR, (lapack_int) itype, jobz, uplo, (lapack_int) n, a, (lapack_int) lda, b, (lapack_int) ldb, w);
}

_DLLAPI int __stdcall mql_dsygvx(int itype, char jobz, char range, char uplo, int n, double *a, int lda, double *b, int ldb, double vl, double vu, int il, int iu, double abstol, int *m, double *w, double *z, int ldz, int *ifail)
{
    return LAPACKE_dsygvx(LAPACK_COL_MAJOR, (lapack_int) itype, jobz, range, uplo, (lapack_int) n, a, (lapack_int) lda, b, (lapack_int) ldb, vl, vu, (lapack_int) il, (lapack_int) iu, abstol, m, w, z, (lapack_int) ldz, ifail);
}

_DLLAPI int __stdcall mql_dsyrfs(char uplo, int n, int nrhs, double *a, int lda, double *af, int ldaf, int *ipiv, double *b, int ldb, double *x, int ldx, double *ferr, double *berr)
{
    return LAPACKE_dsyrfs(LAPACK_COL_MAJOR, uplo, (lapack_int) n, (lapack_int) nrhs, a, (lapack_int) lda, af, (lapack_int) ldaf, ipiv, b, (lapack_int) ldb, x, (lapack_int) ldx, ferr, berr);
}

#if defined USEXBLAS
_DLLAPI int __stdcall mql_dsyrfsx(char uplo, char equed, int n, int nrhs, double *a, int lda, double *af, int ldaf, int *ipiv, double *s, double *b, int ldb, double *x, int ldx, double *rcond, double *berr, int n_err_bnds, double *err_bnds_norm, double *err_bnds_comp, int nparams, double *params)
{
    return LAPACKE_dsyrfsx(LAPACK_COL_MAJOR, uplo, equed, (lapack_int) n, (lapack_int) nrhs, a, (lapack_int) lda, af, (lapack_int) ldaf, ipiv, s, b, (lapack_int) ldb, x, (lapack_int) ldx, rcond, berr, (lapack_int) n_err_bnds, err_bnds_norm, err_bnds_comp, (lapack_int) nparams, params);
}
#endif

_DLLAPI int __stdcall mql_dsysv(char uplo, int n, int nrhs, double *a, int lda, int *ipiv, double *b, int ldb)
{
    return LAPACKE_dsysv(LAPACK_COL_MAJOR, uplo, (lapack_int) n, (lapack_int) nrhs, a, (lapack_int) lda, ipiv, b, (lapack_int) ldb);
}

_DLLAPI int __stdcall mql_dsysvx(char fact, char uplo, int n, int nrhs, double *a, int lda, double *af, int ldaf, int *ipiv, double *b, int ldb, double *x, int ldx, double *rcond, double *ferr, double *berr)
{
    return LAPACKE_dsysvx(LAPACK_COL_MAJOR, fact, uplo, (lapack_int) n, (lapack_int) nrhs, a, (lapack_int) lda, af, (lapack_int) ldaf, ipiv, b, (lapack_int) ldb, x, (lapack_int) ldx, rcond, ferr, berr);
}

#if defined USEXBLAS
_DLLAPI int __stdcall mql_dsysvxx(char fact, char uplo, int n, int nrhs, double *a, int lda, double *af, int ldaf, int *ipiv, char *equed, double *s, double *b, int ldb, double *x, int ldx, double *rcond, double *rpvgrw, double *berr, int n_err_bnds, double *err_bnds_norm, double *err_bnds_comp, int nparams, double *params)
{
    return LAPACKE_dsysvxx(LAPACK_COL_MAJOR, fact, uplo, (lapack_int) n, (lapack_int) nrhs, a, (lapack_int) lda, af, (lapack_int) ldaf, ipiv, equed, s, b, (lapack_int) ldb, x, (lapack_int) ldx, rcond, rpvgrw, berr, (lapack_int) n_err_bnds, err_bnds_norm, err_bnds_comp, (lapack_int) nparams, params);
}
#endif

_DLLAPI int __stdcall mql_dsytrd(char uplo, int n, double *a, int lda, double *d, double *e, double *tau)
{
    return LAPACKE_dsytrd(LAPACK_COL_MAJOR, uplo, (lapack_int) n, a, (lapack_int) lda, d, e, tau);
}

_DLLAPI int __stdcall mql_dsytrf(char uplo, int n, double *a, int lda, int *ipiv)
{
    return LAPACKE_dsytrf(LAPACK_COL_MAJOR, uplo, (lapack_int) n, a, (lapack_int) lda, ipiv);
}

_DLLAPI int __stdcall mql_dsytri(char uplo, int n, double *a, int lda, int *ipiv)
{
    return LAPACKE_dsytri(LAPACK_COL_MAJOR, uplo, (lapack_int) n, a, (lapack_int) lda, ipiv);
}

_DLLAPI int __stdcall mql_dsytrs(char uplo, int n, int nrhs, double *a, int lda, int *ipiv, double *b, int ldb)
{
    return LAPACKE_dsytrs(LAPACK_COL_MAJOR, uplo, (lapack_int) n, (lapack_int) nrhs, a, (lapack_int) lda, ipiv, b, (lapack_int) ldb);
}

_DLLAPI int __stdcall mql_dtbcon(char norm, char uplo, char diag, int n, int kd, double *ab, int ldab, double *rcond)
{
    return LAPACKE_dtbcon(LAPACK_COL_MAJOR, norm, uplo, diag, (lapack_int) n, (lapack_int) kd, ab, (lapack_int) ldab, rcond);
}

_DLLAPI int __stdcall mql_dtbrfs(char uplo, char trans, char diag, int n, int kd, int nrhs, double *ab, int ldab, double *b, int ldb, double *x, int ldx, double *ferr, double *berr)
{
    return LAPACKE_dtbrfs(LAPACK_COL_MAJOR, uplo, trans, diag, (lapack_int) n, (lapack_int) kd, (lapack_int) nrhs, ab, (lapack_int) ldab, b, (lapack_int) ldb, x, (lapack_int) ldx, ferr, berr);
}

_DLLAPI int __stdcall mql_dtbtrs(char uplo, char trans, char diag, int n, int kd, int nrhs, double *ab, int ldab, double *b, int ldb)
{
    return LAPACKE_dtbtrs(LAPACK_COL_MAJOR, uplo, trans, diag, (lapack_int) n, (lapack_int) kd, (lapack_int) nrhs, ab, (lapack_int) ldab, b, (lapack_int) ldb);
}

_DLLAPI int __stdcall mql_dtfsm(char transr, char side, char uplo, char trans, char diag, int m, int n, double alpha, double *a, double *b, int ldb)
{
    return LAPACKE_dtfsm(LAPACK_COL_MAJOR, transr, side, uplo, trans, diag, (lapack_int) m, (lapack_int) n, alpha, a, b, (lapack_int) ldb);
}

_DLLAPI int __stdcall mql_dtftri(char transr, char uplo, char diag, int n, double *a)
{
    return LAPACKE_dtftri(LAPACK_COL_MAJOR, transr, uplo, diag, (lapack_int) n, a);
}

_DLLAPI int __stdcall mql_dtfttp(char transr, char uplo, int n, double *arf, double *ap)
{
    return LAPACKE_dtfttp(LAPACK_COL_MAJOR, transr, uplo, (lapack_int) n, arf, ap);
}

_DLLAPI int __stdcall mql_dtfttr(char transr, char uplo, int n, double *arf, double *a, int lda)
{
    return LAPACKE_dtfttr(LAPACK_COL_MAJOR, transr, uplo, (lapack_int) n, arf, a, (lapack_int) lda);
}

_DLLAPI int __stdcall mql_dtgsja(char jobu, char jobv, char jobq, int m, int p, int n, int k, int l, double *a, int lda, double *b, int ldb, double tola, double tolb, double *alpha, double *beta, double *u, int ldu, double *v, int ldv, double *q, int ldq, int *ncycle)
{
    return LAPACKE_dtgsja(LAPACK_COL_MAJOR, jobu, jobv, jobq, (lapack_int) m, (lapack_int) p, (lapack_int) n, (lapack_int) k, (lapack_int) l, a, (lapack_int) lda, b, (lapack_int) ldb, tola, tolb, alpha, beta, u, (lapack_int) ldu, v, (lapack_int) ldv, q, (lapack_int) ldq, ncycle);
}

_DLLAPI int __stdcall mql_dtgsyl(char trans, int ijob, int m, int n, double *a, int lda, double *b, int ldb, double *c, int ldc, double *d, int ldd, double *e, int lde, double *f, int ldf, double *scale, double *dif)
{
    return LAPACKE_dtgsyl(LAPACK_COL_MAJOR, trans, (lapack_int) ijob, (lapack_int) m, (lapack_int) n, a, (lapack_int) lda, b, (lapack_int) ldb, c, (lapack_int) ldc, d, (lapack_int) ldd, e, (lapack_int) lde, f, (lapack_int) ldf, scale, dif);
}

_DLLAPI int __stdcall mql_dtpcon(char norm, char uplo, char diag, int n, double *ap, double *rcond)
{
    return LAPACKE_dtpcon(LAPACK_COL_MAJOR, norm, uplo, diag, (lapack_int) n, ap, rcond);
}

_DLLAPI int __stdcall mql_dtprfs(char uplo, char trans, char diag, int n, int nrhs, double *ap, double *b, int ldb, double *x, int ldx, double *ferr, double *berr)
{
    return LAPACKE_dtprfs(LAPACK_COL_MAJOR, uplo, trans, diag, (lapack_int) n, (lapack_int) nrhs, ap, b, (lapack_int) ldb, x, (lapack_int) ldx, ferr, berr);
}

_DLLAPI int __stdcall mql_dtptri(char uplo, char diag, int n, double *ap)
{
    return LAPACKE_dtptri(LAPACK_COL_MAJOR, uplo, diag, (lapack_int) n, ap);
}

_DLLAPI int __stdcall mql_dtptrs(char uplo, char trans, char diag, int n, int nrhs, double *ap, double *b, int ldb)
{
    return LAPACKE_dtptrs(LAPACK_COL_MAJOR, uplo, trans, diag, (lapack_int) n, (lapack_int) nrhs, ap, b, (lapack_int) ldb);
}

_DLLAPI int __stdcall mql_dtpttf(char transr, char uplo, int n, double *ap, double *arf)
{
    return LAPACKE_dtpttf(LAPACK_COL_MAJOR, transr, uplo, (lapack_int) n, ap, arf);
}

_DLLAPI int __stdcall mql_dtpttr(char uplo, int n, double *ap, double *a, int lda)
{
    return LAPACKE_dtpttr(LAPACK_COL_MAJOR, uplo, (lapack_int) n, ap, a, (lapack_int) lda);
}

_DLLAPI int __stdcall mql_dtrcon(char norm, char uplo, char diag, int n, double *a, int lda, double *rcond)
{
    return LAPACKE_dtrcon(LAPACK_COL_MAJOR, norm, uplo, diag, (lapack_int) n, a, (lapack_int) lda, rcond);
}

_DLLAPI int __stdcall mql_dtrexc(char compq, int n, double *t, int ldt, double *q, int ldq, int *ifst, int *ilst)
{
    return LAPACKE_dtrexc(LAPACK_COL_MAJOR, compq, (lapack_int) n, t, (lapack_int) ldt, q, (lapack_int) ldq, ifst, ilst);
}

_DLLAPI int __stdcall mql_dtrrfs(char uplo, char trans, char diag, int n, int nrhs, double *a, int lda, double *b, int ldb, double *x, int ldx, double *ferr, double *berr)
{
    return LAPACKE_dtrrfs(LAPACK_COL_MAJOR, uplo, trans, diag, (lapack_int) n, (lapack_int) nrhs, a, (lapack_int) lda, b, (lapack_int) ldb, x, (lapack_int) ldx, ferr, berr);
}

_DLLAPI int __stdcall mql_dtrsyl(char trana, char tranb, int isgn, int m, int n, double *a, int lda, double *b, int ldb, double *c, int ldc, double *scale)
{
    return LAPACKE_dtrsyl(LAPACK_COL_MAJOR, trana, tranb, (lapack_int) isgn, (lapack_int) m, (lapack_int) n, a, (lapack_int) lda, b, (lapack_int) ldb, c, (lapack_int) ldc, scale);
}

_DLLAPI int __stdcall mql_dtrtri(char uplo, char diag, int n, double *a, int lda)
{
    return LAPACKE_dtrtri(LAPACK_COL_MAJOR, uplo, diag, (lapack_int) n, a, (lapack_int) lda);
}

_DLLAPI int __stdcall mql_dtrtrs(char uplo, char trans, char diag, int n, int nrhs, double *a, int lda, double *b, int ldb)
{
    return LAPACKE_dtrtrs(LAPACK_COL_MAJOR, uplo, trans, diag, (lapack_int) n, (lapack_int) nrhs, a, (lapack_int) lda, b, (lapack_int) ldb);
}

_DLLAPI int __stdcall mql_dtrttf(char transr, char uplo, int n, double *a, int lda, double *arf)
{
    return LAPACKE_dtrttf(LAPACK_COL_MAJOR, transr, uplo, (lapack_int) n, a, (lapack_int) lda, arf);
}

_DLLAPI int __stdcall mql_dtrttp(char uplo, int n, double *a, int lda, double *ap)
{
    return LAPACKE_dtrttp(LAPACK_COL_MAJOR, uplo, (lapack_int) n, a, (lapack_int) lda, ap);
}

_DLLAPI int __stdcall mql_dtzrzf(int m, int n, double *a, int lda, double *tau)
{
    return LAPACKE_dtzrzf(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) n, a, (lapack_int) lda, tau);
}

#if defined TESTING
_DLLAPI int __stdcall mql_dlagsy(int n, int k, double *d, double *a, int lda, int *iseed)
{
    return LAPACKE_dlagsy(LAPACK_COL_MAJOR, (lapack_int) n, (lapack_int) k, d, a, (lapack_int) lda, iseed);
}
#endif

_DLLAPI double __stdcall mql_dlapy2(double x, double y)
{
    return LAPACKE_dlapy2(x, y);
}

_DLLAPI double __stdcall mql_dlapy3(double x, double y, double z)
{
    return LAPACKE_dlapy3(x, y, z);
}

_DLLAPI int __stdcall mql_dlartgp(double f, double g, double *cs, double *sn, double *r)
{
    return LAPACKE_dlartgp(f, g, cs, sn, r);
}

_DLLAPI int __stdcall mql_dlartgs(double x, double y, double sigma, double *cs, double *sn)
{
    return LAPACKE_dlartgs(x, y, sigma, cs, sn);
}

_DLLAPI int __stdcall mql_dbbcsd(char jobu1, char jobu2, char jobv1t, char jobv2t, char trans, int m, int p, int q, double *theta, double *phi, double *u1, int ldu1, double *u2, int ldu2, double *v1t, int ldv1t, double *v2t, int ldv2t, double *b11d, double *b11e, double *b12d, double *b12e, double *b21d, double *b21e, double *b22d, double *b22e)
{
    return LAPACKE_dbbcsd(LAPACK_COL_MAJOR, jobu1, jobu2, jobv1t, jobv2t, trans, (lapack_int) m, (lapack_int) p, (lapack_int) q, theta, phi, u1, (lapack_int) ldu1, u2, (lapack_int) ldu2, v1t, (lapack_int) ldv1t, v2t, (lapack_int) ldv2t, b11d, b11e, b12d, b12e, b21d, b21e, b22d, b22e);
}

_DLLAPI int __stdcall mql_dorbdb(char trans, char signs, int m, int p, int q, double *x11, int ldx11, double *x12, int ldx12, double *x21, int ldx21, double *x22, int ldx22, double *theta, double *phi, double *taup1, double *taup2, double *tauq1, double *tauq2)
{
    return LAPACKE_dorbdb(LAPACK_COL_MAJOR, trans, signs, (lapack_int) m, (lapack_int) p, (lapack_int) q, x11, (lapack_int) ldx11, x12, (lapack_int) ldx12, x21, (lapack_int) ldx21, x22, (lapack_int) ldx22, theta, phi, taup1, taup2, tauq1, tauq2);
}

_DLLAPI int __stdcall mql_dorcsd(char jobu1, char jobu2, char jobv1t, char jobv2t, char trans, char signs, int m, int p, int q, double *x11, int ldx11, double *x12, int ldx12, double *x21, int ldx21, double *x22, int ldx22, double *theta, double *u1, int ldu1, double *u2, int ldu2, double *v1t, int ldv1t, double *v2t, int ldv2t)
{
    return LAPACKE_dorcsd(LAPACK_COL_MAJOR, jobu1, jobu2, jobv1t, jobv2t, trans, signs, (lapack_int) m, (lapack_int) p, (lapack_int) q, x11, (lapack_int) ldx11, x12, (lapack_int) ldx12, x21, (lapack_int) ldx21, x22, (lapack_int) ldx22, theta, u1, (lapack_int) ldu1, u2, (lapack_int) ldu2, v1t, (lapack_int) ldv1t, v2t, (lapack_int) ldv2t);
}

_DLLAPI int __stdcall mql_dsyconv(char uplo, char way, int n, double *a, int lda, int *ipiv)
{
    return LAPACKE_dsyconv(LAPACK_COL_MAJOR, uplo, way, (lapack_int) n, a, (lapack_int) lda, ipiv);
}

_DLLAPI int __stdcall mql_dsyswapr(char uplo, int n, double *a, int i1, int i2)
{
    return LAPACKE_dsyswapr(LAPACK_COL_MAJOR, uplo, (lapack_int) n, a, (lapack_int) i1, (lapack_int) i2);
}

_DLLAPI int __stdcall mql_dsytri2(char uplo, int n, double *a, int lda, int *ipiv)
{
    return LAPACKE_dsytri2(LAPACK_COL_MAJOR, uplo, (lapack_int) n, a, (lapack_int) lda, ipiv);
}

_DLLAPI int __stdcall mql_dsytri2x(char uplo, int n, double *a, int lda, int *ipiv, int nb)
{
    return LAPACKE_dsytri2x(LAPACK_COL_MAJOR, uplo, (lapack_int) n, a, (lapack_int) lda, ipiv, (lapack_int) nb);
}

_DLLAPI int __stdcall mql_dsytrs2(char uplo, int n, int nrhs, double *a, int lda, int *ipiv, double *b, int ldb)
{
    return LAPACKE_dsytrs2(LAPACK_COL_MAJOR, uplo, (lapack_int) n, (lapack_int) nrhs, a, (lapack_int) lda, ipiv, b, (lapack_int) ldb);
}

_DLLAPI int __stdcall mql_dgemqrt(char side, char trans, int m, int n, int k, int nb, double *v, int ldv, double *t, int ldt, double *c, int ldc)
{
    return LAPACKE_dgemqrt(LAPACK_COL_MAJOR, side, trans, (lapack_int) m, (lapack_int) n, (lapack_int) k, (lapack_int) nb, v, (lapack_int) ldv, t, (lapack_int) ldt, c, (lapack_int) ldc);
}

_DLLAPI int __stdcall mql_dgeqrt(int m, int n, int nb, double *a, int lda, double *t, int ldt)
{
    return LAPACKE_dgeqrt(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) n, (lapack_int) nb, a, (lapack_int) lda, t, (lapack_int) ldt);
}

_DLLAPI int __stdcall mql_dgeqrt2(int m, int n, double *a, int lda, double *t, int ldt)
{
    return LAPACKE_dgeqrt2(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) n, a, (lapack_int) lda, t, (lapack_int) ldt);
}

_DLLAPI int __stdcall mql_dgeqrt3(int m, int n, double *a, int lda, double *t, int ldt)
{
    return LAPACKE_dgeqrt3(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) n, a, (lapack_int) lda, t, (lapack_int) ldt);
}

_DLLAPI int __stdcall mql_dtpmqrt(char side, char trans, int m, int n, int k, int l, int nb, double *v, int ldv, double *t, int ldt, double *a, int lda, double *b, int ldb)
{
    return LAPACKE_dtpmqrt(LAPACK_COL_MAJOR, side, trans, (lapack_int) m, (lapack_int) n, (lapack_int) k, (lapack_int) l, (lapack_int) nb, v, (lapack_int) ldv, t, (lapack_int) ldt, a, (lapack_int) lda, b, (lapack_int) ldb);
}

_DLLAPI int __stdcall mql_dtpqrt(int m, int n, int l, int nb, double *a, int lda, double *b, int ldb, double *t, int ldt)
{
    return LAPACKE_dtpqrt(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) n, (lapack_int) l, (lapack_int) nb, a, (lapack_int) lda, b, (lapack_int) ldb, t, (lapack_int) ldt);
}

_DLLAPI int __stdcall mql_dtpqrt2(int m, int n, int l, double *a, int lda, double *b, int ldb, double *t, int ldt)
{
    return LAPACKE_dtpqrt2(LAPACK_COL_MAJOR, (lapack_int) m, (lapack_int) n, (lapack_int) l, a, (lapack_int) lda, b, (lapack_int) ldb, t, (lapack_int) ldt);
}

_DLLAPI int __stdcall mql_dtprfb(char side, char trans, char direct, char storev, int m, int n, int k, int l, double *v, int ldv, double *t, int ldt, double *a, int lda, double *b, int ldb)
{
    return LAPACKE_dtprfb(LAPACK_COL_MAJOR, side, trans, direct, storev, (lapack_int) m, (lapack_int) n, (lapack_int) k, (lapack_int) l, v, (lapack_int) ldv, t, (lapack_int) ldt, a, (lapack_int) lda, b, (lapack_int) ldb);
}

/*
_DLLAPI int __stdcall mql_dsysv_rook(char uplo, int n, int nrhs, double *a, int lda, int *ipiv, double *b, int ldb)
{
    return LAPACKE_dsysv_rook(LAPACK_COL_MAJOR, uplo, (lapack_int) n, (lapack_int) nrhs, a, (lapack_int) lda, ipiv, b, (lapack_int) ldb);
}
*/
