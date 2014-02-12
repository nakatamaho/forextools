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

#ifndef _BLAS_DOUBLEONLY_H_
#define _BLAS_DOUBLEONLY_H_

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define F77_FUNC(name,NAME) name ## _

/* As F77_FUNC, but for C identifiers containing underscores. */
#define F77_FUNC_(name,NAME) name ## _

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

#define F77_RET_I int
#define F77_RET_L int
#define F77_RET_D double

#ifdef __cplusplus
extern "C" {
#endif

//BLAS part
#define dasum_f77 F77_FUNC (dasum, DASUM)
#define daxpy_f77 F77_FUNC (daxpy, DAXPY)
#define dcabs1_f77 F77_FUNC (dcabs1, DCABS1)
#define dcopy_f77 F77_FUNC (dcopy, DCOPY)
#define ddot_f77 F77_FUNC (ddot, DDOT)
#define dgbmv_f77 F77_FUNC (dgbmv, DGBMV)
#define dgemm_f77 F77_FUNC (dgemm, DGEMM)
#define dgemv_f77 F77_FUNC (dgemv, DGEMV)
#define dger_f77 F77_FUNC (dger, DGER)
#define dnrm2_f77 F77_FUNC (dnrm2, DNRM2)
#define drot_f77 F77_FUNC (drot, DROT)
#define drotg_f77 F77_FUNC (drotg, DROTG)
#define drotm_f77 F77_FUNC (drotm, DROTM)
#define drotmg_f77 F77_FUNC (drotmg, DROTMG)
#define dsbmv_f77 F77_FUNC (dsbmv, DSBMV)
#define dscal_f77 F77_FUNC (dscal, DSCAL)
#define dsdot_f77 F77_FUNC (dsdot, DSDOT)
#define dspmv_f77 F77_FUNC (dspmv, DSPMV)
#define dspr_f77 F77_FUNC (dspr, DSPR)
#define dspr2_f77 F77_FUNC (dspr2, DSPR2)
#define dswap_f77 F77_FUNC (dswap, DSWAP)
#define dsymm_f77 F77_FUNC (dsymm, DSYMM)
#define dsymv_f77 F77_FUNC (dsymv, DSYMV)
#define dsyr_f77 F77_FUNC (dsyr, DSYR)
#define dsyr2_f77 F77_FUNC (dsyr2, DSYR2)
#define dsyr2k_f77 F77_FUNC (dsyr2k, DSYR2K)
#define dsyrk_f77 F77_FUNC (dsyrk, DSYRK)
#define dtbmv_f77 F77_FUNC (dtbmv, DTBMV)
#define dtbsv_f77 F77_FUNC (dtbsv, DTBSV)
#define dtpmv_f77 F77_FUNC (dtpmv, DTPMV)
#define dtpsv_f77 F77_FUNC (dtpsv, DTPSV)
#define dtrmm_f77 F77_FUNC (dtrmm, DTRMM)
#define dtrmv_f77 F77_FUNC (dtrmv, DTRMV)
#define dtrsm_f77 F77_FUNC (dtrsm, DTRSM)
#define dtrsv_f77 F77_FUNC (dtrsv, DTRSV)
#define dzasum_f77 F77_FUNC (dzasum, DZASUM)
#define dznrm2_f77 F77_FUNC (dznrm2, DZNRM2)
#define idamax_f77 F77_FUNC (idamax, IDAMAX)
#define xerbla_f77 F77_FUNC (xerbla, XERBLA)

F77_RET_D dasum_f77(int *n, double *dx, int *incx);
F77_RET_I daxpy_f77(int *n, double *da, double *dx, int *incx, double *dy, int *incy);
F77_RET_I dcopy_f77(int *n, double *dx, int *incx, double *dy, int *incy);
F77_RET_D ddot_f77(int *n, double *dx, int *incx, double *dy, int *incy);
F77_RET_I dgbmv_f77(const char *trans, int *m, int *n, int *kl, int *ku, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);
F77_RET_I dgemm_f77(const char *transa, const char *transb, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);
F77_RET_I dgemv_f77(const char *trans, int *m, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);
F77_RET_I dger_f77(int *m, int *n, double *alpha, double *x, int *incx, double *y, int *incy, double *a, int *lda);
F77_RET_D dnrm2_f77(int *n, double *x, int *incx);
F77_RET_I drot_f77(int *n, double *dx, int *incx, double *dy, int *incy, double *c, double *s);
F77_RET_I drotg_f77(double *da, double *db, double *c, double *s);
F77_RET_I drotm_f77(int *n, double *dx, int *incx, double *dy, int *incy, double *dparam);
F77_RET_I drotmg_f77(double *dd1, double *dd2, double *dx1, double *dy1, double *dparam);
F77_RET_I dsbmv_f77(const char *uplo, int *n, int *k, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);
F77_RET_I dscal_f77(int *n, double *da, double *dx, int *incx);
F77_RET_I dspmv_f77(const char *uplo, int *n, double *alpha, double *ap, double *x, int *incx, double *beta, double *y, int *incy);
F77_RET_I dspr_f77(const char *uplo, int *n, double *alpha, double *x, int *incx, double *ap);
F77_RET_I dspr2_f77(const char *uplo, int *n, double *alpha, double *x, int *incx, double *y, int *incy, double *ap);
F77_RET_I dswap_f77(int *n, double *dx, int *incx, double *dy, int *incy);
F77_RET_I dsymm_f77(const char *side, const char *uplo, int *m, int *n, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);
F77_RET_I dsymv_f77(const char *uplo, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);
F77_RET_I dsyr_f77(const char *uplo, int *n, double *alpha, double *x, int *incx, double *a, int *lda);
F77_RET_I dsyr2_f77(const char *uplo, int *n, double *alpha, double *x, int *incx, double *y, int *incy, double *a, int *lda);
F77_RET_I dsyr2k_f77(const char *uplo, const char *trans, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);
F77_RET_I dsyrk_f77(const char *uplo, const char *trans, int *n, int *k, double *alpha, double *a, int *lda, double *beta, double *c, int *ldc);
F77_RET_I dtbmv_f77(const char *uplo, const char *trans, const char *diag, int *n, int *k, double *a, int *lda, double *x, int *incx);
F77_RET_I dtbsv_f77(const char *uplo, const char *trans, const char *diag, int *n, int *k, double *a, int *lda, double *x, int *incx);
F77_RET_I dtpmv_f77(const char *uplo, const char *trans, const char *diag, int *n, double *ap, double *x, int *incx);
F77_RET_I dtpsv_f77(const char *uplo, const char *trans, const char *diag, int *n, double *ap, double *x, int *incx);
F77_RET_I dtrmm_f77(const char *side, const char *uplo, const char *transa, const char *diag, int *m, int *n, double *alpha, double *a, int *lda, double *b, int *ldb);
F77_RET_I dtrmv_f77(const char *uplo, const char *trans, const char *diag, int *n, double *a, int *lda, double *x, int *incx);
F77_RET_I dtrsm_f77(const char *side, const char *uplo, const char *transa, const char *diag, int *m, int *n, double *alpha, double *a, int *lda, double *b, int *ldb);
F77_RET_I dtrsv_f77(const char *uplo, const char *trans, const char *diag, int *n, double *a, int *lda, double *x, int *incx);
F77_RET_I idamax_f77(int *n, double *dx, int *incx);

#ifdef __cplusplus
}
#endif
#endif
