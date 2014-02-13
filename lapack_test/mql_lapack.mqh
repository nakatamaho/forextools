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

#import "lapackinterface.dll"
double mql_ddot(int n, double &x[], int incx, double &y[], int incy);
double mql_dasum(int n, double &x[], int incx);
double mql_dnrm2(int N, double &X[], int incX);
double mql_dznrm2(int N, double &X[], int incX);
void mql_daxpy(int n, double alpha, double &x[], int incx, double &y[], int incy);
void mql_dcopy(int n, double &x[], int incx, double &y[], int incy);
void mql_dswap(int n, double &x[], int incx, double &y[], int incy);
void mql_drot(int N, double &X[], int incX, double &Y[], int incY, double c, double s);
void mql_drotg(double &a[], double &b[], double &c[], double &s[]);
void mql_drotm(int N, double &X[], int incX, double &Y[], int incY, double &P[]);
void mql_drotmg(double &d1, double &d2, double &b1, double b2, double &P);
void mql_dscal(int N, double alpha, double &X[], int incX);
void mql_dgemv(char trans, int m, int n, double alpha, double &a[], int lda, double &x[], int incx, double beta, double &y[], int incy);
void mql_dger(int M, int N, double alpha, double &X[], int incX, double &Y[], int incY, double &A[], int lda);
void mql_dtrsv(char Uplo, char TransA, char Diag, int N, double &A[], int lda, double &X[], int incX);
void mql_dtrmv(char Uplo, char TransA, char Diag, int N, double &A[], int lda, double &X[], int incX);
void mql_dsyr(char Uplo, int N, double alpha, double &X[], int incX, double &A[], int lda);
void mql_dgbmv(char TransA, int M, int N, int KL, int KU, double alpha, double &A[], int lda, double &X[], int incX, double beta, double &Y[], int incY);
void mql_dsbmv(char Uplo, int N, int K, double alpha, double &A, int lda, double &X[], int incX, double beta, double &Y[], int incY);
void mql_dtbmv(char Uplo, char TransA, char Diag, int N, int K, double &A[], int lda, double &X[], int incX);
void mql_dtbsv(char Uplo, char TransA, char Diag, int N, int K, double &A[], int lda, double &X[], int incX);
void mql_dtpmv(char Uplo, char TransA, char Diag, int N, double &Ap[], double &X[], int incX);
void mql_dtpsv(char Uplo, char TransA, char Diag, int N, double &Ap[], double &X[], int incX);
void mql_dsymv(char Uplo, int N, double alpha, double &A[], int lda, double &X[], int incX, double beta, double &Y[], int incY);
void mql_dspmv(char Uplo, int N, double alpha, double &Ap[], double &X[], int incX, double beta, double &Y[], int incY);
void mql_dspr(char Uplo, int N, double alpha, double &X, int incX, double &Ap[]); 
void mql_dspr2(char Uplo, int N, double alpha, double &X[], int incX, double &Y[], int incY, double &A[]);
void mql_dgemm(char TransA, char TransB, int M, int N, int K,  double alpha, double &A[], int lda, double &B[], int ldb, double beta, double &C[], int ldc);
void mql_dsymm(char Side, char Uplo, int M, int N, double alpha, double &A[], int lda, double &B[], int ldb, double beta, double &C[], int ldc);
void mql_dsyrk(char Uplo, char Trans, int N, int K, double alpha, double &A[], int lda, double beta, double &C[], int ldc);
void mql_dsyr2k(char Uplo, char Trans, int N, int K, double alpha, double &A[], int lda, double &B[], int ldb, double beta, double &C[], int ldc);
void mql_dtrmm(char Side, char Uplo, char TransA, char Diag, int M, int N, double alpha, double &A[], int lda, double &B[], int ldb);
void mql_dtrsm(char Side, char Uplo, char TransA, char Diag, int M, int N, double alpha, double &A[], int lda, double &B[], int ldb);

//LAPACK
int mql_dbdsdc(char uplo, char compq, int n, double &d[], double &e[], double &u[], int ldu, double &vt[], int ldvt, double &q[], int &iq[]);
int mql_dbdsqr(char uplo, int n, int ncvt, int nru, int ncc, double &d[], double &e[], double &vt[], int ldvt, double &u[], int ldu, double &c[], int ldc);
int mql_ddisna(char job, int m, int n, double &d[], double &sep[]);
int mql_dgbbrd(char vect, int m, int n, int ncc, int kl, int ku, double &ab[], int ldab, double &d[], double &e[], double &q[], int ldq, double &pt[], int ldpt, double &c[], int ldc);
int mql_dgbcon(char norm, int n, int kl, int ku, double &ab[], int ldab, int &ipiv[], double anorm, double &rcond);
int mql_dgbequ(int m, int n, int kl, int ku, double &ab[], int ldab, double &r[], double &c[], double &rowcnd, double &colcnd, double &amax);
int mql_dgbequb(int m, int n, int kl, int ku, double &ab[], int ldab, double &r[], double &c[], double &rowcnd, double &colcnd, double &amax);
int mql_dgbrfs(char trans, int n, int kl, int ku, int nrhs, double &ab[], int ldab, double &afb[], int ldafb, int &ipiv[], double &b[], int ldb, double &x[], int ldx, double &ferr[], double &berr[]);
int mql_dgbrfsx(char trans, char equed, int n, int kl, int ku, int nrhs, double &ab[], int ldab, double &afb[], int ldafb, int &ipiv[], double &r[], double &c[], double &b[], int ldb, double &x[], int ldx, double &rcond, double &berr[], int n_err_bnds, double &err_bnds_norm[], double &err_bnds_comp[], int nparams, double &params[]);
int mql_dgbsv(int n, int kl, int ku, int nrhs, double &ab[], int ldab, int &ipiv[], double &b[], int ldb);
int mql_dgbsvx(char fact, char trans, int n, int kl, int ku, int nrhs, double &ab[], int ldab, double &afb[], int ldafb, int &ipiv[], char equed, double &r[], double &c[], double &b[], int ldb, double &x[], int ldx, double &rcond, double &ferr[], double &berr[], double &rpivot[]);
int mql_dgbsvxx(char fact, char trans, int n, int kl, int ku, int nrhs, double &ab[], int ldab, double &afb[], int ldafb, int &ipiv[], char equed, double &r[], double &c[], double &b[], int ldb, double &x[], int ldx, double &rcond, double &rpvgrw, double &berr[], int n_err_bnds, double &err_bnds_norm[], double &err_bnds_comp[], int nparams, double &params[]);
int mql_dgbtrf(int m, int n, int kl, int ku, double &ab[], int ldab, int &ipiv[]);
int mql_dgbtrs(char trans, int n, int kl, int ku, int nrhs, double &ab[], int ldab, int &ipiv[], double &b[], int ldb);
int mql_dgebak(char job, char side, int n, int ilo, int ihi, double &scale[], int m, double &v[], int ldv);
int mql_dgebal(char job, int n, double &a[], int lda, int &ilo, int &ihi, double &scale[]);
int mql_dgebrd(int m, int n, double &a[], int lda, double &d[], double &e[], double &tauq[], double &taup[]);
int mql_dgecon(char norm, int n, double &a[], int lda, double anorm, double &rcond);
int mql_dgeequ(int m, int n, double &a[], int lda, double &r[], double &c[], double &rowcnd, double &colcnd, double &amax);
int mql_dgeequb(int m, int n, double &a[], int lda, double &r[], double &c[], double &rowcnd, double &colcnd, double &amax);
int mql_dgeev(char jobvl, char jobvr, int n, double &a[], int lda, double &wr[], double &wi[], double &vl[], int ldvl, double &vr[], int ldvr);
int mql_dgeevx(char balanc, char jobvl, char jobvr, char sense, int n, double &a[], int lda, double &wr[], double &wi[], double &vl[], int ldvl, double &vr[], int ldvr, int &ilo, int &ihi, double &scale[], double &abnrm, double &rconde, double &rcondv);
int mql_dgehrd(int n, int ilo, int ihi, double &a[], int lda, double &tau[]);
int mql_dgejsv(char joba, char jobu, char jobv, char jobr, char jobt, char jobp, int m, int n, double &a[], int lda, double &sva[], double &u[], int ldu, double &v[], int ldv, double &stat[], int &istat[]);
int mql_dgelq2(int m, int n, double &a[], int lda, double &tau[]);
int mql_dgelqf(int m, int n, double &a[], int lda, double &tau[]);
int mql_dgels(char trans, int m, int n, int nrhs, double &a[], int lda, double &b[], int ldb);
int mql_dgelsd(int m, int n, int nrhs, double &a[], int lda, double &b[], int ldb, double &s[], double rcond, int &rank);
int mql_dgelss(int m, int n, int nrhs, double &a[], int lda, double &b[], int ldb, double &s[], double rcond, int &rank);
int mql_dgelsy(int m, int n, int nrhs, double &a[], int lda, double &b[], int ldb, int &jpvt[], double rcond, int &rank);
int mql_dgeqlf(int m, int n, double &a[], int lda, double &tau[]);
int mql_dgeqp3(int m, int n, double &a[], int lda, int &jpvt[], double &tau[]);
int mql_dgeqpf(int m, int n, double &a[], int lda, int &jpvt[], double &tau[]);
int mql_dgeqr2(int m, int n, double &a[], int lda, double &tau[]);
int mql_dgeqrf(int m, int n, double &a[], int lda, double &tau[]);
int mql_dgeqrfp(int m, int n, double &a[], int lda, double &tau[]);
int mql_dgerfs(char trans, int n, int nrhs, double &a[], int lda, double &af[], int ldaf, int &ipiv[], double &b[], int ldb, double &x[], int ldx, double &ferr[], double &berr[]);
int mql_dgerfsx(char trans, char equed, int n, int nrhs, double &a[], int lda, double &af[], int ldaf, int &ipiv[], double &r[], double &c[], double &b[], int ldb, double &x[], int ldx, double &rcond, double &berr[], int n_err_bnds, double &err_bnds_norm[], double &err_bnds_comp[], int nparams, double &params[]);
int mql_dgerqf(int m, int n, double &a[], int lda, double &tau[]);
int mql_dgesdd(char jobz, int m, int n, double &a[], int lda, double &s[], double &u[], int ldu, double &vt[], int ldvt);
int mql_dgesv(int n, int nrhs, double &a[], int lda, int &ipiv[], double &b[], int ldb);
int mql_dsgesv(int n, int nrhs, double &a[], int lda, int &ipiv[], double &b[], int ldb, double &x[], int ldx, int &iter[]);
int mql_dgesvd(char jobu, char jobvt, int m, int n, double &a[], int lda, double &s[], double &u[], int ldu, double &vt[], int ldvt, double &superb[]);
int mql_dgesvj(char joba, char jobu, char jobv, int m, int n, double &a[], int lda, double &sva[], int mv, double &v[], int ldv, double &stat[]);
int mql_dgesvx(char fact, char trans, int n, int nrhs, double &a[], int lda, double &af[], int ldaf, int &ipiv[], char equed, double &r[], double &c[], double &b[], int ldb, double &x[], int ldx, double &rcond, double &ferr[], double &berr[], double &rpivot[]);
int mql_dgesvxx(char fact, char trans, int n, int nrhs, double &a[], int lda, double &af[], int ldaf, int &ipiv[], char equed, double &r[], double &c[], double &b[], int ldb, double &x[], int ldx, double &rcond, double &rpvgrw, double &berr[], int n_err_bnds, double &err_bnds_norm[], double &err_bnds_comp[], int nparams, double &params[]);
int mql_dgetf2(int m, int n, double &a[], int lda, int &ipiv[]);
int mql_dgetrf(int m, int n, double &a[], int lda, int &ipiv[]);
int mql_dgetri(int n, double &a[], int lda, int &ipiv[]);
int mql_dgetrs(char trans, int n, int nrhs, double &a[], int lda, int &ipiv[], double &b[], int ldb);
int mql_dggbak(char job, char side, int n, int ilo, int ihi, double &lscale, double &rscale, int m, double &v[], int ldv);
int mql_dggbal(char job, int n, double &a[], int lda, double &b[], int ldb, int &ilo, int &ihi, double &lscale, double &rscale);
int mql_dggev(char jobvl, char jobvr, int n, double &a[], int lda, double &b[], int ldb, double &alphar[], double &alphai[], double &beta[], double &vl[], int ldvl, double &vr[], int ldvr);
int mql_dggevx(char balanc, char jobvl, char jobvr, char sense, int n, double &a[], int lda, double &b[], int ldb, double &alphar[], double &alphai[], double &beta[], double &vl[], int ldvl, double &vr[], int ldvr, int &ilo, int &ihi, double &lscale, double &rscale, double &abnrm, double &bbnrm, double &rconde, double &rcondv);
int mql_dggglm(int n, int m, int p, double &a[], int lda, double &b[], int ldb, double &d[], double &x[], double &y[]);
int mql_dgghrd(char compq, char compz, int n, int ilo, int ihi, double &a[], int lda, double &b[], int ldb, double &q[], int ldq, double &z[], int ldz);
int mql_dgglse(int m, int n, int p, double &a[], int lda, double &b[], int ldb, double &c[], double &d[], double &x[]);
int mql_dggqrf(int n, int m, int p, double &a[], int lda, double &taua[], double &b[], int ldb, double &taub[]);
int mql_dggrqf(int m, int p, int n, double &a[], int lda, double &taua[], double &b[], int ldb, double &taub[]);
int mql_dggsvd(char jobu, char jobv, char jobq, int m, int n, int p, int &k, int &l, double &a[], int lda, double &b[], int ldb, double &alpha[], double &beta[], double &u[], int ldu, double &v[], int ldv, double &q[], int ldq, int &iwork[]);
int mql_dggsvp(char jobu, char jobv, char jobq, int m, int p, int n, double &a[], int lda, double &b[], int ldb, double tola, double tolb, int &k, int &l, double &u[], int ldu, double &v[], int ldv, double &q[], int ldq);
int mql_dgtcon(char norm, int n, double &dl[], double &d[], double &du[], double &du2[], int &ipiv[], double anorm, double &rcond);
int mql_dgtrfs(char trans, int n, int nrhs, double &dl[], double &d[], double &du[], double &dlf[], double &df[], double &duf[], double &du2[], int &ipiv[], double &b[], int ldb, double &x[], int ldx, double &ferr[], double &berr[]);
int mql_dgtsv(int n, int nrhs, double &dl[], double &d[], double &du[], double &b[], int ldb);
int mql_dgtsvx(char fact, char trans, int n, int nrhs, double &dl[], double &d[], double &du[], double &dlf[], double &df[], double &duf[], double &du2[], int &ipiv[], double &b[], int ldb, double &x[], int ldx, double &rcond, double &ferr[], double &berr[]);
int mql_dgttrf(int n, double &dl[], double &d[], double &du[], double &du2[], int &ipiv[]);
int mql_dgttrs(char trans, int n, int nrhs, double &dl[], double &d[], double &du[], double &du2[], int &ipiv[], double &b[], int ldb);
int mql_dhgeqz(char job, char compq, char compz, int n, int ilo, int ihi, double &h[], int ldh, double &t[], int ldt, double &alphar[], double &alphai[], double &beta[], double &q[], int ldq, double &z[], int ldz);
int mql_dhseqr(char job, char compz, int n, int ilo, int ihi, double &h[], int ldh, double &wr[], double &wi[], double &z[], int ldz);
int mql_dlacn2(int n, double &v[], double &x[], int &isgn[], double &est, int &kase, int &isave[]);
int mql_dlacpy(char uplo, int m, int n, double &a[], int lda, double &b[], int ldb);
int mql_dlagge(int m, int n, int kl, int ku, double &d[], double &a[], int lda, int &iseed[]);
double mql_dlamch(char cmach);
double mql_dlange(char norm, int m, int n, double &a[], int lda);
double mql_dlansy(char norm, char uplo, int n, double &a[], int lda);
double mql_dlantr(char norm, char uplo, char diag, int m, int n, double &a[], int lda);
int mql_dlarfb(char side, char trans, char direct, char storev, int m, int n, int k, double &v[], int ldv, double &t[], int ldt, double &c[], int ldc);
int mql_dlarfg(int n, double &alpha[], double &x[], int incx, double &tau[]);
int mql_dlarft(char direct, char storev, int n, int k, double &v[], int ldv, double &tau[], double &t[], int ldt);
int mql_dlarfx(char side, int m, int n, double &v[], double tau, double &c[], int ldc, double &work[]);
int mql_dlarnv(int idist, int &iseed[], int n, double &x[]);
int mql_dlaset(char uplo, int m, int n, double alpha, double beta, double &a[], int lda);
int mql_dlasrt(char id, int n, double &d[]);
int mql_dlaswp(int n, double &a[], int lda, int k1, int k2, int &ipiv[], int incx);
int mql_dlatms(int m, int n, char dist, int &iseed[], char sym, double &d[], int mode, double cond, double dmax, int kl, int ku, char pack, double &a[], int lda);
int mql_dlauum(char uplo, int n, double &a[], int lda);
int mql_dopgtr(char uplo, int n, double &ap[], double &tau[], double &q[], int ldq);
int mql_dopmtr(char side, char uplo, char trans, int m, int n, double &ap[], double &tau[], double &c[], int ldc);
int mql_dorgbr(char vect, int m, int n, int k, double &a[], int lda, double &tau[]);
int mql_dorghr(int n, int ilo, int ihi, double &a[], int lda, double &tau[]);
int mql_dorglq(int m, int n, int k, double &a[], int lda, double &tau[]);
int mql_dorgql(int m, int n, int k, double &a[], int lda, double &tau[]);
int mql_dorgqr(int m, int n, int k, double &a[], int lda, double &tau[]);
int mql_dorgrq(int m, int n, int k, double &a[], int lda, double &tau[]);
int mql_dorgtr(char uplo, int n, double &a[], int lda, double &tau[]);
int mql_dormbr(char vect, char side, char trans, int m, int n, int k, double &a[], int lda, double &tau[], double &c[], int ldc);
int mql_dormhr(char side, char trans, int m, int n, int ilo, int ihi, double &a[], int lda, double &tau[], double &c[], int ldc);
int mql_dormlq(char side, char trans, int m, int n, int k, double &a[], int lda, double &tau[], double &c[], int ldc);
int mql_dormql(char side, char trans, int m, int n, int k, double &a[], int lda, double &tau[], double &c[], int ldc);
int mql_dormqr(char side, char trans, int m, int n, int k, double &a[], int lda, double &tau[], double &c[], int ldc);
int mql_dormrq(char side, char trans, int m, int n, int k, double &a[], int lda, double &tau[], double &c[], int ldc);
int mql_dormrz(char side, char trans, int m, int n, int k, int l, double &a[], int lda, double &tau[], double &c[], int ldc);
int mql_dormtr(char side, char uplo, char trans, int m, int n, double &a[], int lda, double &tau[], double &c[], int ldc);
int mql_dpbcon(char uplo, int n, int kd, double &ab[], int ldab, double anorm, double &rcond);
int mql_dpbequ(char uplo, int n, int kd, double &ab[], int ldab, double &s[], double &scond, double &amax);
int mql_dpbrfs(char uplo, int n, int kd, int nrhs, double &ab[], int ldab, double &afb[], int ldafb, double &b[], int ldb, double &x[], int ldx, double &ferr[], double &berr[]);
int mql_dpbstf(char uplo, int n, int kb, double &bb[], int ldbb);
int mql_dpbsv(char uplo, int n, int kd, int nrhs, double &ab[], int ldab, double &b[], int ldb);
int mql_dpbsvx(char fact, char uplo, int n, int kd, int nrhs, double &ab[], int ldab, double &afb[], int ldafb, char equed, double &s[], double &b[], int ldb, double &x[], int ldx, double &rcond, double &ferr[], double &berr[]);
int mql_dpbtrf(char uplo, int n, int kd, double &ab[], int ldab);
int mql_dpbtrs(char uplo, int n, int kd, int nrhs, double &ab[], int ldab, double &b[], int ldb);
int mql_dpftrf(char transr, char uplo, int n, double &a[]);
int mql_dpftri(char transr, char uplo, int n, double &a[]);
int mql_dpftrs(char transr, char uplo, int n, int nrhs, double &a[], double &b[], int ldb);
int mql_dpocon(char uplo, int n, double &a[], int lda, double anorm, double &rcond);
int mql_dpoequ(int n, double &a[], int lda, double &s[], double &scond, double &amax);
int mql_dpoequb(int n, double &a[], int lda, double &s[], double &scond, double &amax);
int mql_dporfs(char uplo, int n, int nrhs, double &a[], int lda, double &af[], int ldaf, double &b[], int ldb, double &x[], int ldx, double &ferr[], double &berr[]);
int mql_dporfsx(char uplo, char equed, int n, int nrhs, double &a[], int lda, double &af[], int ldaf, double &s[], double &b[], int ldb, double &x[], int ldx, double &rcond, double &berr[], int n_err_bnds, double &err_bnds_norm[], double &err_bnds_comp[], int nparams, double &params[]);
int mql_dposv(char uplo, int n, int nrhs, double &a[], int lda, double &b[], int ldb);
int mql_dsposv(char uplo, int n, int nrhs, double &a[], int lda, double &b[], int ldb, double &x[], int ldx, int &iter[]);
int mql_dposvx(char fact, char uplo, int n, int nrhs, double &a[], int lda, double &af[], int ldaf, char equed, double &s[], double &b[], int ldb, double &x[], int ldx, double &rcond, double &ferr[], double &berr[]);
int mql_dposvxx(char fact, char uplo, int n, int nrhs, double &a[], int lda, double &af[], int ldaf, char equed, double &s[], double &b[], int ldb, double &x[], int ldx, double &rcond, double &rpvgrw, double &berr[], int n_err_bnds, double &err_bnds_norm[], double &err_bnds_comp[], int nparams, double &params[]);
int mql_dpotrf(char uplo, int n, double &a[], int lda);
int mql_dpotri(char uplo, int n, double &a[], int lda);
int mql_dpotrs(char uplo, int n, int nrhs, double &a[], int lda, double &b[], int ldb);
int mql_dppcon(char uplo, int n, double &ap[], double anorm, double &rcond);
int mql_dppequ(char uplo, int n, double &ap[], double &s[], double &scond, double &amax);
int mql_dpprfs(char uplo, int n, int nrhs, double &ap[], double &afp[], double &b[], int ldb, double &x[], int ldx, double &ferr[], double &berr[]);
int mql_dppsv(char uplo, int n, int nrhs, double &ap[], double &b[], int ldb);
int mql_dppsvx(char fact, char uplo, int n, int nrhs, double &ap[], double &afp[], char equed, double &s[], double &b[], int ldb, double &x[], int ldx, double &rcond, double &ferr[], double &berr[]);
int mql_dpptrf(char uplo, int n, double &ap[]);
int mql_dpptri(char uplo, int n, double &ap[]);
int mql_dpptrs(char uplo, int n, int nrhs, double &ap[], double &b[], int ldb);
int mql_dpstrf(char uplo, int n, double &a[], int lda, int &piv[], int &rank, double tol);
int mql_dptcon(int n, double &d[], double &e[], double anorm, double &rcond);
int mql_dpteqr(char compz, int n, double &d[], double &e[], double &z[], int ldz);
int mql_dptrfs(int n, int nrhs, double &d[], double &e[], double &df[], double &ef[], double &b[], int ldb, double &x[], int ldx, double &ferr[], double &berr[]);
int mql_dptsv(int n, int nrhs, double &d[], double &e[], double &b[], int ldb);
int mql_dptsvx(char fact, int n, int nrhs, double &d[], double &e[], double &df[], double &ef[], double &b[], int ldb, double &x[], int ldx, double &rcond, double &ferr[], double &berr[]);
int mql_dpttrf(int n, double &d[], double &e[]);
int mql_dpttrs(int n, int nrhs, double &d[], double &e[], double &b[], int ldb);
int mql_dsbev(char jobz, char uplo, int n, int kd, double &ab[], int ldab, double &w[], double &z[], int ldz);
int mql_dsbevd(char jobz, char uplo, int n, int kd, double &ab[], int ldab, double &w[], double &z[], int ldz);
int mql_dsbevx(char jobz, char range, char uplo, int n, int kd, double &ab[], int ldab, double &q[], int ldq, double vl, double vu, int il, int iu, double abstol, int &m, double &w[], double &z[], int ldz, int &ifail[]);
int mql_dsbgst(char vect, char uplo, int n, int ka, int kb, double &ab[], int ldab, double &bb[], int ldbb, double &x[], int ldx);
int mql_dsbgv(char jobz, char uplo, int n, int ka, int kb, double &ab[], int ldab, double &bb[], int ldbb, double &w[], double &z[], int ldz);
int mql_dsbgvd(char jobz, char uplo, int n, int ka, int kb, double &ab[], int ldab, double &bb[], int ldbb, double &w[], double &z[], int ldz);
int mql_dsbgvx(char jobz, char range, char uplo, int n, int ka, int kb, double &ab[], int ldab, double &bb[], int ldbb, double &q[], int ldq, double vl, double vu, int il, int iu, double abstol, int &m, double &w[], double &z[], int ldz, int &ifail[]);
int mql_dsbtrd(char vect, char uplo, int n, int kd, double &ab[], int ldab, double &d[], double &e[], double &q[], int ldq);
int mql_dsfrk(char transr, char uplo, char trans, int n, int k, double alpha, double &a[], int lda, double beta, double &c[]);
int mql_dspcon(char uplo, int n, double &ap[], int &ipiv[], double anorm, double &rcond);
int mql_dspev(char jobz, char uplo, int n, double &ap[], double &w[], double &z[], int ldz);
int mql_dspevd(char jobz, char uplo, int n, double &ap[], double &w[], double &z[], int ldz);
int mql_dspevx(char jobz, char range, char uplo, int n, double &ap[], double vl, double vu, int il, int iu, double abstol, int &m, double &w[], double &z[], int ldz, int &ifail[]);
int mql_dspgst(int itype, char uplo, int n, double &ap[], double &bp[]);
int mql_dspgv(int itype, char jobz, char uplo, int n, double &ap[], double &bp[], double &w[], double &z[], int ldz);
int mql_dspgvd(int itype, char jobz, char uplo, int n, double &ap[], double &bp[], double &w[], double &z[], int ldz);
int mql_dspgvx(int itype, char jobz, char range, char uplo, int n, double &ap[], double &bp[], double vl, double vu, int il, int iu, double abstol, int &m, double &w[], double &z[], int ldz, int &ifail[]);
int mql_dsprfs(char uplo, int n, int nrhs, double &ap[], double &afp[], int &ipiv[], double &b[], int ldb, double &x[], int ldx, double &ferr[], double &berr[]);
int mql_dspsv(char uplo, int n, int nrhs, double &ap[], int &ipiv[], double &b[], int ldb);
int mql_dspsvx(char fact, char uplo, int n, int nrhs, double &ap[], double &afp[], int &ipiv[], double &b[], int ldb, double &x[], int ldx, double &rcond, double &ferr[], double &berr[]);
int mql_dsptrd(char uplo, int n, double &ap[], double &d[], double &e[], double &tau[]);
int mql_dsptrf(char uplo, int n, double &ap[], int &ipiv[]);
int mql_dsptri(char uplo, int n, double &ap[], int &ipiv[]);
int mql_dsptrs(char uplo, int n, int nrhs, double &ap[], int &ipiv[], double &b[], int ldb);
int mql_dstebz(char range, char order, int n, double vl, double vu, int il, int iu, double abstol, double &d[], double &e[], int &m, int &nsplit, double &w[], int &iblock[], int &isplit[]);
int mql_dstedc(char compz, int n, double &d[], double &e[], double &z[], int ldz);
int mql_dstegr(char jobz, char range, int n, double &d[], double &e[], double vl, double vu, int il, int iu, double abstol, int &m, double &w[], double &z[], int ldz, int &isuppz);
int mql_dstein(int n, double &d[], double &e[], int m, double &w[], int &iblock[], int &isplit[], double &z[], int ldz, int &ifailv[]);
int mql_dsteqr(char compz, int n, double &d[], double &e[], double &z[], int ldz);
int mql_dsterf(int n, double &d[], double &e[]);
int mql_dstev(char jobz, int n, double &d[], double &e[], double &z[], int ldz);
int mql_dstevd(char jobz, int n, double &d[], double &e[], double &z[], int ldz);
int mql_dstevr(char jobz, char range, int n, double &d[], double &e[], double vl, double vu, int il, int iu, double abstol, int &m, double &w[], double &z[], int ldz, int &isuppz);
int mql_dstevx(char jobz, char range, int n, double &d[], double &e[], double vl, double vu, int il, int iu, double abstol, int &m, double &w[], double &z[], int ldz, int &ifail[]);
int mql_dsycon(char uplo, int n, double &a[], int lda, int &ipiv[], double anorm, double &rcond);
int mql_dsyequb(char uplo, int n, double &a[], int lda, double &s[], double &scond, double &amax);
int mql_dsyev(char jobz, char uplo, int n, double &a[], int lda, double &w[]);
int mql_dsyevd(char jobz, char uplo, int n, double &a[], int lda, double &w[]);
int mql_dsyevr(char jobz, char range, char uplo, int n, double &a[], int lda, double vl, double vu, int il, int iu, double abstol, int &m, double &w[], double &z[], int ldz, int &isuppz);
int mql_dsyevx(char jobz, char range, char uplo, int n, double &a[], int lda, double vl, double vu, int il, int iu, double abstol, int &m, double &w[], double &z[], int ldz, int &ifail[]);
int mql_dsygst(int itype, char uplo, int n, double &a[], int lda, double &b[], int ldb);
int mql_dsygv(int itype, char jobz, char uplo, int n, double &a[], int lda, double &b[], int ldb, double &w[]);
int mql_dsygvd(int itype, char jobz, char uplo, int n, double &a[], int lda, double &b[], int ldb, double &w[]);
int mql_dsygvx(int itype, char jobz, char range, char uplo, int n, double &a[], int lda, double &b[], int ldb, double vl, double vu, int il, int iu, double abstol, int &m, double &w[], double &z[], int ldz, int &ifail[]);
int mql_dsyrfs(char uplo, int n, int nrhs, double &a[], int lda, double &af[], int ldaf, int &ipiv[], double &b[], int ldb, double &x[], int ldx, double &ferr[], double &berr[]);
int mql_dsyrfsx(char uplo, char equed, int n, int nrhs, double &a[], int lda, double &af[], int ldaf, int &ipiv[], double &s[], double &b[], int ldb, double &x[], int ldx, double &rcond, double &berr[], int n_err_bnds, double &err_bnds_norm[], double &err_bnds_comp[], int nparams, double &params[]);
int mql_dsysv(char uplo, int n, int nrhs, double &a[], int lda, int &ipiv[], double &b[], int ldb);
int mql_dsysvx(char fact, char uplo, int n, int nrhs, double &a[], int lda, double &af[], int ldaf, int &ipiv[], double &b[], int ldb, double &x[], int ldx, double &rcond, double &ferr[], double &berr[]);
int mql_dsysvxx(char fact, char uplo, int n, int nrhs, double &a[], int lda, double &af[], int ldaf, int &ipiv[], char equed, double &s[], double &b[], int ldb, double &x[], int ldx, double &rcond, double &rpvgrw, double &berr[], int n_err_bnds, double &err_bnds_norm[], double &err_bnds_comp[], int nparams, double &params[]);
int mql_dsytrd(char uplo, int n, double &a[], int lda, double &d[], double &e[], double &tau[]);
int mql_dsytrf(char uplo, int n, double &a[], int lda, int &ipiv[]);
int mql_dsytri(char uplo, int n, double &a[], int lda, int &ipiv[]);
int mql_dsytrs(char uplo, int n, int nrhs, double &a[], int lda, int &ipiv[], double &b[], int ldb);
int mql_dtbcon(char norm, char uplo, char diag, int n, int kd, double &ab[], int ldab, double &rcond);
int mql_dtbrfs(char uplo, char trans, char diag, int n, int kd, int nrhs, double &ab[], int ldab, double &b[], int ldb, double &x[], int ldx, double &ferr[], double &berr[]);
int mql_dtbtrs(char uplo, char trans, char diag, int n, int kd, int nrhs, double &ab[], int ldab, double &b[], int ldb);
int mql_dtfsm(char transr, char side, char uplo, char trans, char diag, int m, int n, double alpha, double &a[], double &b[], int ldb);
int mql_dtftri(char transr, char uplo, char diag, int n, double &a[]);
int mql_dtfttp(char transr, char uplo, int n, double &arf[], double &ap[]);
int mql_dtfttr(char transr, char uplo, int n, double &arf[], double &a[], int lda);
int mql_dtgsja(char jobu, char jobv, char jobq, int m, int p, int n, int k, int l, double &a[], int lda, double &b[], int ldb, double tola, double tolb, double &alpha[], double &beta[], double &u[], int ldu, double &v[], int ldv, double &q[], int ldq, int &ncycle[]);
int mql_dtgsyl(char trans, int ijob, int m, int n, double &a[], int lda, double &b[], int ldb, double &c[], int ldc, double &d[], int ldd, double &e[], int lde, double &f[], int ldf, double &scale[], double &dif);
int mql_dtpcon(char norm, char uplo, char diag, int n, double &ap[], double &rcond);
int mql_dtprfs(char uplo, char trans, char diag, int n, int nrhs, double &ap[], double &b[], int ldb, double &x[], int ldx, double &ferr[], double &berr[]);
int mql_dtptri(char uplo, char diag, int n, double &ap[]);
int mql_dtptrs(char uplo, char trans, char diag, int n, int nrhs, double &ap[], double &b[], int ldb);
int mql_dtpttf(char transr, char uplo, int n, double &ap[], double &arf[]);
int mql_dtpttr(char uplo, int n, double &ap[], double &a[], int lda);
int mql_dtrcon(char norm, char uplo, char diag, int n, double &a[], int lda, double &rcond);
int mql_dtrexc(char compq, int n, double &t[], int ldt, double &q[], int ldq, int &ifst, int &ilst);
int mql_dtrrfs(char uplo, char trans, char diag, int n, int nrhs, double &a[], int lda, double &b[], int ldb, double &x[], int ldx, double &ferr[], double &berr[]);
int mql_dtrsyl(char trana, char tranb, int isgn, int m, int n, double &a[], int lda, double &b[], int ldb, double &c[], int ldc, double &scale[]);
int mql_dtrtri(char uplo, char diag, int n, double &a[], int lda);
int mql_dtrtrs(char uplo, char trans, char diag, int n, int nrhs, double &a[], int lda, double &b[], int ldb);
int mql_dtrttf(char transr, char uplo, int n, double &a[], int lda, double &arf[]);
int mql_dtrttp(char uplo, int n, double &a[], int lda, double &ap[]);
int mql_dtzrzf(int m, int n, double &a[], int lda, double &tau[]);
int mql_dlagsy(int n, int k, double &d[], double &a[], int lda, int &iseed[]);
double mql_dlapy2(double x, double y);
double mql_dlapy3(double x, double y, double z);
int mql_dlartgp(double f, double g, double &cs, double &sn, double &r[]);
int mql_dlartgs(double x, double y, double sigma, double &cs, double &sn);
int mql_dbbcsd(char jobu1, char jobu2, char jobv1t, char jobv2t, char trans, int m, int p, int q, double &theta[], double &phi[], double &u1[], int ldu1, double &u2[], int ldu2, double &v1t[], int ldv1t, double &v2t[], int ldv2t, double &b11d[], double &b11e[], double &b12d[], double &b12e[], double &b21d[], double &b21e[], double &b22d[], double &b22e[]);
int mql_dorbdb(char trans, char signs, int m, int p, int q, double &x11[], int ldx11, double &x12[], int ldx12, double &x21[], int ldx21, double &x22[], int ldx22, double &theta[], double &phi[], double &taup1[], double &taup2[], double &tauq1[], double &tauq2[]);
int mql_dorcsd(char jobu1, char jobu2, char jobv1t, char jobv2t, char trans, char signs, int m, int p, int q, double &x11[], int ldx11, double &x12[], int ldx12, double &x21[], int ldx21, double &x22[], int ldx22, double &theta[], double &u1[], int ldu1, double &u2[], int ldu2, double &v1t[], int ldv1t, double &v2t[], int ldv2t);
int mql_dsyconv(char uplo, char way, int n, double &a[], int lda, int &ipiv[]);
int mql_dsyswapr(char uplo, int n, double &a[], int i1, int i2);
int mql_dsytri2(char uplo, int n, double &a[], int lda, int &ipiv[]);
int mql_dsytri2x(char uplo, int n, double &a[], int lda, int &ipiv[], int nb);
int mql_dsytrs2(char uplo, int n, int nrhs, double &a[], int lda, int &ipiv[], double &b[], int ldb);
int mql_dgemqrt(char side, char trans, int m, int n, int k, int nb, double &v[], int ldv, double &t[], int ldt, double &c[], int ldc);
int mql_dgeqrt(int m, int n, int nb, double &a[], int lda, double &t[], int ldt);
int mql_dgeqrt2(int m, int n, double &a[], int lda, double &t[], int ldt);
int mql_dgeqrt3(int m, int n, double &a[], int lda, double &t[], int ldt);
int mql_dtpmqrt(char side, char trans, int m, int n, int k, int l, int nb, double &v[], int ldv, double &t[], int ldt, double &a[], int lda, double &b[], int ldb);
int mql_dtpqrt(int m, int n, int l, int nb, double &a[], int lda, double &b[], int ldb, double &t[], int ldt);
int mql_dtpqrt2(int m, int n, int l, double &a[], int lda, double &b[], int ldb, double &t[], int ldt);
int mql_dtprfb(char side, char trans, char direct, char storev, int m, int n, int k, int l, double &v[], int ldv, double &t[], int ldt, double &a[], int lda, double &b[], int ldb);
int mql_dsysv_rook(char uplo, int n, int nrhs, double &a[], int lda, int &ipiv[], double &b[], int ldb);

#import

