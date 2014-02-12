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

int mql_dsyev(char jobz, char uplo, int n, double &a[], int lda, double &w[]);
int mql_dgesvd(char jobu, char jobvt, int m, int n, double &a[], int lda, double &s[], double &u[], int ldu, double &vt[], int ldvt, double &superb[]);
int mql_dgesdd(char jobz, int m, int n, double &a[], int lda, double &s[], double &u[], int ldu, double &vt[], int ldvt);
int mql_dgels(char trans, int m, int n, int nrhs, double &a[], int lda, double &b[], int ldb);
#import

