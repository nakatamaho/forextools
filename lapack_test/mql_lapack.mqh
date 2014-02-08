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

enum CBLAS_ORDER     {CblasRowMajor=101, CblasColMajor=102};
enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113, CblasConjNoTrans=114};
enum CBLAS_UPLO      {CblasUpper=121, CblasLower=122};
enum CBLAS_DIAG      {CblasNonUnit=131, CblasUnit=132};
enum CBLAS_SIDE      {CblasLeft=141, CblasRight=142};

#define LAPACK_ROW_MAJOR               101
#define LAPACK_COL_MAJOR               102
#define LAPACK_WORK_MEMORY_ERROR       -1010
#define LAPACK_TRANSPOSE_MEMORY_ERROR  -1011

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
void mql_dgemv(CBLAS_ORDER order, CBLAS_TRANSPOSE trans, int m, int n, double alpha, double &a[], int lda, double &x[], int incx, double beta, double &y[], int incy);
void mql_dger(CBLAS_ORDER order, int M, int N, double alpha, double &X[], int incX, double &Y[], int incY, double &A[], int lda);
void mql_dtrsv(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, int N, double &A[], int lda, double &X[], int incX);
void mql_dtrmv(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, int N, double &A[], int lda, double &X[], int incX);
void mql_dsyr(CBLAS_ORDER order, CBLAS_UPLO Uplo, int N, double alpha, double &X[], int incX, double &A[], int lda);
void mql_dgbmv(CBLAS_ORDER order, CBLAS_TRANSPOSE TransA, int M, int N, int KL, int KU, double alpha, double &A[], int lda, double &X[], int incX, double beta, double &Y[], int incY);
void mql_dsbmv(CBLAS_ORDER order, CBLAS_UPLO Uplo, int N, int K, double alpha, double &A, int lda, double &X[], int incX, double beta, double &Y[], int incY);
void mql_dtbmv(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, int N, int K, double &A[], int lda, double &X[], int incX);
void mql_dtbsv(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, int N, int K, double &A[], int lda, double &X[], int incX);
void mql_dtpmv(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, int N, double &Ap[], double &X[], int incX);
void mql_dtpsv(CBLAS_ORDER order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, int N, double &Ap[], double &X[], int incX);
void mql_dsymv(CBLAS_ORDER order, CBLAS_UPLO Uplo, int N, double alpha, double &A[], int lda, double &X[], int incX, double beta, double &Y[], int incY);
void mql_dspmv(CBLAS_ORDER order, CBLAS_UPLO Uplo, int N, double alpha, double &Ap[], double &X[], int incX, double beta, double &Y[], int incY);
void mql_dspr(CBLAS_ORDER order, CBLAS_UPLO Uplo, int N, double alpha, double &X, int incX, double &Ap[]); 
void mql_dspr2(CBLAS_ORDER order, CBLAS_UPLO Uplo, int N, double alpha, double &X[], int incX, double &Y[], int incY, double &A[]);
void mql_dgemm(CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB, int M, int N, int K,  double alpha, double &A[], int lda, double &B[], int ldb, double beta, double &C[], int ldc);
void mql_dsymm(CBLAS_ORDER Order, CBLAS_SIDE Side, CBLAS_UPLO Uplo, int M, int N, double alpha, double &A[], int lda, double &B[], int ldb, double beta, double &C[], int ldc);
void mql_dsyrk(CBLAS_ORDER Order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans, int N, int K, double alpha, double &A[], int lda, double beta, double &C[], int ldc);
void mql_dsyr2k(CBLAS_ORDER Order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans, int N, int K, double alpha, double &A[], int lda, double &B[], int ldb, double beta, double &C[], int ldc);
void mql_dtrmm(CBLAS_ORDER Order, CBLAS_SIDE Side,  CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, int M, int N, double alpha, double &A[], int lda, double &B[], int ldb);
void mql_dtrsm(CBLAS_ORDER Order, CBLAS_SIDE Side, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, int M, int N, double alpha, double &A[], int lda, double &B[], int ldb);

int mql_dsyev(int matrix_order, char jobz, char uplo, int n, double &a[], int lda, double &w[]);
int mql_dgesvd(int matrix_order, char jobu, char jobvt, int m, int n, double &a[], int lda, double &s[], double &u[], int ldu, double &vt[], int ldvt, double &superb[]);
int mql_dgesdd(int matrix_order, char jobz, int m, int n, double &a[], int lda, double &s[], double &u[], int ldu, double &vt[], int ldvt);
int mql_dgels(int matrix_order, char trans, int m, int n, int nrhs, double &a[], int lda, double &b[], int ldb);
#import

