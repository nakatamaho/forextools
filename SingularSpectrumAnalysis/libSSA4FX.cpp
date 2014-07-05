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

//Nina Golyandina, Anatoly Zhigljavsky "Singular Spectrum Analysis for Time Series"
//Springer 2013
//DOI 10.1007/978-3-642-34913-3 ISBN 978-3-642-34913-3
//Chapter 2 can be downloaded for free
//http://www.springer.com/cda/content/document/cda_downloaddocument/9783642349126-c2.pdf

//Other documents
//http://www.gistatgroup.com/cat/book3/index.html
//http://cran.r-project.org/web/packages/Rssa/
//http://ssa.cf.ac.uk/a_brief_introduction_to_ssa.pdf
//http://www.jds-online.com/file_download/133/JDS-396.pdf
//http://arxiv.org/pdf/0911.4498v2.pdf

#define WIN32_LEAN_AND_MEAN	// Exclude rarely-used stuff from Windows headers
#include <windows.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <complex>
#include <fftw3.h>
#include <cblas.h>
#include <lapacke.h>

#define F77_FUNC(name,NAME) name ## _
#define F77_FUNC_(name,NAME) name ## _
#define dlansvd_irl_f77 F77_FUNC (dlansvd_irl, DLANDSVD_IRL)
#define dlansvd_f77 F77_FUNC (dlansvd, DLANDSVD)

typedef void (*daprod) (char *transa, int *m, int *n, double *x, double *y, double *dparm, int *iparm);
extern "C" int dlansvd_irl_f77(char *which, char *jobu, char *jobv, int *m, int *n, int *dim, int *p, int *neig, int *maxiter, daprod aprod, double *u, int *ldu, double *sigma, double *bnd, double *v, int *ldv, double *tolin, double *work, int *lwork, int *iwork, int *liwork, double *doption, int *option, int *info, double *dparm, int *iparm);
extern "C" int dlansvd_f77(char *jobu, char *jobv, int *m, int *n, int *k, int *kmax, daprod aprod, double *u, int *ldu, double *sigma, double *bnd, double *v, int *ldv, double *tolin, double *work, int *lwork, int *iwork, int *liwork, double *doption, int *ioption, int *info, double *dparm, int *iparm);

#define _DLLAPI extern "C" __declspec(dllexport)

_DLLAPI void __stdcall BasicSSA(double *x, int N, int L, int Rmax, double *xtilde);
_DLLAPI void __stdcall BasicSSA_LAPACK(double *x, int N, int L, int Rmax, double *xtilde);
//for compatibility
_DLLAPI void __stdcall fastsingular(double *x, int N, int L, int Rmax, double *xtilde);
_DLLAPI void __stdcall fastSingular(double *x, int N, int L, int Rmax, double *xtilde);

void FastSSAMatVecMult(double *F, int N, int L, double *v, double *p);
void FastSSAMatTransVecMult(double *F, int N, int L, double *v, double *p);
void RankOneHankelization(double *u, int L, double *v, int K, double sigma, double *g, int N);

int Mlsame(const char *a, const char *b)
{
    if (toupper(*a) == toupper(*b))
	return 1;
    return 0;
}

void matvecprod_dense(char *transa, int *m, int *n, double *x, double *y, double *dparm, int *iparm)
{
    int i, j;
    int ldx = iparm[0];
    double rtmp;
    double alpha = 1.0, beta = 0.0;
    if (Mlsame(transa, "N")) {
	cblas_dgemv(CblasColMajor, CblasNoTrans, (*m), (*n), alpha, dparm, ldx, x, 1, beta, y, 1);
	return;
    }
    if (Mlsame(transa, "T")) {
	cblas_dgemv(CblasColMajor, CblasTrans, (*m), (*n), alpha, dparm, ldx, x, 1, beta, y, 1);
	return;
    }
}

// Do Basic Singular Spectrum Analysis for a time series using PROPACK
// x is real-valued time series (x_1, x_2, \cdots, x_N) of length N.
// L is the window length, K = L - N + 1;
// xtilde is reconstructed series.
_DLLAPI void __stdcall BasicSSA(double *x, int N, int L, int Rmax, double *xtilde)
{
    int K = N - L + 1;
    int i, j, k, p, q;

    int ldx = L, ldystar = L, ldu = L, ldv = K;
    double tolin = 1e-9;
    int nb = 128;		//block
    int kmax = K * 3;		//maximum itration number
    int lwork = L + K + 9 * kmax + 5 * kmax * kmax + 4 + std::max(3 * kmax * kmax + 4 * kmax + 4, nb * std::max(L, K));

    double *work = new double[lwork];
    int liwork = 8 * kmax;
    int *iwork = new int[liwork];
    int info;
    int *iparm = new int[1];
    iparm[0] = L;		//lda;
    double *doption = new double[3];
    int *ioption = new int[2];
    double eps = 1e-15;		//dlamch('e')
    char jobu = 'Y';		//left singular vectors
    char jobv = 'Y';		//right singular vectors

    doption[0] = sqrt(eps);
    doption[1] = pow(eps, 0.75);
    doption[2] = 0.0;
    ioption[0] = 0;
    ioption[1] = 1;

    double *X = new double[K * L];
    double *Ystar = new double[K * L];
    double *U = new double[L * kmax];
    double *V = new double[K * kmax];
    double *S = new double[Rmax];
    double *bnd = new double[Rmax];
    double rtmp;

    //construct L-trajectory matrix [eq. (2.1)]
    for (j = 1; j <= K; j++) {
	for (i = 1; i <= L; i++) {
	    X[(i - 1) + (j - 1) * ldx] = x[(i - 1) + (j - 1)];
	}
    }

    dlansvd_f77(&jobu, &jobv, &L, &K, &Rmax, &kmax, matvecprod_dense, U, &ldu, S, bnd, V, &ldv, &tolin, work, &lwork, iwork, &liwork, doption, ioption, &info, X, iparm);

    RankOneHankelization(U, L, V, K, S[0], xtilde, N);
    double *g = new double[N];
    for (i = 2; i <= Rmax; i++) {
	RankOneHankelization(&U[ldu * (i - 1)], L, &V[ldv * (i - 1)], K, S[i - 1], g, N);
	cblas_daxpy(N, 1.0, g, 1, xtilde, 1);
    }
    delete[]g;
    delete[]ioption;
    delete[]doption;
    delete[]iparm;
    delete[]iwork;
    delete[]work;
    delete[]bnd;
    delete[]S;
    delete[]V;
    delete[]U;
    delete[]Ystar;
    delete[]X;
}

void Embedding(double *F, int N, int L, double *X, int ldx)
{
    int i, j;
    int K = N - L + 1;
    for (j = 1; j <= K; j++) {
	for (i = 1; i <= L; i++) {
	    X[(i - 1) + (j - 1) * ldx] = F[(i - 1) + (j - 1)];
	}
    }
}

void FastSSAMatVecMult(double *F, int N, int L, double *v, double *p)
{
    double *c = new double[N];
    double *w = new double[N];
    double *pp = new double[N];
    std::complex <double> *ctilde, *wtilde, *ptilde;
    ctilde = (std::complex <double>*) fftw_malloc(sizeof(std::complex <double>) * N);
    wtilde = (std::complex <double>*) fftw_malloc(sizeof(std::complex <double>) * N);
    ptilde = (std::complex <double>*) fftw_malloc(sizeof(std::complex <double>) * N);
    int i, j;

    for (i = N - L + 1; i <= N; i++) c[i - (N - L + 1)] = F[i - 1];
    for (i = 1; i <= N - L; i++) c[i + (L - 1)] = F[i - 1];
    for (i = 1; i <= N; i++) w[i - 1] = 0.0;
    for (i = N - L + 1, j = 1; i >= 1; i--, j++) w[j - 1] = v[i - 1];

    fftw_plan x = fftw_plan_dft_r2c_1d(N, c, reinterpret_cast <fftw_complex *>(ctilde), FFTW_ESTIMATE);
    fftw_plan y = fftw_plan_dft_r2c_1d(N, w, reinterpret_cast <fftw_complex *>(wtilde), FFTW_ESTIMATE);
    fftw_execute(x);
    fftw_execute(y);
    for (i = 1; i <= N; i++) ptilde[i - 1] = ctilde[i - 1] * wtilde[i - 1];
    fftw_plan z = fftw_plan_dft_c2r_1d(N, reinterpret_cast <fftw_complex *>(ptilde), pp, FFTW_ESTIMATE);
    fftw_execute(z);
    for (i = 1; i <= L; i++) p[i - 1] = pp[i - 1] / N;

    fftw_destroy_plan(z);
    fftw_destroy_plan(y);
    fftw_destroy_plan(x);
    fftw_free(ptilde);
    fftw_free(wtilde);
    fftw_free(ctilde);
    delete[]pp;
    delete[]w;
    delete[]c;
}

void FastSSAMatTransVecMult(double *F, int N, int L, double *v, double *p)
{
    double *c = new double[N];
    double *w = new double[N];
    double *pp = new double[N];
    std::complex <double>*ctilde, *wtilde, *ptilde;
    ctilde = (std::complex <double>*) fftw_malloc(sizeof(std::complex <double>) * N);
    wtilde = (std::complex <double>*) fftw_malloc(sizeof(std::complex <double>) * N);
    ptilde = (std::complex <double>*) fftw_malloc(sizeof(std::complex <double>) * N);
    int i, j;
    int K = N - L + 1;

    for (i = N - L + 1; i <= N; i++) c[i - (N - L) - 1] = F[i - 1];
    for (i = 1; i <= N - L; i++) c[i + L - 1] = F[i - 1];
    for (i = 1; i <= N - L; i++) w[i - 1] = 0.0;
    for (i = N - L + 1, j = L; i <= N; i++, j--) w[i - 1] = v[j - 1];

    fftw_plan x = fftw_plan_dft_r2c_1d(N, c, reinterpret_cast <fftw_complex *>(ctilde), FFTW_ESTIMATE);
    fftw_plan y = fftw_plan_dft_r2c_1d(N, w, reinterpret_cast <fftw_complex *>(wtilde), FFTW_ESTIMATE);
    fftw_execute(x);
    fftw_execute(y);
    for (i = 1; i <= N; i++) ptilde[i - 1] = ctilde[i - 1] * wtilde[i - 1];
    fftw_plan z = fftw_plan_dft_c2r_1d(N, reinterpret_cast <fftw_complex *>(ptilde), pp, FFTW_ESTIMATE);
    fftw_execute(z);

    for (i = 1; i <= K; i++) p[i - 1] = pp[i + L - 2] / N;

    fftw_destroy_plan(z);
    fftw_destroy_plan(y);
    fftw_destroy_plan(x);
    fftw_free(ptilde);
    fftw_free(wtilde);
    fftw_free(ctilde);

    delete[]pp;
    delete[]w;
    delete[]c;
}

void RankOneHankelization(double *u, int L, double *v, int K, double sigma, double *g, int N)
{
    double *utilde, *vtilde, *w, *gprime;
    utilde = (double *) fftw_malloc(sizeof(double) * N);
    vtilde = (double *) fftw_malloc(sizeof(double) * N);
    w = (double *) fftw_malloc(sizeof(double) * N);
    gprime = (double *) fftw_malloc(sizeof(double) * N);

    std::complex <double>*uhat, *vhat, *ghat;
    uhat = (std::complex <double>*) fftw_malloc(sizeof(std::complex <double>) * N);
    vhat = (std::complex <double>*) fftw_malloc(sizeof(std::complex <double>) * N);
    ghat = (std::complex <double>*) fftw_malloc(sizeof(std::complex <double>) * N);
    int i, j;

    for (i = 1; i <= L; i++) utilde[i - 1] = u[i - 1];
    for (i = L + 1; i <= N; i++) utilde[i - 1] = 0.0;
    for (i = 1; i <= K; i++) vtilde[i - 1] = v[i - 1];
    for (i = K + 1; i <= N; i++) vtilde[i - 1] = 0.0;

    fftw_plan x = fftw_plan_dft_r2c_1d(N, utilde, reinterpret_cast <fftw_complex *>(uhat), FFTW_ESTIMATE | FFTW_UNALIGNED);
    fftw_plan y = fftw_plan_dft_r2c_1d(N, vtilde, reinterpret_cast <fftw_complex *>(vhat), FFTW_ESTIMATE | FFTW_UNALIGNED);
    fftw_execute(x);
    fftw_execute(y);
    for (i = 1; i <= N; i++) ghat[i - 1] = vhat[i - 1] * uhat[i - 1] / double (N);
    fftw_plan z = fftw_plan_dft_c2r_1d(N, reinterpret_cast <fftw_complex *>(ghat), gprime, FFTW_ESTIMATE | FFTW_UNALIGNED);
    fftw_execute(z);

    for (i = 1; i <= L; i++) w[i - 1] = 1.0 / double (i);
    for (i = L + 1; i <= N - L; i++) w[i - 1] = 1.0 / double (L);
    for (i = N - L + 1, j = L; i <= N; i++, j--) w[i - 1] = 1.0 / double (j);
    for (i = 1; i <= N; i++) g[i - 1] = sigma * w[i - 1] * gprime[i - 1];

    fftw_destroy_plan(z);
    fftw_destroy_plan(y);
    fftw_destroy_plan(x);
    fftw_free(ghat);
    fftw_free(vhat);
    fftw_free(uhat);

    fftw_free(gprime);
    fftw_free(w);
    fftw_free(vtilde);
    fftw_free(utilde);
}

void Hankelization(double *Y, int ldy, int N, int L, double *ytilde)
{
    int K = N - L + 1;
    int j, k;
    double rtmp;
// Hankelization (or diagonal averaging) [eq. (2.4)]
    for (k = 1; k <= L; k++) {
	rtmp = 0.0;
	for (j = 1; j <= k; j++) {
	    rtmp = rtmp + Y[j - 1 + (k - j) * ldy];
	}
	ytilde[k - 1] = rtmp / double (k);
    }
    for (k = L; k < K; k++) {
	rtmp = 0.0;
	for (j = 1; j <= L; j++) {
	    rtmp = rtmp + Y[j - 1 + (k - j) * ldy];
	}
	ytilde[k - 1] = rtmp / double (L);
    }
    for (k = K; k <= N; k++) {
	rtmp = 0.0;
	for (j = k - K + 1; j <= L; j++) {
	    rtmp = rtmp + Y[j - 1 + (k - j) * ldy];
	}
	ytilde[k - 1] = rtmp / double (N - k + 1);
    }
}

_DLLAPI void __stdcall fastsingular(double *x, int N, int L, int Rmax, double *xtilde)
{
    BasicSSA(x, N, L, Rmax, xtilde);
}

_DLLAPI void __stdcall fastSingular(double *x, int N, int L, int Rmax, double *xtilde)
{
    BasicSSA(x, N, L, Rmax, xtilde);
}

// Do Basic Singular Spectrum Analysis for a time series using LAPACK
// x is real-valued time series (x_1, x_2, \cdots, x_N) of length N.
// L is the window length, K = N - L + 1;
// xtilde is reconstructed series.
// this BasicSSA is not fast; do full SVD via LAPACK and do not employ Hankel structure of matrix.
_DLLAPI void __stdcall BasicSSA_LAPACK(double *x, int N, int L, int Rmax, double *xtilde)
{
    int K = N - L + 1;
    int i, j, k;
    int p, q;
    double *X = new double[L * K];
    double *Ystar = new double[L * K];
    double *U = new double[L * L];
    double *VT = new double[K * K];
    double *S = new double[L];
    double rtmp;
    int ldx = L, ldystar = L, ldu = L, ldvt = K;

    //construct L-trajectory matrix [eq. (2.1)]
    for (j = 1; j <= K; j++) {
	for (i = 1; i <= L; i++) {
	    X[(i - 1) + (j - 1) * ldx] = x[(i - 1) + (j - 1)];
	}
    }

    // Do singular value decomposition of trajectory matrix.
    // implicitly correspond to [eq. (2.2)]
    LAPACKE_dgesdd(LAPACK_COL_MAJOR, 'A', L, K, X, L, S, U, L, VT, K);

    // Reconstruction using eigentriples
    // (\lambda_i, Ui, V_i), 1 <= i <= Rmax
    for (q = 1; q <= K; q++) {
	for (p = 1; p <= L; p++) {
	    rtmp = 0.0;
	    for (i = 1; i <= Rmax; i++) {
		rtmp = rtmp + U[(p - 1) + (i - 1) * ldu] * S[i - 1] * VT[(i - 1) + (q - 1) * ldvt];
	    }
	    Ystar[(p - 1) + (q - 1) * ldystar] = rtmp;
	}
    }

    // Hankelization (or diagonal averaging) [eq. (2.4)]
    for (k = 1; k <= L; k++) {
	rtmp = 0.0;
	for (j = 1; j <= k; j++) {
	    rtmp = rtmp + Ystar[j - 1 + (k - j) * ldystar];
	}
	xtilde[k - 1] = rtmp / double (k);
    }
    for (k = L; k < K; k++) {
	rtmp = 0.0;
	for (j = 1; j <= L; j++) {
	    rtmp = rtmp + Ystar[j - 1 + (k - j) * ldystar];
	}
	xtilde[k - 1] = rtmp / double (L);
    }
    for (k = K; k <= N; k++) {
	rtmp = 0.0;
	for (j = k - K + 1; j <= L; j++) {
	    rtmp = rtmp + Ystar[j - 1 + (k - j) * ldystar];
	}
	xtilde[k - 1] = rtmp / double (N - k + 1);
    }

    delete[]S;
    delete[]VT;
    delete[]U;
    delete[]Ystar;
    delete[]X;
}
