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
#include <lapacke.h>

#define _DLLAPI extern "C" __declspec(dllexport)

_DLLAPI void __stdcall BasicSSA(double *x, int N, int L, int Rmax, double *xtilde);
_DLLAPI void __stdcall BasicSSA_2(double *x, int N, int L, double threshold, double *xtilde);

// Do Basic Singular Spectrum Analysis for a time series.
// x is real-valued time series (x_1, x_2, \cdots, x_N) of length N.
// L is the window length, K = N - N + 1;
// xtilde is reconstructed series.
// this BasicSSA is not fast; do full SVD via LAPACK and do not employ Hankel structure of matrix.
_DLLAPI void __stdcall BasicSSA(double *x, int N, int L, int Rmax, double *xtilde)
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
    // (\lambda_i, Ui, V_i), 1 <= i <= R
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

// Do Basic Singular Spectrum Analysis for a time series.
// Maxmium number of singular vaule is determined by thershold,
// 90% of trace -> 1, 99% of trace -> 2, 99.9%  of trace -> 3 etc.
_DLLAPI void __stdcall BasicSSA_2(double *x, int N, int L, double threshold, double *xtilde)
{
    int K = N - L + 1;
    int i, j, k;
    int p, q;
    int Rmax;
    double *X = new double[L * K];
    double *Ystar = new double[L * K];
    double *U = new double[L * L];
    double *VT = new double[K * K];
    double *S = new double[L];
    double rtmp, trace;
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

    // LAPACK dgesdd (dgesvd) calculates its singular values in descending order
    // and determine Rmax which eigentriples to be concerned by threadshold.
    // (\lambda_i, Ui, V_i), 1 <= i <= rank
    // implicitly correspond to [eq. (2.3)]
    trace = 0.0;
    for (i = 1; i <= L; i++)
	trace = trace + S[i - 1];
    trace = trace * (1.0 - pow(10.0, -threshold));
    rtmp = 0.0;
    for (i = 1; i <= L; i++) {
	rtmp = rtmp + S[i - 1];
	Rmax = i;
	if (rtmp > trace)
	    break;
    }

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
    double *utilde = new double[N];
    double *vtilde = new double[N];
    double *w = new double[N];
    double *gprime = new double[N];

    std::complex <double>*uhat, *vhat, *ghat;
    uhat = (std::complex <double>*) fftw_malloc(sizeof(std::complex <double>) * N);
    vhat = (std::complex <double>*) fftw_malloc(sizeof(std::complex <double>) * N);
    ghat = (std::complex <double>*) fftw_malloc(sizeof(std::complex <double>) * N);
    int i, j;

    for (i = 1; i <= L; i++)	 utilde[i - 1] = u[i - 1];
    for (i = L + 1; i <= N; i++) utilde[i - 1] = 0.0;
    for (i = 1; i <= K; i++)  	 vtilde[i - 1] = v[i - 1];
    for (i = K + 1; i <= N; i++) vtilde[i - 1] = 0.0;

    fftw_plan x = fftw_plan_dft_r2c_1d(N, utilde, reinterpret_cast <fftw_complex *>(uhat), FFTW_ESTIMATE);
    fftw_plan y = fftw_plan_dft_r2c_1d(N, vtilde, reinterpret_cast <fftw_complex *>(vhat), FFTW_ESTIMATE);
    fftw_execute(x);
    fftw_execute(y);
    for (i = 1; i <= N; i++) ghat[i - 1] = vhat[i - 1] * uhat[i - 1] / double (N);
    fftw_plan z = fftw_plan_dft_c2r_1d(N, reinterpret_cast <fftw_complex *>(ghat), gprime, FFTW_ESTIMATE);
    fftw_execute(z);

    for (i = 1; i <= L; i++)                     w[i - 1] = 1.0 / double (i);
    for (i = L + 1; i <= N - L; i++)             w[i - 1] = 1.0 / double (L);
    for (i = N - L + 1, j = L; i <= N; i++, j--) w[i - 1] = 1.0 / double (j);
    for (i = 1; i <= N; i++) g[i - 1] = sigma * w[i - 1] * gprime[i - 1];

    fftw_destroy_plan(z);
    fftw_destroy_plan(y);
    fftw_destroy_plan(x);
    fftw_free(ghat);
    fftw_free(vhat);
    fftw_free(uhat);

    delete[]gprime;
    delete[]w;
    delete[]vtilde;
    delete[]utilde;
}

void Bidiagonalization(double *X, int N, int ldx, double *U, int L, int ldu, double *V, int K, int ldv, double *alpha, double *beta)
{
    double *p = new double[L];
    double *r = new double[K];
    double _a = 1.0, _b = 0.0;
    double rtmp;
    int i, j, k;

    for (i = 0; i < L * L; i++)	U[i] = 0.0;
    for (i = 0; i < K * K; i++)	V[i] = 0.0;

    //initial matrix
    for (i = 1; i <= L; i++) p[i - 1] = 1.0;
    beta[0] = cblas_dnrm2(L, p, 1);
    for (i = 1; i <= L; i++) U[(i - 1)] = p[i - 1] / beta[0];
    for (i = 1; i <= K; i++) v0[i - 1] = 0.0;
    for (j = 1; j <= L; j++) {
	cblas_dgemv(CblasColMajor, CblasTrans, L, K, _a, X, ldx, &U[(j - 1) * ldu], 1, _b, r, 1);
	if (j == 1) cblas_daxpy(K, -beta[j - 1], V, 1, r, 1); //actually zero vector should be input.
	else        cblas_daxpy(K, -beta[j - 1], &V[(j - 2) * ldv], 1, r, 1);
	alpha[j - 1] = cblas_dnrm2(K, r, 1);
	rtmp = 1.0 / alpha[j - 1];
	cblas_dcopy(K, r, 1, &V[(j - 1) * ldv], 1);
	cblas_dscal(K, rtmp, &V[(j - 1) * ldv], 1);
	cblas_dgemv(CblasColMajor, CblasNoTrans, L, K, _a, X, ldx, &V[(j - 1) * ldv], 1, _b, p, 1);
	cblas_daxpy(L, -alpha[j - 1], &U[(j - 1) * ldu], 1, p, 1);
	beta[j] = cblas_dnrm2(L, p, 1);
	rtmp = 1.0 / beta[j];
	cblas_dcopy(L, p, 1, &U[j * ldu], 1);
	cblas_dscal(L, rtmp, &U[j * ldu], 1);
    }
    delete []p;
    delete []r;
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
