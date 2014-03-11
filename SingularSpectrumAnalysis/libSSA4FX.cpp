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

#include <lapacke.h>

#define _DLLAPI extern "C" __declspec(dllexport)

_DLLAPI void __stdcall BasicSSA(double *x, int N, int L, int Rmax, double *xtilde);

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
    double rtmp, trace;
    int ldx = L, ldystar = L, ldxt = K, ldu = L, ldvt = K, ret;

    //construct L-trajectory matrix [eq. (2.1)]
    for (j = 1; j <= K; j++) {
	for (i = 1; i <= L; i++) {
	    X[(i - 1) + (j - 1) * ldx] = x[(i - 1) + (j - 1)];
	}
    }

    // Do singular value decomposition of trajectory matrix.
    // implicitly correspond to [eq. (2.2)]
    ret = LAPACKE_dgesdd(LAPACK_COL_MAJOR, 'A', L, K, X, L, S, U, L, VT, K);
    trace = 0.0;

/*
    // LAPACK dgesdd (dgesvd) calculates its singular values in descending order
    // and determine Lmax which eigentriples to be concerned by threadshold.
    // (\lambda_i, Ui, V_i), 1 <= i <= rank
    // implicitly correspond to [eq. (2.3)]
    for (i = 1; i <= L; i++) trace = trace + S[i - 1];
    trace = trace * (1.0 - pow(10.0, -threshold));
    rtmp = 0.0;
    for (i = 1; i <= L; i++) {
	rtmp = rtmp + S[i - 1];
	Rmax = i;
	if (rtmp > trace)
	    break;
    }
*/

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
