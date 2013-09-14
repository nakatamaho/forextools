/*
Copyright 2013 NAKATA Maho. All rights reserved.

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
#include <fftw3.h>

//----
#define MT4_EXPFUNC __declspec(dllexport)

BOOL APIENTRY DllMain(HANDLE hModule, DWORD ul_reason_for_call, LPVOID lpReserved)
{
//----
    switch (ul_reason_for_call) {
    case DLL_PROCESS_ATTACH:
    case DLL_THREAD_ATTACH:
    case DLL_THREAD_DETACH:
    case DLL_PROCESS_DETACH:
	break;
    }
//----
    return (TRUE);
}

#ifdef __cplusplus
extern "C"
{
#endif

MT4_EXPFUNC void __stdcall fastcosinetransform(double *a, int length, BOOL inversefct)
{
	fftw_plan p;
	double *in;
	int i;

	 in = (double *) fftw_malloc(sizeof(double) * length);
	for (i = 0; i < length; i++)
	     in[i] = a[i];

	//Note: you can use DCT-I and IV but may break consistency
	//with fastsinetransform, dct1stderivative, dct2ndderivative
	//and hilberttransformation.
	if (inversefct)
	     p = fftw_plan_r2r_1d(length, in, in, FFTW_REDFT01, FFTW_ESTIMATE);
	else
	     p = fftw_plan_r2r_1d(length, in, in, FFTW_REDFT10, FFTW_ESTIMATE);

	 fftw_execute(p);
	if (inversefct)
	{
	    for (i = 0; i < length; i++)
		a[i] = in[i] / length / 2.0;
	} else
	{
	    for (i = 0; i < length; i++)
		a[i] = in[i];
	}
	fftw_destroy_plan(p);
	fftw_free(in);
	return;
}

MT4_EXPFUNC void __stdcall fastsinetransform(double *a, int length, BOOL inversefct)
{
	fftw_plan p;
	double *in;
	int i;
	in = (double *) fftw_malloc(sizeof(double) * length);
	for (i = 0; i < length; i++)
	    in[i] = a[i];

	//DST-II and IDST-II (DST-III)
	if (inversefct)
	    p = fftw_plan_r2r_1d(length, in, in, FFTW_RODFT01, FFTW_ESTIMATE);
	else
	    p = fftw_plan_r2r_1d(length, in, in, FFTW_RODFT10, FFTW_ESTIMATE);
	fftw_execute(p);
	if (inversefct) {
	    for (i = 0; i < length; i++)
		a[i] = in[i] / length / 2.0;
	} else {
	    for (i = 0; i < length; i++)
		a[i] = in[i];
	}
	fftw_destroy_plan(p);
	fftw_free(in);
	return;
}

MT4_EXPFUNC void __stdcall dct1stderivative(double *a, int length)
{
	fftw_plan p;
	double *in;
	int i;
	in = (double *) fftw_malloc(sizeof(double) * length);
	for (i = 0; i < length - 1; i++)
	    in[i] = -a[i + 1] * M_PI * (i + 1) / length;
	in[length - 1] = 0.0;
	p = fftw_plan_r2r_1d(length, in, in, FFTW_RODFT01, FFTW_ESTIMATE);
	fftw_execute(p);
	for (i = 0; i < length; i++)
	    a[i] = in[i] / length / 2.0;
	fftw_destroy_plan(p);
	fftw_free(in);
	return;
}

MT4_EXPFUNC void __stdcall dct2ndderivative(double *a, int length)
{
	fftw_plan p;
	double *in;
	int i;
	in = (double *) fftw_malloc(sizeof(double) * length);
	for (i = 0; i < length; i++)
	    in[i] = -a[i] * (M_PI * i / length) * (M_PI * i / length);
	in[0] = 0.0;
	p = fftw_plan_r2r_1d(length, in, in, FFTW_REDFT01, FFTW_ESTIMATE);
	fftw_execute(p);
	for (i = 0; i < length; i++)
	    a[i] = in[i] / length / 2.0;
	fftw_destroy_plan(p);
	fftw_free(in);
	return;
}

MT4_EXPFUNC void __stdcall hilberttransform(double *a, double *b, int length)
{
	for (int i = 0; i < length - 1; i++)
	    b[i] = a[i + 1];
	b[length - 1] = 0.0;
	fastsinetransform(b, length, true);
}

#ifdef __cplusplus
}
#endif
