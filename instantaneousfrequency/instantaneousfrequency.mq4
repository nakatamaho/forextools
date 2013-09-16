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

#property copyright "NAKATA Maho"
#property link "https://github.com/nakatamaho/hilberttransform"

#property indicator_separate_window
#property indicator_buffers 3

#import "fftwinterface.dll"
void fastcosinetransform(double &a[], int dctlength, int inversefct);
void fastsinetransform(double &a[], int dctlength, int inversefct);
void hilberttransform(double &a[], double &b[], int length);
void dct1stderivative(double &a[], int length);
#import

//extern double n = 10;
//extern int cutoff_frequency_short = 128;      //cut-off frequency
//extern int cutoff_frequency_medium = 32;

extern double n = 12;
extern int cutoff_frequency_short = 1024;	//cut-off frequency
extern int cutoff_frequency_medium = 256;

double DCT_buffer[];
double DCT_buffer2[];
double freq_buffer[];

int dctlength;
bool ErrorFlag = FALSE;

double aa[];
double bb[];
double cc[];
double dd[];
double ee[];

#define pi 3.14159265358979323846

int init()
{
    Print("initialized");
    SetIndexStyle(0, DRAW_LINE, EMPTY, 2, Red);
    SetIndexBuffer(0, DCT_buffer);

    SetIndexStyle(1, DRAW_LINE, EMPTY, 2, Blue);
    SetIndexBuffer(1, DCT_buffer2);

    SetIndexStyle(2, DRAW_LINE, EMPTY, 2, Green);
    SetIndexBuffer(2, freq_buffer);

    dctlength = MathPow(2, n);

    int M = ArrayResize(aa, dctlength + 1);
    if (M < dctlength) {
	Print("Cannot allocate memory", dctlength, M);
	ErrorFlag = TRUE;
	return (0);
    }
    M = ArrayResize(bb, dctlength + 1);
    M = ArrayResize(cc, dctlength + 1);
    M = ArrayResize(dd, dctlength + 1);
    M = ArrayResize(ee, dctlength + 1);
    SetIndexDrawBegin(0, Bars - M);
    SetIndexDrawBegin(1, Bars - M);
    SetIndexDrawBegin(2, Bars - M);
    return (0);
}

double instantaneousfrequency(double &realpart[], double &imagpart[], double &realpartd[], double &imagpartd[], int shift)
{
    double x, v, denominator;
    double xdot, vdot, numerator;
    double instfreq;

    x = realpart[shift];
    v = imagpart[shift];
    xdot = realpartd[shift];
    vdot = imagpartd[shift];
    numerator = x * vdot - v * xdot;
    denominator = x * x + v * v;
    if (denominator != 0.)
	instfreq = numerator / denominator;
    else
	instfreq = 0.0;
    return (instfreq);
}

int deinit()
{
    return (0);
}

int start()
{
    int i, j;
    double phase;

    for (i = 0; i < dctlength; i++) {
	aa[i] = Close[i];
    }
    fastcosinetransform(aa, dctlength, false);

/*
//for debugging   
    for (i = 0; i < dctlength; i++) {
        if (i > 512 && i < 1024 ) aa[i] = MathSin(2 * i * pi * 2 * 64 / dctlength);
	else aa[i] = 0.0;
    }
*/

//bandpass filter by a rectangular window.
    for (i = 0; i < dctlength; i++) {
	if (i >= cutoff_frequency_short)
	    aa[i] = .0;
	if (i <= cutoff_frequency_medium)
	    aa[i] = .0;
    }
    aa[0] = 0;			//move origin to average of the period

    ArrayCopy(bb, aa);

//Hilbert transformation
    hilberttransform(aa, bb, dctlength);

//1st derivative of real part of analytic signal via dct
    ArrayCopy(cc, aa);
    dct1stderivative(cc, dctlength);

//1st derivative of imaginary part of analytic signal via dct
    ArrayCopy(dd, bb);
    fastcosinetransform(dd, dctlength, false);
    dct1stderivative(dd, dctlength);

//reverse transformation
    fastcosinetransform(aa, dctlength, true);

    for (i = 0; i < dctlength - 1; i++) {
	ee[i] = 0.0;
    }
    for (i = 0; i < dctlength - 1; i++) {
	ee[i] = instantaneousfrequency(aa, bb, cc, dd, i);
    }

//back transformation is not necessary for hilbert transformed one
    for (i = 0; i < dctlength; i++) {
	DCT_buffer[i] = 0.0;
    }
    for (i = 0; i < dctlength; i++) {
	DCT_buffer2[i] = 0.0;
    }
    for (i = 0; i < dctlength; i++) {
	freq_buffer[i] = ee[i];
    }

    return (0);
}
