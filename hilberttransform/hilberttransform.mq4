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
#property indicator_buffers 2

#import "fftwinterface.dll"
void fastcosinetransform(double &a[], int dctlength, int inversefct);
void fastsinetransform(double &a[], int dctlength, int inversefct);
void hilberttransform(double &a[], double &b[], int length);
#import

extern double n = 10;
extern int cutoff_frequency_short = 128;	//cut-off frequency
extern int cutoff_frequency_medium = 32;

//extern double n = 12;
//extern int cutoff_frequency_short = 1024;       //cut-off frequency
//extern int cutoff_frequency_medium = 256;

double DCT_buffer[];
double DCT_buffer2[];

int dctlength;
bool ErrorFlag = FALSE;

double aa[];
double bb[];

#define LEN 5
double a[];
double b[];

#define pi 3.14159265358979323846

int init()
{
    Print("initialized");
    SetIndexStyle(0, DRAW_LINE, EMPTY, 2, Red);
    SetIndexBuffer(0, DCT_buffer);

    SetIndexStyle(1, DRAW_LINE, EMPTY, 2, Blue);
    SetIndexBuffer(1, DCT_buffer2);

    dctlength = MathPow(2, n);

    int M = ArrayResize(aa, dctlength + 1);
    if (M < dctlength) { Print("Cannot allocate memory", dctlength, M); ErrorFlag = TRUE; return (0); }
    SetIndexDrawBegin(0, Bars - M);

    ArrayResize(a, LEN);
    ArrayResize(b, LEN);
    for (int i=0; i<LEN;i++) a[i] = i+1;
    for (i=0; i<LEN;i++) Print("ORG : ", a[i]);
    fastcosinetransform(a, LEN, false);
    ArrayCopy(b, a);
    for (i=0; i<LEN;i++) Print("DCT : ", b[i]);
    hilberttransform(a, b, LEN);
    for (i=0; i<LEN;i++) Print("HIL : ", b[i]);

    return (0);
}

int deinit()
{
    return (0);
}

int start()
{
    int i;
    double phase;

    for (i = 0; i < dctlength; i++) {
	aa[i] = Close[i];
    }

/*
//for debugging   
    for (i = 0; i < dctlength; i++) {
	aa[i] = MathSin(10.0 * i * pi / 180);
    }
*/
    fastcosinetransform(aa, dctlength, false);
    ArrayCopy(bb, aa);

//bandpass filter by a rectangular window.
    for (i = 0; i < dctlength; i++) {
        if (i >= cutoff_frequency_short)  aa[i] = .0;
        if (i <= cutoff_frequency_medium) aa[i] = .0;
    }
    aa[0] = 0; //move origin to average of the period

//Hilbert transformation
    hilberttransform_mql(aa, bb, dctlength);

//reverse transformation
    fastcosinetransform(aa, dctlength, true);

//back transformation is not necessary for hilbert transformed one
    for (i = 0; i < dctlength; i++) {
	DCT_buffer[i] = aa[i];
    }
    for (i = 0; i < dctlength; i++) {
	DCT_buffer2[i] = bb[i];
    }
    return (0);
}

//textbook implementation of sine part of Hilbert transform via DCT
//Note: it scales like O(N^2)
void hilberttransform_mql(double &y[], double &x[], int length)
{
    double theta, tmp;
    for (int n = 0; n < length; n++) {
	tmp = 0.0;
	for (int k = 1; k < length; k++) {
	    theta = (pi * k * (2.0 * n + 1) / (2.0 * length));
	    tmp = tmp + (y[k] / length) * MathSin(theta);
	}
	x[n] = tmp;
    }
}

/*
void hilberttransform_mql(double &y[], double &x[], int length)
{
    double theta, tmp;
    y[0] = y[0] / MathSqrt(length) / 2.0;
    for (int n = 1; n < length; n++) {
        y[n] = y[n] / MathSqrt(2.0 * length);
    }

    double a[]; ArrayResize(a, length + 1);
    a[0] = 1.0 / MathSqrt(length);
    for (int i = 1; i < length; i++)
        a[i] = MathSqrt(2.0 / length);

    for (n = 0; n < length; n++) {
        tmp = 0.0;
        for (int k = 0; k < length; k++) {
            theta = (pi * k * (2.0 * n + 1) / (2.0 * length));
            tmp = tmp + (a[k] * y[k] * MathSin(theta));
        }
        x[n] = tmp;
    }
}
*/
