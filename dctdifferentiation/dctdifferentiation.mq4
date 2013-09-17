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
#property link "https://github.com/nakatamaho"

#property indicator_separate_window
#property indicator_buffers 3

#import "fftwinterface.dll"
void fastcosinetransform(double &a[], int dctlength, int inversefct);
void fastsinetransform(double &a[], int dctlength, int inversefct);
void dct1stderivative(double &a[], int dctlength);
void dct2ndderivative(double &a[], int dctlength);
#import

extern double n = 12;
extern int cutoff_frequency_high = 1024;	//cut-off frequency
extern int cutoff_frequency_low = 256;	//cut-off frequency
double DCT_buffer[];
double deriv1_buffer[];
double deriv2_buffer[];
double tmp_buffer1[];
double tmp_buffer2[];
double tmp_buffer3[];
int dctlength;
bool ErrorFlag = FALSE;
int init()
{
    Print("initialized");
    SetIndexStyle(0, DRAW_LINE, EMPTY, 2, Red);
    SetIndexBuffer(0, DCT_buffer);
    SetIndexStyle(1, DRAW_LINE, EMPTY, 2, Blue);
    SetIndexBuffer(1, deriv1_buffer);
    SetIndexStyle(2, DRAW_LINE, EMPTY, 2, Green);
    SetIndexBuffer(2, deriv2_buffer);

    dctlength = MathPow(2, n);
    int M = ArrayResize(tmp_buffer1, dctlength);
    if (M < dctlength) {
	Print("Cannot allocate memory", dctlength, M);
	ErrorFlag = TRUE;
	return (0);
    }

    M = ArrayResize(tmp_buffer2, dctlength);
    if (M < dctlength) {
	Print("Cannot allocate memory", dctlength, M);
	ErrorFlag = TRUE;
	return (0);
    }

    M = ArrayResize(tmp_buffer3, dctlength);
    if (M < dctlength) {
	Print("Cannot allocate memory", dctlength, M);
	ErrorFlag = TRUE;
	return (0);
    }


    SetIndexDrawBegin(0, Bars - M);
    SetIndexDrawBegin(1, Bars - M);
    SetIndexDrawBegin(2, Bars - M);
    return (0);
}

int deinit()
{
    return (0);
}


//this function can be used to avoid transient phenomena
void shift_cyclic_array(double &a[], int length, int shift)
{
    if (shift == 0)
	return;
    int p;
    double tmp[];
    ArrayResize(tmp, length);
    for (p = 0; p < length; p++)
	tmp[p] = a[p];
    if (shift > 0) {
	for (p = shift; p < length; p++) {
	    a[p - shift] = tmp[p];
	}
	for (p = 0; p < shift; p++) {
	    a[length - shift + p] = tmp[p];
	}
    }
    if (shift < 0) {
	shift = -shift - 1;	//sorry, not a good coding, please fix!
	for (p = shift + 1; p < length; p++) {
	    a[p] = tmp[p - shift - 1];
	}
	for (p = 0; p <= shift; p++) {
	    a[p] = tmp[p + (length - shift) - 1];
	}
    }
}

int start()
{
    int i;
    for (i = dctlength - 1; i >= 0; i--) {
	tmp_buffer1[i] = Close[i];
    }

//the discrete cosine (DCT-II) transformation
    fastcosinetransform(tmp_buffer1, dctlength, false);

//cutoff by rectangular window
    for (i = 0; i < dctlength; i++) {
	if (i >= cutoff_frequency_high)
	    tmp_buffer1[i] = .0;
	if (i <= cutoff_frequency_low)
	    tmp_buffer1[i] = .0;
    }

    ArrayCopy(tmp_buffer2, tmp_buffer1);
    ArrayCopy(tmp_buffer3, tmp_buffer1);

//1st and 2nd derivative. No backtransformation is needed.
    dct1stderivative(tmp_buffer2, dctlength);
    dct2ndderivative(tmp_buffer3, dctlength);

    fastcosinetransform(tmp_buffer1, dctlength, true);
    for (i = 0; i <= dctlength - 1; i++) {
	DCT_buffer[i] = tmp_buffer1[i];
	deriv1_buffer[i] = tmp_buffer2[i];
	deriv2_buffer[i] = tmp_buffer3[i];
    }
    return (0);
}
