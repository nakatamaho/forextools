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
#property link "https://github.com/nakatamaho/filteroversampler"

//#property indicator_separate_window
#property indicator_chart_window
#property indicator_buffers 1

#import "fftwinterface.dll"
void fastcosinetransform(double &a[], int dctlength, int inversefct);
void fastsinetransform(double &a[], int dctlength, int inversefct);
void dct1stderivative(double &a[], int dctlength);
void dct2ndderivative(double &a[], int dctlength);

#import

//extern double n = 9;
//extern int cutoff_frequency_short = 128;      //cut-off frequency
//extern int cutoff_frequency_medium = 32;
extern double n = 10;
extern int cutoff_frequency_short = 1024;	//cut-off frequency
extern int cutoff_frequency_medium = 0;
extern double Gamma = 0.7;	//parameter for LagurreFilter
extern int oversampleperiod = PERIOD_M15;
int dctlength;
bool ErrorFlag = FALSE;
double DCT_buffer[];
double DCTderiv_buffer[];
double tmp_buffer[];
double tmp_buffer2[];
double LaguerreFilterBuffer[];
int init()
{
    Print("initialized");
    SetIndexStyle(0, DRAW_LINE, EMPTY, 2, Red);
    SetIndexBuffer(0, DCT_buffer);
    SetIndexStyle(1, DRAW_LINE, EMPTY, 2, Blue);
    SetIndexBuffer(1, DCTderiv_buffer);
    dctlength = MathPow(2, n);
    int M = ArrayResize(tmp_buffer, dctlength);
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
    M = ArrayResize(LaguerreFilterBuffer, dctlength);
    if (M < dctlength) {
	Print("Cannot allocate memory", dctlength, M);
	ErrorFlag = TRUE;
	return (0);
    }
    SetIndexDrawBegin(0, Bars - M);
    return (0);
}

int deinit()
{
    return (0);
}

int start()
{
    int i;
    LaguerreFilter(dctlength, Gamma, LaguerreFilterBuffer, oversampleperiod);
    for (i = dctlength - 1; i >= 0; i--) {
	tmp_buffer[i] = LaguerreFilterBuffer[i];
    }
    fastcosinetransform(tmp_buffer, dctlength, false);

//bandpass filter by a rectangular window.
    for (i = 0; i < dctlength; i++) {
	if (i >= cutoff_frequency_short)
	    tmp_buffer[i] = .0;
	if (i <= cutoff_frequency_medium)
	    tmp_buffer[i] = .0;
    }

//1st derivative
    ArrayCopy(tmp_buffer2, tmp_buffer);
    dct1stderivative(tmp_buffer2, dctlength);

//reverse transformation
    fastcosinetransform(tmp_buffer, dctlength, true);
    for (i = 0; i < dctlength; i++) {
	DCT_buffer[i] = tmp_buffer[i];
    }
    for (i = 0; i < dctlength; i++) {
	DCTderiv_buffer[i] = tmp_buffer2[i];
    }
    return (0);
}

double P(int index, int timeframe)
{
    return ((iHigh(NULL, timeframe, index) + iLow(NULL, timeframe, index)) / 2.0);
}

void LaguerreFilter(int length, double Gamma, double &LaguerreFilterBuffer[], int timeframe)
{
    double L0[];
    double L1[];
    double L2[];
    double L3[];
    double tmpfilterbuffer[];
    int s;
    int currentperiod = Period();
    int oversampleratio = currentperiod / timeframe;
    int oversamplelength = length * oversampleratio;
    int M = ArrayResize(L0, oversamplelength);
    M = ArrayResize(L1, oversamplelength);
    M = ArrayResize(L2, oversamplelength);
    M = ArrayResize(L3, oversamplelength);
    M = ArrayResize(tmpfilterbuffer, oversamplelength);
    for (s = oversamplelength; s >= 0; s--) {
	if (s > oversamplelength - 4) {
	    L0[s] = P(s, timeframe);
	    L1[s] = P(s, timeframe);
	    L2[s] = P(s, timeframe);
	    L3[s] = P(s, timeframe);
	    tmpfilterbuffer[s] = P(s, timeframe);
	}
	L0[s] = (1.0 - Gamma) * P(s, timeframe) + Gamma * L0[s + 1];
	L1[s] = -Gamma * L0[s] + L0[s + 1] + Gamma * L1[s + 1];
	L2[s] = -Gamma * L1[s] + L1[s + 1] + Gamma * L2[s + 1];
	L3[s] = -Gamma * L2[s] + L2[s + 1] + Gamma * L3[s + 1];
	tmpfilterbuffer[s] = (L0[s] + 2.0 * L1[s] + 2.0 * L2[s] + L3[s]) / 6.0;
    }
    for (s = 0; s < length; s++) {
	LaguerreFilterBuffer[s] = (tmpfilterbuffer[s * oversampleratio]
				   + 2.0 * tmpfilterbuffer[s * oversampleratio + 1]
				   + 2.0 * tmpfilterbuffer[s * oversampleratio + 2]
				   + tmpfilterbuffer[s * oversampleratio + 3]) / 6.0;
    }
}
