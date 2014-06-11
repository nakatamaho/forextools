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

#property copyright   "NAKATA Maho 2014"
#property link        "https://github.com/nakatamaho/forextools"
#property description "Casual SSA (Last-Point SSA)"
#property strict

#property indicator_chart_window
#property indicator_buffers 3
#property indicator_color1 DeepPink
#property indicator_color2 Blue
#property indicator_color3 Red
#property indicator_style1 STYLE_SOLID
#property indicator_style2 STYLE_SOLID
#property indicator_style3 STYLE_SOLID
#property indicator_width1 3
#property indicator_width2 4
#property indicator_width3 4

#import "libSSA4FX.dll"
void BasicSSA(double &x[], int N, int L, int Rmax, double &xtilde[]);
void BasicSSA_2(double &x[], int N, int L, double threshold, double &xtilde[]);
#import

input int TotalLength = 300;
input int WindowLength = 30;
input int Rmax = 3;
input double threshold = 2.0;
input int SSAMethod = 1;

double SSABuffer[];
double ExtBuffer[];
double PriceBuffer[];
double UpLine[];
double DnLine[];

int OnInit(void)
{
    IndicatorBuffers(3);
    IndicatorDigits(Digits);

    SetIndexStyle(0, DRAW_LINE);
    SetIndexBuffer(0, ExtBuffer);
    SetIndexShift(0, 0);
    SetIndexLabel(0, "Basic SSA");
    SetIndexDrawBegin(0, TotalLength);

    SetIndexStyle(1, DRAW_LINE);
    SetIndexBuffer(1, UpLine);
    SetIndexShift(1, 0);
    SetIndexDrawBegin(1, TotalLength);

    SetIndexStyle(2, DRAW_LINE);
    SetIndexBuffer(2, DnLine);
    SetIndexShift(2, 0);
    SetIndexDrawBegin(2, TotalLength);

    return (INIT_SUCCEEDED);
}

int OnCalculate(const int rates_total, const int prev_calculated, const datetime & time[],
		const double &open[], const double &high[], const double &low[], const double &close[], const long &tick_volume[],
		const long &volume[], const int &spread[])
{
    int i, j, k;
    int limit = prev_calculated;
    if (rates_total <= TotalLength * 2)
	return (0);

    if(prev_calculated > 0) limit--;

    ArraySetAsSeries(ExtBuffer, true);
    ArraySetAsSeries(SSABuffer, true);
    ArraySetAsSeries(PriceBuffer, true);
    ArraySetAsSeries(UpLine, true);
    ArraySetAsSeries(DnLine, true);
    ArraySetAsSeries(close, true);

    ArrayResize(SSABuffer, TotalLength);
    ArrayResize(PriceBuffer, TotalLength);

    for (i = TotalLength - 1 ; i > limit ; i--) {
        k = TotalLength - 1;
	for (j = i + TotalLength - 1 ; j >= i ; j--) {
	    PriceBuffer[k] = close[j];
            k--; 
	}
	if (SSAMethod)
	    BasicSSA(PriceBuffer, TotalLength, WindowLength, Rmax, SSABuffer);
	else
	    BasicSSA_2(PriceBuffer, TotalLength, WindowLength, threshold, SSABuffer);
	ExtBuffer[i] = SSABuffer[0];
    }

    for (i = 0; i <= TotalLength; i++) {
	if (ExtBuffer[i] >= ExtBuffer[i + 1]) {
	    UpLine[i] = ExtBuffer[i];
	    DnLine[i] = EMPTY_VALUE;
	} else {
	    UpLine[i] = EMPTY_VALUE;
	    DnLine[i] = ExtBuffer[i];
	}
	if (UpLine[i] != EMPTY_VALUE && UpLine[i + 1] == EMPTY_VALUE)
	    UpLine[i + 1] = ExtBuffer[i + 1];
	if (DnLine[i] != EMPTY_VALUE && DnLine[i + 1] == EMPTY_VALUE)
	    DnLine[i + 1] = ExtBuffer[i + 1];
    }
    return (rates_total);
}
