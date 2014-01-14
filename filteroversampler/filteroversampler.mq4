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

#property indicator_chart_window
#property indicator_buffers 2

extern double Gamma = 0.5;
extern int oversampleperiod = PERIOD_M15;
extern int length = 128;

bool ErrorFlag = FALSE;
double LaguerreFilterBuffer[];
double EhlersFilterBuffer[];

int init()
{
    Print("initialized");
    SetIndexStyle(0, DRAW_LINE, EMPTY, 4, Red);
    SetIndexBuffer(0, LaguerreFilterBuffer);
    SetIndexDrawBegin(0, length);

    SetIndexStyle(1, DRAW_LINE, EMPTY, 4, Blue);
    SetIndexBuffer(1, EhlersFilterBuffer);
    SetIndexDrawBegin(1, length);

    return (0);
}

int deinit()
{
    return (0);
}

int start()
{
    int i;
//    LaguerreFilter(length, Gamma, LaguerreFilterBuffer, oversampleperiod);
    EhlersFilter(length, EhlersFilterBuffer, oversampleperiod);
    return (0);
}

double P(int index, int timeframe)
{
    return ((iHigh(NULL, timeframe, index) + iLow(NULL, timeframe, index)) / 2.0);
}

void EhlersFilter(int period, double &EhlersFilterBuffer[], int timeframe)
{
    double Smooth[], SumCoef, Num;
    double Coef[], Distance2;
    double tmpfilterbuffer[];

    int M, s, count, lookback, length = 10;
    int currentperiod = Period();
    int oversampleratio = currentperiod / timeframe;
    int oversamplelength = period * oversampleratio;

    M = ArrayResize(Smooth, oversamplelength);
    M = ArrayResize(Coef, oversamplelength);
    M = ArrayResize(tmpfilterbuffer, oversamplelength);

    for (s = oversamplelength; s >= 0; s--) {
	Smooth[s] = (P(s, timeframe) + 2.0 * P(s, timeframe) + 2.0 * P(s, timeframe) + P(s, timeframe)) / 6.0;
    }
    for (s = oversamplelength; s >= 0; s--) {
	Distance2 = 0.0;
	for (lookback = 0; lookback < length; lookback++) {
	    Distance2 = Distance2 + (Smooth[s] - Smooth[s + lookback]) * (Smooth[s] - Smooth[s + lookback]);
	}
	Coef[s] = Distance2;

	Num = 0.0;
	SumCoef = 0.0;
	for (count = 0; count < length; count++) {
	    Num = Num + Coef[count + s] * Smooth[count + s];
	    SumCoef = SumCoef + Coef[count + s];
	}
	if (MathAbs(SumCoef) != 0.0) {
	    tmpfilterbuffer[s] = Num / SumCoef;
	}
    }
    for (s = 0; s < period; s++) {
	EhlersFilterBuffer[s] = (tmpfilterbuffer[s * oversampleratio]
				   + 2.0 * tmpfilterbuffer[s * oversampleratio + 1]
				   + 2.0 * tmpfilterbuffer[s * oversampleratio + 2]
				   + tmpfilterbuffer[s * oversampleratio + 3]) / 6.0;
    }
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
