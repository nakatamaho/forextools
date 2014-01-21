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

#property indicator_chart_window
#property indicator_buffers 1
#property indicator_color1 DeepPink
#property indicator_style1 STYLE_SOLID
#property indicator_width1 3
double ExtMapBuffer[];
double array[][6];

int init()
{
    SetIndexStyle(0, DRAW_LINE);
    SetIndexBuffer(0, ExtMapBuffer);
    return (0);
}

int start()
{
    int counted_bars = IndicatorCounted();
    if (counted_bars > 0) counted_bars--;
    int limit = Bars - counted_bars;
    OverSamplingMA(0, limit - 1, ExtMapBuffer, Symbol(), PERIOD_M5);
    return (0);
}

void OverSamplingMA(int start, int end, double &retbuffer[], string symbol, int timeframe = 0)
{
    double tmpbuffer[];
    int CurrentPeriod = NULL;
    if (timeframe == 0)	timeframe = Period();
    int OverSamplingRatio = Period();
    OverSamplingRatio = OverSamplingRatio / timeframe;
    int AvailableBarsForOverSampling =ArrayCopyRates(array, symbol, timeframe);
    int RequiredBarsForOversampling = (end + 1) * OverSamplingRatio - 1;
    int BarsForOversampling = (end - start + 1) * OverSamplingRatio;

    if (AvailableBarsForOverSampling < RequiredBarsForOversampling)
    end = ( AvailableBarsForOverSampling + 1 ) / OverSamplingRatio - 1; 

    int start_oversampling = start * OverSamplingRatio;
    int end_oversampling = end * OverSamplingRatio;

    ArrayResize(tmpbuffer, BarsForOversampling);
    int j = 0;
//Warning! It is very easily to violate memory access.
//Take care for over sampling ratio.
    for (int i = start_oversampling; i <= end_oversampling; i++) {
	tmpbuffer[j] = (array[i][1] + array[i + 1][1] + array[i + 2][1] + array[i + 3][1]) / 4.0;
	j++;
    }
    for (i = start; i <= end; i++) {
	double tmp = 0.0;
	for (j = 0; j < OverSamplingRatio; j++) tmp = tmp + tmpbuffer[i * OverSamplingRatio + j];
	tmp = tmp / OverSamplingRatio;
	retbuffer[i] = tmp;
    }
}
