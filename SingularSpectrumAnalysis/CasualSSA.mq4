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
#property indicator_buffers 5
#property indicator_color1 Black
#property indicator_color2 Blue
#property indicator_color3 Red
#property indicator_color4 LimeGreen
#property indicator_color5 DarkOrange
#property indicator_style1 STYLE_SOLID
#property indicator_style2 STYLE_SOLID
#property indicator_style3 STYLE_SOLID
#property indicator_width1 1
#property indicator_width2 4
#property indicator_width3 4
#property indicator_width4 2
#property indicator_width5 2

#import "libSSA4FX.dll"
void BasicSSA(double &x[], int N, int L, int Rmax, double &xtilde[]);
void BasicSSA_LAPACK(double &x[], int N, int L, int Rmax, double &xtilde[]);
#import

enum IS_USE_MA{
  USE_MA,
  DONOT_USE_MA,
};

enum SVD_METHOD{
  PROPACK,
  LAPACK,
};

enum PUSH_NOTIFICATION{
  USE_PUSH_NOTIFICATION,
  DONOT_USE_PUSH_NOTIFICATION,
};

enum ADAPTATION{
  USE_ADAPTATION,
  DONOT_USE_ADAPTATION,
};

input int TotalLength = 512;
input int WindowLength = 60;
input int CheckLength = 12;

input int Rmax = 2;
input IS_USE_MA use_ma = USE_MA;
input SVD_METHOD svd_method = PROPACK;
input PUSH_NOTIFICATION push_notification = USE_PUSH_NOTIFICATION;
input ADAPTATION adaptation = DONOT_USE_ADAPTATION;
input int MA_Period = 2;
input int STEALTH_MODE = 0;

double SSABuffer[];
double ExtBuffer[];
double PriceBuffer[];
double UpLine[];
double DnLine[];
double UpArrow[];
double DnArrow[];
int barcounter4calculation;
int _WindowLength;

int OnInit(void)
{
    barcounter4calculation = 0;
    _WindowLength = WindowLength;
    IndicatorBuffers(5);
    IndicatorDigits(Digits);

    SetIndexStyle(0, DRAW_LINE);
    SetIndexBuffer(0, ExtBuffer);
    SetIndexShift(0, 0);
    SetIndexLabel(0, "Basic SSA");
    SetIndexDrawBegin(0, TotalLength);
    ArraySetAsSeries(ExtBuffer, true);

    SetIndexStyle(1, DRAW_LINE);
    SetIndexBuffer(1, UpLine);
    SetIndexShift(1, 0);
    SetIndexDrawBegin(1, TotalLength);
    ArraySetAsSeries(UpLine, true);

    SetIndexStyle(2, DRAW_LINE);
    SetIndexBuffer(2, DnLine);
    SetIndexShift(2, 0);
    SetIndexDrawBegin(2, TotalLength);
    ArraySetAsSeries(DnLine, true);

    SetIndexStyle(3, DRAW_ARROW);
    SetIndexBuffer(3, UpArrow);
    SetIndexShift(3, 0);
    SetIndexDrawBegin(3, TotalLength);
    ArraySetAsSeries(UpArrow, true);
    SetIndexArrow(3, 233);

    SetIndexStyle(4, DRAW_ARROW);
    SetIndexBuffer(4, DnArrow);
    SetIndexShift(4, 0);
    SetIndexDrawBegin(4, TotalLength);
    ArraySetAsSeries(DnArrow, true);
    SetIndexArrow(4, 234);

    return (INIT_SUCCEEDED);
}

void OnDeinit(const int reason)
{
    //IndicatorRelease(MA_handle);
}

int OnCalculate(const int rates_total, const int prev_calculated, const datetime & time[], const double &open[], const double &high[], const double &low[], const double &close[], const long &tick_volume[], const long &volume[], const int &spread[])
{
    int i, j;
    int limit = rates_total - prev_calculated;
    string notification_message;
    //Calculation shold be done on the first tick of the new bar.
    if (rates_total == barcounter4calculation) return (rates_total);

    if (rates_total <= TotalLength * 2) {
	printf("Error # of rates_total is too small (%d).", rates_total);
	return (0);
    }
    if (rates_total - TotalLength < limit) limit = rates_total - TotalLength;
    if (limit == 0) limit = 1;

    ArraySetAsSeries(SSABuffer, true);
    ArraySetAsSeries(PriceBuffer, true);
    ArraySetAsSeries(close, true);

    ArrayResize(SSABuffer, TotalLength);
    ArrayResize(PriceBuffer, TotalLength);

    for (i = limit - 1; i >= 0; i--) {
	for (j = TotalLength - 1; j >= 0; j--) {
	    if (use_ma == USE_MA) {
		PriceBuffer[j] = iMA(NULL, 0, MA_Period, 0, MODE_EMA, PRICE_CLOSE, i + j);
	    } else {
		PriceBuffer[j] = close[j + i];
	    }
	}
	if (svd_method == PROPACK) {
	    BasicSSA(PriceBuffer, TotalLength, _WindowLength, Rmax, SSABuffer);
        } else {
	    BasicSSA_LAPACK(PriceBuffer, TotalLength, _WindowLength, Rmax, SSABuffer);
        }
	ExtBuffer[i] = SSABuffer[0];
   
    }
/*
        if (adaptation == USE_ADAPTATION) {
        double RMSE = 0.0;
        for (i = CheckLength - 1; i >= 0; i--) {
            RMSE = RMSE + (close [i] - ExtBuffer[i]) * (close [i] - ExtBuffer[i]);
        }
        RMSE = MathSqrt (RMSE / (double)CheckLength) * 1000.0;
        if (RMSE > 0.8)  if (_WindowLength > 30) _WindowLength = _WindowLength - 5;
        if (RMSE <= 0.8) if (_WindowLength < 120) _WindowLength = _WindowLength + 5;
        printf ("WindowLength : %g RMSE %g", _WindowLength, RMSE);
    }
*/
    for (i = limit - 1; i >= 0; i--) {
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
    for (i = limit - 1; i >= 0; i--) {
	UpArrow[i] = EMPTY_VALUE;
	DnArrow[i] = EMPTY_VALUE;
	if (UpLine[i + 1] != EMPTY_VALUE && UpLine[i] == EMPTY_VALUE)
	    DnArrow[i] = ExtBuffer[i];
	if (DnLine[i + 1] != EMPTY_VALUE && DnLine[i] == EMPTY_VALUE)
	    UpArrow[i] = ExtBuffer[i];
    }
    if (UpLine[1] != EMPTY_VALUE && UpLine[2] == EMPTY_VALUE) {
	if (rates_total > barcounter4calculation) {
	    if (STEALTH_MODE == 0) {
                double _UpLine = iCustom (NULL, PERIOD_H1, "CasualSSA", 256, 30, CheckLength, Rmax, USE_MA, PROPACK, DONOT_USE_PUSH_NOTIFICATION, DONOT_USE_ADAPTATION, MA_Period, 1, 1, 1);
	        if (_UpLine != EMPTY_VALUE) 
	           notification_message = StringFormat("CasualSSA: %s, %d Timeframe, STRONG BUY @%g %s",Symbol(), Period(), Bid, TimeToStr(TimeGMT()));
                else
  	           notification_message = StringFormat("CasualSSA: %s, %d Timeframe, BUY @%g %s",Symbol(), Period(), Bid, TimeToStr(TimeGMT()));
	        Print(notification_message);
	        if (push_notification == USE_PUSH_NOTIFICATION) SendNotification(notification_message);
            }
	    barcounter4calculation = rates_total;
	}
    }
    if (DnLine[1] != EMPTY_VALUE && DnLine[2] == EMPTY_VALUE) {
	if (rates_total > barcounter4calculation) {
	    if (STEALTH_MODE == 0) {
                double _DnLine = iCustom (NULL, PERIOD_H1, "CasualSSA", 256, 30, CheckLength, Rmax, USE_MA, PROPACK, DONOT_USE_PUSH_NOTIFICATION, DONOT_USE_ADAPTATION, MA_Period, 1, 2, 1);
	        if (_DnLine != EMPTY_VALUE) 
                    notification_message = StringFormat("CasualSSA: %s, %d Timeframe, STRONG SELL @%g %s", Symbol(), Period(), Ask, TimeToStr(TimeGMT()));
                else
                    notification_message = StringFormat("CasualSSA: %s, %d Timeframe, SELL @%g %s", Symbol(), Period(), Ask, TimeToStr(TimeGMT()));
	        Print(notification_message);
   	        if (push_notification == USE_PUSH_NOTIFICATION) SendNotification(notification_message);
            }
	    barcounter4calculation = rates_total;
	}
    }
    barcounter4calculation = rates_total;
    return (rates_total);
}
