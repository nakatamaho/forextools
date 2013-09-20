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
#property link "https://github.com/nakatamaho/forextools/metatrader4library"

// This expert advisor obtains 
//1) Offset of UTC time from server time
//    int GetUTCOffset();
//2) GetCurrent Time in UTC 
//    datetime GetCurrentTimeInUTC();

#define ExpertAdvisorName "Time Utility Expert Advisor 1.0"

#import "Kernel32.dll"
void GetSystemTime(int &WinAPI_SYSTEMTIME[]);
void GetLocalTime(int &WinAPI_SYSTEMTIME[]);
#import

int init()
{
    return (0);
}

int deinit()
{
    return (0);
}

void start()
{
    while (!IsStopped() && IsExpertEnabled()) {
	Comment("-------------------------------------------------\n",
		"Server Time    : ", TimeToStr(TimeCurrent()), "\n",
		"Local Time     : ", TimeToStr(TimeLocal()), "\n",
		"UTC Time       : ", TimeToStr(GetCurrentTimeInUTC()), "\n",
		"UTC(GMT) Offset: ", GetUTCOffset(), "\n", 
		"-------------------------------------------------\n");
	WindowRedraw();
	Sleep(1000);
    }
}

datetime GetCurrentTimeInUTC()
{
    int WinAPI_SYSTEMTIME[4];
    GetSystemTime(WinAPI_SYSTEMTIME);	//Obtain current time in UTC using WinAPI
    int wYear = WinAPI_SYSTEMTIME[0] & 65535;
    int wMonth = WinAPI_SYSTEMTIME[0] >> 16;
    int wDay = WinAPI_SYSTEMTIME[1] >> 16;
    int wHour = WinAPI_SYSTEMTIME[2] & 65535;
    int wMinute = WinAPI_SYSTEMTIME[2] >> 16;
    int wSecond = WinAPI_SYSTEMTIME[3] & 65535;
//Dirty casts. May not work for newer systems.
    string sYear = wYear;
    string sMonth = wMonth;
    string sDay = wDay;
    string sHour = wHour;
    string sMinute = wMinute;
    string sSecond = wSecond;
    datetime CurrentTimeInUTC = StrToTime(StringConcatenate(sYear, ".", sMonth, ".", sDay, " ", sHour, ":", sMinute, ":", sSecond));
    return (CurrentTimeInUTC);
}

int GetUTCOffset()
{
    datetime CurrentTimeInUTC = GetCurrentTimeInUTC();
    int UTCOffset = (MathRound(TimeCurrent() - CurrentTimeInUTC + 600.0) / 3600.0);
    return (UTCOffset);
}
