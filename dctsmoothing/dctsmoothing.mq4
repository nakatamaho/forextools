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
#property link      "https://github.com/nakatamaho"
#property indicator_chart_window
#property indicator_buffers 1

#import "fftwinterface.dll"
void fastcosinetransform(double& a[], int dctlength, int inversefct);
#import

extern double n = 12;
extern int cutoff_frequency = 1024;

double DCT_buffer[];
double tmp_buffer[];
int dctlength;
bool ErrorFlag = FALSE;

int init()
{
   Print ("initialized");
   SetIndexStyle(0, DRAW_LINE, EMPTY, 2, Red);
   SetIndexBuffer(0, DCT_buffer);

   dctlength = MathPow(2, n);

   int M = ArrayResize(tmp_buffer, dctlength);
   if (M < dctlength) { 
      Print("Cannot allocate memory", dctlength, M);
      ErrorFlag=TRUE;
      return (0);
   }
   SetIndexDrawBegin(0, Bars - M);
   return(0);
}

int deinit()
{
   return(0);
}

int start()
{
    int i;
    for(i = dctlength -1; i >= 0; i--) {
       tmp_buffer[i] = Close[i];
    }

//the discrete cosine (DCT-II) transformation    
    fastcosinetransform(tmp_buffer, dctlength, false);

//cutoff by rectangular window
    for(i = 0; i < dctlength; i++) {
        if(i >= cutoff_frequency) tmp_buffer[i] = .0;
    }

//reverse discrete cosine (DCT-III) transformation
    fastcosinetransform(tmp_buffer, dctlength, true);

    for(i = 0; i <= dctlength - 1; i++) {
        DCT_buffer[i] = tmp_buffer[i];
    }
    return(0);
}

