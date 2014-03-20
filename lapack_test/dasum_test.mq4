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

#property copyright   "2014- NAKATA Maho"
#property link        "https://github.com/nakatamaho/forextools"

#property indicator_chart_window

#include "mql_lapack.mqh"

//Matlab/Octave format
void printvec(int N, double &x[])
{
    double mtmp;
    string STR;
    STR = "[ ";
    for (int i = 0; i < N; i++) {
        mtmp = x[i];
        STR = STR + StringFormat("%5.2e", mtmp);
        if (i < N-1) STR = STR + (", ");
    }
    STR = STR + ("]");
    Print(STR);
}

int init()
{
    int n = 4;
    double x[], dasum_ret;
    ArrayResize_8byte_aligned_double(x, n);

    x[0] = 1.0;    x[1] = 8.0;    x[2] = 3.0;    x[3] = 4.0;

    printf("x =");    printvec(n, x);    printf("\n");

    dasum_ret = mql_dasum(n, x, 1);
    printf("ans = %f", dasum_ret);  
    return (0);
}

int deinit()
{
    return (0);
}

int start()
{
    return (0);
}
