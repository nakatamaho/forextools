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

//Using this expert advisor, you can open a miniscreen on your chart.
//This screen can be used casual "print" debug.

#define ExpertAdvisorName "MiniScreen Expert Advisor"

#define MINISCREENMAX_WIDTH_X 60
#define MINISCREENMAX_WIDTH_Y 40
#define LINEPROCESS_DONOTHING 0
#define LINEPROCESS_SCROLL_UP 1
#define LINEPROCESS_CRLF 2
string MiniScreenFrameBuffer[MINISCREENMAX_WIDTH_X][MINISCREENMAX_WIDTH_Y];
int TotalObjectStrings = 0;
int MiniScreenCursorX = 0;
int MiniScreenCursorY = 0;
int MiniScreenOffsetX = 0;
int MiniScreenOffsetY = 0;
int MiniScreenWidthX;
int MiniScreenWidthY;
string FontName = "Helvetica";
int FontSize = 12;
bool MiniScreenInitialized = False;

int init()
{
    InitMiniScreen(30,10);
    UpdateScreen();
    return (0);
}

int deinit()
{
    DeInitMiniScreen();
    return (0);
}


void start()
{
    CommentEx("TEST TEST TEST TEST");
}

void InitMiniScreen(int widthx, int widthy)
{
    MiniScreenInitialized = False;
    CleanMiniScreen();
    if (widthx >= 60 && widthx < 5) {
	Print("Invalid width X");
	return;
    }
    if (widthy >= 40 && widthy < 5) {
	Print("Invalid width Y");
	return;
    }
    MiniScreenWidthX = widthx;
    MiniScreenWidthY = widthy;
    MiniScreenInitialized = True;
}

void DeInitMiniScreen()
{
    CleanMiniScreen();
    MiniScreenInitialized = False;
}

int LineProcessingEvaluation()
{
    if (MiniScreenCursorX == MiniScreenWidthX) {
	if (MiniScreenCursorY == MiniScreenWidthY - 1) {
	    return (LINEPROCESS_SCROLL_UP);
	} else {
	    return (LINEPROCESS_CRLF);
	}
    }
    return (LINEPROCESS_DONOTHING);
}

void FrameBufferPutChar(string character)
{
    string CRLFString = "\n";
    if (StringFind(character, CRLFString) != -1) {
	if (MiniScreenCursorY == MiniScreenWidthY - 1) {
	    MiniScreenScrollUp();
	    MiniScreenCursorX = 0;
	} else {
	    MiniScreenCursorX = 0;
	    MiniScreenCursorY++;
	}
	return;
    }
    MiniScreenFrameBuffer[MiniScreenCursorX][MiniScreenCursorY] = character;
    MiniScreenCursorX++;
    return;
}

void CommentEx(string CommentString = "")
{
    int len = StringLen(CommentString);
    int currentposition = 0;
    int cursorstatus;

    while (currentposition < len) {
	cursorstatus = LineProcessingEvaluation();
	switch (cursorstatus) {
	case LINEPROCESS_SCROLL_UP:
	    MiniScreenScrollUp();
	    MiniScreenCursorX = 0;
	    FrameBufferPutChar(StringSubstr(CommentString, currentposition, 1));
	    break;
	case LINEPROCESS_CRLF:
	    MiniScreenCursorY++;
	    MiniScreenCursorX = 0;
	    FrameBufferPutChar(StringSubstr(CommentString, currentposition, 1));
	    break;
	default:
	    FrameBufferPutChar(StringSubstr(CommentString, currentposition, 1));
	    break;
	}
	currentposition++;
    }
}

void MiniScreenScrollUp()
{
    int x, y;
    string NullString = "\x00";
    for (y = 1; y < MiniScreenWidthY; y++) {
	for (x = 0; x < MiniScreenWidthX; x++) {
	    MiniScreenFrameBuffer[x][y - 1] = MiniScreenFrameBuffer[x][y];
	}
    }
    for (x = 0; x < MiniScreenWidthX; x++) {
	    MiniScreenFrameBuffer[x][MiniScreenWidthY-1] = NullString;
    }
    UpdateScreen();
/*
    for (j = MiniScreenWidthY; j > 0; j--) {
	ClearOneLine(j - 1);
	WriteOneLineToWin(j - 1);
    }
*/
    return;
}

void CleanMiniScreen()
{
    MiniScreenCursorX = 0;
    MiniScreenCursorY = 0;
    for (int x = 0; x < MINISCREENMAX_WIDTH_X; x++) {
	for (int y = 0; y < MINISCREENMAX_WIDTH_Y; y++) {
	    MiniScreenFrameBuffer[x][y] = "\x00";
	}
    }
}

void UpdateScreen()
{
    string StringForScreen = "";
    string NullString = "\x00";
    string CRLFString = "\n";
    for (int y = 0; y < MiniScreenWidthY; y++) {
	for (int x = 0; x < MiniScreenWidthX; x++) {
		StringForScreen = StringConcatenate(StringForScreen, MiniScreenFrameBuffer[x][y]);
	}
	StringForScreen = StringConcatenate(StringForScreen, MiniScreenFrameBuffer[x][y], CRLFString);
    }
    Comment(StringForScreen);
/* For debugging
    for (y = 0; y < MiniScreenWidthY; y++) {
	for (x = 0; x < MiniScreenWidthX; x++) {
	    Print("MiniScreenFrameBuffer[", x, "][", y, "] ", "\"", MiniScreenFrameBuffer[x][y], "\"", " ", StringGetChar(MiniScreenFrameBuffer[x][y],0));
	}
    }
*/
}
