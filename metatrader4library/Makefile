MT4DIR=c:/Program\ Files/Alpari\ DEMO\ 1
MSYSHOME=c:/mingw/msys/1.0
MINGWHOME=c:/mingw
MT4EXPERTDIR=$(MT4DIR)/experts
MT4INDICATORDIR=$(MT4DIR)/experts/indicators
MT4LIBRARYDIR=$(MT4DIR)/experts/libraries
EA=TimeUtility
LIBGCCDLL=libgcc_s_dw2-1.dll
GXX=mingw32-c++

all: $(EA).ex4

$(EA).ex4: $(EA).mq4
	$(MT4DIR)/MetaLang.exe $(EA).mq4

install:
	rm -f $(MT4EXPERTDIR)/$(EA).*
	rm -f $(MT4DIR)/$(LIBGCCDLL)

	cp $(EA).ex4 $(MT4EXPERTDIR)

clean:
	rm -f *~ *def *a *.o *ex4 *.log
