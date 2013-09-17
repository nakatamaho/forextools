MT4DIR=c:/Program\ Files/Alpari\ DEMO\ 1
MSYSHOME=c:/mingw/msys/1.0
MINGWHOME=c:/mingw
MT4EXPERTDIR=$(MT4DIR)/experts
MT4INDICATORDIR=$(MT4DIR)/experts/indicators
MT4LIBRARYDIR=$(MT4DIR)/experts/libraries
INDICATOR=dctdifferentiation
LIBGCCDLL=libgcc_s_dw2-1.dll
FFTWINTERFACE=fftwinterface
GXX=mingw32-c++

all: $(INDICATOR).ex4 $(FFTWINTERFACE).dll

$(FFTWINTERFACE).dll: ../$(FFTWINTERFACE)/$(FFTWINTERFACE).cpp
	$(GXX) -I$(MSYSHOME)/local/include -march=corei7 -O2 -shared -c ../$(FFTWINTERFACE)/$(FFTWINTERFACE).cpp
	$(GXX) -shared -o $(FFTWINTERFACE).dll $(FFTWINTERFACE).o -Wl,--kill-at,--output-def,$(FFTWINTERFACE).def,--out-implib,$(FFTWINTERFACE).a -L$(MSYSHOME)/local/lib -Wl,-Bstatic -lfftw3
	#objdump.exe -p $(FFTWINTERFACE).dll | grep "\.dll"

$(INDICATOR).ex4: $(INDICATOR).mq4
	$(MT4DIR)/MetaLang.exe $(INDICATOR).mq4

install:
	rm -f $(MT4INDICATORDIR)/$(INDICATOR).*
	rm -f $(MT4LIBRARYDIR)/$(FFTWINTERFACE).*
	rm -f $(MT4DIR)/$(LIBGCCDLL)

	cp $(INDICATOR).ex4 $(MT4INDICATORDIR)
	cp $(FFTWINTERFACE).dll $(MT4LIBRARYDIR)
	cp $(MINGWHOME)/bin/$(LIBGCCDLL) $(MT4DIR)

clean:
	rm -f *~ *def *a *.o *ex4 *.log $(FFTWINTERFACE).dll $(FFTWINTERFACE).o
