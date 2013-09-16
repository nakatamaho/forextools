MT4DIR=/c/Program\ Files/Alpari\ Japan\ DEMO\ 1
MT4INDICATORDIR=$(MT4DIR)/experts/indicators/
MT4LIBRARYDIR=$(MT4DIR)/experts/libraries/
INDICATOR=instantaneousfrequency
MINGWDLL=libgcc_s_dw2-1.dll
INTERFACE=fftwinterface

all: $(INDICATOR).ex4 $(INTERFACE).dll

$(INTERFACE).dll: $(INTERFACE).cpp
	g++ -I/usr/local/include -march=corei7 -O2 -shared -c "$(INTERFACE)".cpp
	g++ -shared -o "$(INTERFACE)".dll "$(INTERFACE)".o -Wl,--kill-at,--output-def,"$(INTERFACE)".def,--out-implib,"$(INTERFACE).a" -L/usr/local/lib -Wl,-Bstatic -lfftw3
	objdump.exe -p "$(INTERFACE)".dll | grep "\.dll"

$(INDICATOR).ex4: $(INDICATOR).mq4
	$(MT4DIR)/MetaLang.exe $(INDICATOR).mq4

install:
	rm -f $(MT4INDICATORDIR)/$(INDICATOR).*
	rm -f $(MT4LIBRARYDIR)/$(INTERFACE).*
	rm -f $(MT4DIR)/$(MINGWDLL) 

	cp $(INDICATOR).ex4 $(MT4INDICATORDIR)
	cp $(INTERFACE).dll $(MT4LIBRARYDIR)
	cp /mingw/bin/$(MINGWDLL) $(MT4DIR)/

clean:
	rm -f *~ *def *a *.o *ex4 *.log $(INTERFACE).dll
