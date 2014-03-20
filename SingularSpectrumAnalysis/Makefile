MT4DIR0 = c:/Program\ Files/Alpari\ UK\ DEMO\ 1
MT4DIR1 = c:/Program\ Files/Alpari\ JP\ DEMO\ 1
MT4DIR2 = c:/Program\ Files/Pepperstone\ DEMO\ 1 
MT4DIR3 = c:/Program\ Files/Thinkforex\ DEMO\ 1
MT4DIR4 = c:/Program\ Files/IC\ Markets\ DEMO\ 1
MT4DIR5 = c:/Program\ Files/Citi\ FX\ DEMO\ 1
MT4DIR6 = c:/Program\ Files/FXDD\ DEMO\ 1
MT4DIR7 = c:/Program\ Files/FOREX.com\ UK\ DEMO\ 1
MT4DIR9 = c:/Program\ Files/InstaForex\ DEMO\ 1
MT4DIR10 = c:/Program\ Files/MRC\ Market\ DEMO\ 1
MT4DIR11 = c:/Program\ Files/OANDA\ DEMO\ 1
MT4DIR12 = c:/Program\ Files/JFD\ DEMO\ 1
MT4DIRs = $(MT4DIR0) $(MT4DIR1) $(MT4DIR2) $(MT4DIR3) $(MT4DIR4) $(MT4DIR5) $(MT4DIR6) $(MT4DIR7) $(MT4DIR8) $(MT4DIR9) $(MT4DIR10) $(MT4DIR11) $(MT4DIR12)

MQL = $(MT4DIR0)/mql.exe

MSYSHOME        = c:/mingw/msys/1.0
MINGWHOME       = c:/mingw
MT4EXPERTDIR    = MQL4/Experts
MT4INDICATORDIR = MQL4/Indicators
MT4LIBRARYDIR   = MQL4/Libraries

INDICATORs_SRC  = BasicSSA.mq4 CasualSSA.mq4
INDICATORs      = $(INDICATORs_SRC:%.mq4=%.ex4)
LIBSSA		= libSSA4FX
LAPACKLIB       = $(MSYSHOME)/local/lib/libopenblas.a
#LAPACKLIB       = $(MSYSHOME)/local/lib/liblapacke.a $(MSYSHOME)/local/lib/liblapack.a $(MSYSHOME)/local/lib/libblas.a -lgfortran
FFTWLIB         = -L$(MSYSHOME)/local/lib -Wl,-Bstatic -lfftw3 -lm
STATICALIZER    = -static -lgfortran -static-libgfortran -static-libgcc

.SUFFIXES: .mq4 .ex4
.mq4.ex4:
	unix2dos < $< > $<_
	$(MQL) //mql4 $<_ 
	rm $<_

all: $(INDICATORs) $(LIBSSA).dll

$(LIBSSA).dll: $(LIBSSA).cpp
	$(CXX) -I$(MSYSHOME)/local/include -O2 -shared -c $(LIBSSA).cpp
	$(CXX) -shared -o $(LIBSSA).dll $(LIBSSA).o $(LAPACKLIB) $(FFTWLIB) $(STATICALIZER) -Wl,--kill-at,--output-def,$(LIBSSA).def,--out-implib,$(LIBSSA).a
	objdump.exe -p $(LIBSSA).dll | grep "\.dll"

install:
	for mt4dir in $(MT4DIRs) ; do \
	( \
		TMPDIR=`echo $$mt4dir | sed 's/ /\\ /g' ` ;  \
		cp $(LIBSSA).dll "$${TMPDIR}/$(MT4LIBRARYDIR)" ; \
		for indicator in $(INDICATORs) ; do \
		( \
			rm -f "$${TMPDIR}/$(MT4INDICATORDIR)/$$indicator" ; \
			cp $$indicator "$${TMPDIR}/$(MT4INDICATORDIR)/" ; \
		) ;\
		done \
	) ;\
	done
clean:
	rm -f *~ *.o *.log \#* *_ $(INDICATORs) ll mql_lapack.mqh $(LIBSSA).dll $(LIBSSA).a $(LIBSSA).def
