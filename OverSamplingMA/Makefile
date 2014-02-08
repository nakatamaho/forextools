MT4DIR0=c:/Program\ Files/Alpari\ DEMO\ 1
MT4DIR1=c:/Program\ Files/Alpari\ DEMO\ 2
MT4DIR2=c:/Program\ Files/Alpari\ DEMO\ 3
MT4DIR3=C:/Program\ Files/FOREX.com\ UK\ DEMO\ 1
MT4DIR4=C:/Program\ Files/IC\ Markets\ DEMO\ 1
MT4DIR5=C:/Program\ Files/InstaForex\ DEMO\ 1
77;10103;0cMT4DIR6=C:/Program\ Files/MRC\ DEMO\ 1
MT4DIR7=C:/Program\ Files/Pepperstone\ DEMO\ 1 
MT4DIR8=C:/Program\ Files/Saxo\ Banque\ DEMO\ 1
MT4DIR9=C:/Program\ Files/Thinkforex\ DEMO\ 1
MT4DIRs=$(MT4DIR0) $(MT4DIR1) $(MT4DIR2) $(MT4DIR3) $(MT4DIR4) $(MT4DIR5) $(MT4DIR6) $(MT4DIR7) $(MT4DIR8) $(MT4DIR9)
MSYSHOME=c:/mingw/msys/1.0
MINGWHOME=c:/mingw
MT4EXPERTDIR=experts
MT4INDICATORDIR=experts/indicators
MT4LIBRARYDIR=experts/libraries
GXX=mingw32-c++

INDICATORs_SRC=OverSamplingMA.mq4

INDICATORs=$(INDICATORs_SRC:%.mq4=%.ex4)
EAs=$(EAs_SRC:%.mq4=%.ex4)

.SUFFIXES: .mq4 .ex4
.mq4.ex4:
	$(MT4DIR0)/MetaLang.exe $<

all: $(INDICATORs) #$(EAs)

install:
	for mt4dir in $(MT4DIRs) ; do \
	( \
		for indicator in $(INDICATORs) ; do \
		( \
			TMPDIR=`echo $$mt4dir | sed 's/ /\\ /g' ` ;  \
			rm -f "$${TMPDIR}/$(MT4INDICATORDIR)/$$indicator" ; \
			cp $$indicator "$${TMPDIR}/$(MT4INDICATORDIR)/" ; \
		) ;\
		done \
	) ;\
	done
clean:
	rm -f *~ *.o *ex4 $(INDICATORs) $(EAs)


