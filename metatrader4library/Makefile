MT4DIR=c:/Program\ Files/Alpari\ DEMO\ 1
MSYSHOME=c:/mingw/msys/1.0
MINGWHOME=c:/mingw
MT4EXPERTDIR=$(MT4DIR)/experts
MT4INDICATORDIR=$(MT4DIR)/experts/indicators
MT4LIBRARYDIR=$(MT4DIR)/experts/libraries
EA1=TimeUtility
EA2=MiniScreen
EA3=EAtemplate
LIBGCCDLL=libgcc_s_dw2-1.dll
GXX=mingw32-c++

all: $(EA1).ex4 $(EA2).ex4 $(EA3).ex4

$(EA1).ex4: $(EA1).mq4
	$(MT4DIR)/MetaLang.exe $<

$(EA2).ex4: $(EA2).mq4
	$(MT4DIR)/MetaLang.exe $<

$(EA3).ex4: $(EA3).mq4
	$(MT4DIR)/MetaLang.exe $<

install:
	rm -f $(MT4EXPERTDIR)/$(EA1).*
	rm -f $(MT4EXPERTDIR)/$(EA2).*
	rm -f $(MT4EXPERTDIR)/$(EA3).*
	rm -f $(MT4DIR)/$(LIBGCCDLL)

	cp $(EA1).ex4 $(MT4EXPERTDIR)
	cp $(EA2).ex4 $(MT4EXPERTDIR)
	cp $(EA3).ex4 $(MT4EXPERTDIR)

clean:
	rm -f *~ *def *a *.o *ex4 *.log
