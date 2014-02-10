UNAME      = $(shell uname |  sed 's/_.*//g')
USE_MKL    = yes
ifeq ($(UNAME),$(filter $(UNAME),CYGWIN MINGW32))
CXX        = g++
FC         = gfortran
CXXFLAGS   = -Ic:/mingw/msys/1.0/local/include -Wall
LDFLAGS    = -Lc:/mingw/msys/1.0/local/lib -static-libgcc -static-libgfortran -static-libstdc++
LIBS       = -lopenblas -lgfortran
endif
ifeq ($(UNAME),Darwin)
CXX        = g++-mp-4.8
FC         = gfortran-mp-4.8
CXXFLAGS   = -I/opt/local/include -DUSE_MACPORTS
LDFLAGS    = -L/opt/local/lib
LIBS       = -lopenblas -lgfortran
endif
ifeq ($(UNAME),Linux)
ifeq ($(USE_MKL),yes)
CXX        = icpc
FC         = ifort
CXXFLAGS   = -mkl -I/opt/intel/composerxe/mkl/include -DUSE_MKL
LDFLAGS    = -mkl
else
CXX        = g++
FC         = gfortran
CXXFLAGS   =
LDFLAGS    =
LIBS       = -lopenblas -lgfortran
endif
endif

TESTSRCS   = dgemm_test.cpp dgesvd_test.cpp dsyev_test.cpp
BENCHSRCS  = dgemm_bench.cpp dgesvd_bench.cpp dgesdd_bench.cpp dsyev_bench.cpp
TESTOBJS   = $(TESTSRCS:%.cpp=%.o)
TESTPROGS  = $(TESTOBJS:%.o=%)
BENCHOBJS  = $(BENCHSRCS:%.cpp=%.o)
BENCHPROGS = $(BENCHOBJS:%.o=%)
MTSRC	   = mt19937ar.c
MTOBJ	   = mt19937ar.o

all: $(TESTPROGS) $(BENCHPROGS)

$(MTOBJ): $(MTSRC)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(BENCHPROGS): %: %.o $(MTOBJ)
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS) $(MTOBJ)

$(TESTPROGS): %: %.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)

$(BENCHOBJS),$(TESTOBJS): %.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(TESTPROGS) $(TESTOBJS) $(BENCHPROGS) $(BENCHOBJS) $(MTOBJ) *~ \#* *.core
