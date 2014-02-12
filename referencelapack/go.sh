#!/bin/bash
LAPACKVERSION=3.5.0
rm -rf lapack-$LAPACKVERSION #lapack-$LAPACKVERSION.tgz 
#wget http://www.netlib.org/lapack/lapack-$LAPACKVERSION.tgz
tar xvfz lapack-$LAPACKVERSION.tgz
cd lapack-$LAPACKVERSION
patch -p0 < ../patch.for_lapack-$LAPACKVERSION
cp INSTALL/make.inc.gfortran make.inc

echo Building reference BLAS
cd BLAS/SRC ; make ; cd ../../
echo Building reference LAPACK
make
echo Building LAPACKe
cd lapacke
make
echo All done
cd ..
cp lapack-$LAPACKVERSION/liblapack.a  c:/Mingw/msys/1.0/local/lib
cp lapack-$LAPACKVERSION/liblapacke.a c:/Mingw/msys/1.0/local/lib
cp lapack-$LAPACKVERSION/liblblas.a c:/Mingw/msys/1.0/local/lib/
cp lapack-$LAPACKVERSION/lapacke/include/lapacke.h c:/Mingw/msys/1.0/local/include/
cp lapack-$LAPACKVERSION/lapacke/include/lapacke_config.h c:/Mingw/msys/1.0/local/include/
cp lapack-$LAPACKVERSION/lapacke/include/lapacke_mangling.h c:/Mingw/msys/1.0/local/include/
cp lapack-$LAPACKVERSION/lapacke/include/lapacke_utils.h c:/Mingw/msys/1.0/local/include/
