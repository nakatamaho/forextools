#!/bin/bash
LAPACKVERSION=3.5.0
rm -rf lapack-$LAPACKVERSION #lapack-$LAPACKVERSION.tgz 
wget http://www.netlib.org/lapack/lapack-$LAPACKVERSION.tgz
tar xvfz lapack-$LAPACKVERSION.tgz
cd lapack-$LAPACKVERSION
patch -p0 < ../patch.for_lapack-$LAPACKVERSION
cp make.inc.example make.inc

echo Building reference BLAS
cd BLAS/SRC ; make ; cd ../../
echo Building reference LAPACK
make
echo Building LAPACKe
cd lapacke
make
echo All done
cd ..
cp liblapack.a  c:/Mingw/msys/1.0/local/lib
cp liblapacke.a c:/Mingw/msys/1.0/local/lib
cp libblas.a c:/Mingw/msys/1.0/local/lib/
cp lapacke/include/lapacke.h c:/Mingw/msys/1.0/local/include/
cp lapacke/include/lapacke_config.h c:/Mingw/msys/1.0/local/include/
cp lapacke/include/lapacke_mangling.h c:/Mingw/msys/1.0/local/include/
cp lapacke/include/lapacke_utils.h c:/Mingw/msys/1.0/local/include/
