#!/bin/bash
LAPACKVERSION=3.5.0
rm -rf lapack-$LAPACKVERSION #lapack-$LAPACKVERSION.tgz 
#wget http://www.netlib.org/lapack/lapack-$LAPACKVERSION.tgz
tar xvfz lapack-$LAPACKVERSION.tgz
cd lapack-$LAPACKVERSION
patch -p1 < ../patch.for_lapack-$LAPACKVERSION
cp INSTALL/make.inc.gfortran make.inc

echo Building reference BLAS
cd BLAS/SRC ; make ; cd ../../
echo Building reference LAPACK
make
echo Building LAPACKe
cd lapacke
make
echo All done
