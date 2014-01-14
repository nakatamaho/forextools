#!/bin/sh
uname -a
OPENBLASVERSION=0.2.8
rm -rf OpenBLAS-v$OPENBLASVERSION.tar.gz OpenBLAS-$OPENBLASVERSION
wget --no-check-certificate http://github.com/xianyi/OpenBLAS/archive/v$OPENBLASVERSION.tar.gz -O OpenBLAS-v$OPENBLASVERSION.tar.gz
tar xvfz OpenBLAS-v$OPENBLASVERSION.tar.gz
cd OpenBLAS-0.2.8
make BINARY=32 CC=gcc FC=gfortran SMP_SERVER=yes NUM_THREADS=8 DYNAMIC_ARCH=1
make PREFIX=/usr/local install BINARY=32 CC=gcc FC=gfortran SMP_SERVER=yes NUM_THREADS=8 DYNAMIC_ARCH=1
