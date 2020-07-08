#!/bin/sh
# Chamara Rajapakshe (2020)
# Based on DISORT example
#----------------------------------------------
# This script compiles and executes DISORT
#
#DISORT path
dpath="/umbc/xfs1/zzbatmos/users/charaj1/taki/ACRS/DISORT/disort4.0.98"
rm -f a.out
cat READ_VAL.f90 $dpath/DISOTESTAUX.f $dpath/DISORT.f $dpath/BDREF.f $dpath/DISOBRDF.f $dpath/ERRPACK.f $dpath/LINPAK.f $dpath/LAPACK.f $dpath/RDI1MACH.f > code.f
#gfortran -O0 -g -fcheck=all -fdump-core -ffpe-trap=invalid,zero,overflow,underflow,denormal -Wall rte_driver.f90 code.f 
#gfortran -O0 -g -fcheck=all -fdump-core -ffpe-trap=invalid,zero,overflow -Wall rte_driver.f90 code.f 
#To run the example code
#gfortran -O3 -g -fcheck=all -fdump-core -fbounds-check -Wall disotest.f90 code.f -o disort4_unit_tests.exe
#chmod u+x ./disort4_unit_tests.exe
#time ./disort4_unit_tests.exe
#rm -f code.f
#Trying singel-layer cloud example
rm single_lay_cld.exe INTENSITY.dat
gfortran -O3 -g -fcheck=all -fdump-core -fbounds-check -Wall single_lay_cld.f90 code.f -o single_lay_cld.exe
chmod u+x ./single_lay_cld.exe
time ./single_lay_cld.exe > single_lay_cld.out
