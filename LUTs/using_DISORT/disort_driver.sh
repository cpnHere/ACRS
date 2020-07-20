#!/bin/sh
# Chamara Rajapakshe (2020)
#----------------------------------------------
jobid=$1
inp=$2
outp=$3
cp ${inp}${jobid}_inputFile.dat inputFile.dat
rm INTENSITY.dat
./single_lay_cld.exe > single_lay_cld.out
cp -i INTENSITY.dat ${outp}${jobid}_INTENSITY.dat
cp -i single_lay_cld.out ${outp}${jobid}_single_lay_cld.out 
