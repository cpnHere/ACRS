#!/bin/sh
# Chamara Rajapakshe (2020)
#----------------------------------------------
jobid=$1
cp ${jobid}_inputFile.dat inputFile.dat
rm INTENSITY.dat
time ./single_lay_cld.exe > single_lay_cld.out
cp -i INTENSITY.dat results/${jobid}_INTENSITY.dat
cp -i single_lay_cld.out results/${jobid}_single_lay_cld.out 
