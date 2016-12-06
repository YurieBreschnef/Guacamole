#!/bin/bash -e
# usage: call bash script with general plotting arguments

# gnuplot XDIM YDIM LX LY ASPRATIO SINDEX EINDEX
#XDIM     = $1
#YDIM     = $2
#LX       = $3
#LY       = $4
#ASPRATIO = $5
#SINDEX   = $6  # staring file index  for datafile numbering SINDEX.type.dat
#EINDEX   = $7  # end index
#MAXSPEC  = $8
#MINSPEC  = $9

echo "bash: ---plot all---"
  echo "  -bash: plot combo"
  gnuplot -e "XDIM=$1;YDIM=$2;LX=$3;LY=$4;ASPRATIO=$5;SINDEX=$6;EINDEX=$7;MAXSPEC=$8;MINSPEC=$9;" ./output/plot_combo.gp
  echo "  -bash: plot spec_inspect"
  gnuplot -e "XDIM=$1;YDIM=$2;LX=$3;LY=$4;ASPRATIO=$5;SINDEX=$6;EINDEX=$7;MAXSPEC=$8;MINSPEC=$9;" ./output/plot_spec.gp
  echo "  -bash: plot_stat"
  gnuplot -e "XDIM=$1;YDIM=$2;LX=$3;LY=$4;ASPRATIO=$5;SINDEX=$6;EINDEX=$7;MAXSPEC=$8;MINSPEC=$9;" ./output/plot_stat.gp



