#!/bin/sh
set -e
set -o pipefail

usage()
{
echo "WIG modifying script using RPKM value"
echo ""
echo "USAGE (two input parameters):"
echo "/home/jlever/bin/WIG_RPKM.sh <wig_file> <rpkm> <wig_out>"
echo ""
echo "OUTPUT: wig file multipled by rpkm (as GZ archive)"
echo ""
}

wig=$1
rpkm=$2
out=$3

if [ $# -ne 3 ]; then
	usage
	exit 1
fi

if [ ! -f $wig ]; then
	usage
	echo "ERROR: Unable to access file: $wig"
	exit 1
fi

less $wig | awk -v rpkm=$rpkm '{if($1=="fixedStep") {print $0;} else {print rpkm*$1}}' | gzip -9 > $out

