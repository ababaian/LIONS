#!/bin/bash
set -e
set -o pipefail
set -x

db=$1

echo "Generating RepeatMasker data from db $db"

# chr1    3000001 3000156 -       L1_Mur2 LINE    L1

if [ "$db" == "hg19" ]; then
	mysql -h ucsc -u ucsc -pucsc -e "use $db -A; SELECT genoName,genoStart,genoEnd,strand,repName,repClass,repFamily FROM rmsk" | tail -n +2 > RepeatMaskerRaw_$db
else
	rm -f RepeatMaskerRaw_$db
	for chr in `seq 1 22` "X" "Y"
	do
		table=chr$chr"_rmsk"
		mysql -h ucsc -u ucsc -pucsc -e "use $db -A; SELECT genoName,genoStart,genoEnd,strand,repName,repClass,repFamily FROM $table" | tail -n +2 >> RepeatMaskerRaw_$db
	done
	#exit
fi

awk '{ if ( $6=="SINE" || $6=="LINE" || $6=="LTR" || $6=="SINE?" || $6=="LINE?" || $6=="LTR?") print $0 }' RepeatMaskerRaw_$db > SINES_LINES_LTRS_$db

awk '{ if ( $6=="DNA" || $6=="DNA?" || $6=="SINE" || $6=="LINE" || $6=="LTR" || $6=="SINE?" || $6=="LINE?" || $6=="LTR?" || $6=="Unknown" || $6=="Unknown?" || $6=="Other") print $0 }' RepeatMaskerRaw_$db > ForChimericSearch_$db

awk '{ print $1"\t"$2"\t"$3"\t"$5":"$6":"$7 }' RepeatMaskerRaw_$db > ForChimericSearch_$db.bed