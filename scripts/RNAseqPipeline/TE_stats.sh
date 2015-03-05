#!/bin/bash
set -e
set -o pipefail

SHELL_BASE="/home/ababaian/software/RNAseqPipeline"
#SHELL_BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

J=/gsc/software/linux-x86_64/jre-1.6.0_16/bin/java
#J=/gsc/software/linux-x86_64/jre1.7.0_03/bin/java

#RegionsCoverageFromWigCalculator=/projects/03/genereg/projects/SOLEXA/lib/RegionsCoverageFromWigCalculator_May6_2011.jar


RegionsCoverageFromWigCalculator=$SHELL_BASE/bin/RegionsCoverageFromWigCalculator_May6_2011.jar


# INPUT Parameters
	name=$1
	rpkm=$2
	species=$3
	sr=$4

	if [ $sr == 'R' ]; then
		wig=$5
	else
		wig_pos=$5
		wig_neg=$6
	fi

# Resources
	reschrs=$(sh $SHELL_BASE/RNAgetRes.sh $species)
	res=$(echo $reschrs | cut -f 1 -d ',')
	chrs=$(echo $reschrs | cut -f 2 -d ',')
	TEdata=$(echo $reschrs | cut -f 4 -d ',')
	TEout=`basename $TEdata`

echo "RPKM=$rpkm"

if [ "$sr" == "R" ]; then

	# run TE RPKM code
	

	$J -jar -Xmx15G $RegionsCoverageFromWigCalculator -w $wig -r $TEdata -o ./$name/ -s $chrs -n $name

	# From RPKM log:
	# RPKM coefficent:   (1000000000/100)/13540148=0.7385

	#echo "Command1"
	#awk '{print $1":"$2":"$3"\t"($3-$2)"\t"$5":"$4":"$6"\t"$7"\t"$8"\t"$9}' $name/$TEout.$name.coverage | sort -k3 | awk '{print $3"\t"$5"\t"$2}' | awk -v r=$rpkm '{if (I==$1 || count==0) {SumC=SumC+$2; count=count+1; SumD=SumD+$3;I=$1} else {v=r*SumC/SumD; print I"\t"count"\t"v; SumC=$2; SumD=$3; I=$1;count=1;}}' I="LINE:CR1:CR1_Mam" SumD=0 SumC=0 count=0 > $name/$TEout.$name.agglomorated.RPKM

	#echo "Command2"
	#awk '{print $1":"$2":"$3"\t"($3-$2)"\t"$5":"$4":"$6"\t"$7"\t"$10}' $name/$TEout.$name.coverage | sort -k3 | awk '{print $3"\t"$5"\t"$2}' | awk '{if (I==$1) {SumC=SumC+$2; count=count+1; SumD=SumD+$3;if ($3>maxVal) maxVal=$3;I=$1} else {print I"\t"count"\t"maxVal; SumC=$2; SumD=$3; I=$1;count=1; maxVal=0}}' I="LINE:CR1:CR1_Mam" SumD=0 SumC=0 count=0 maxVal=0 > $name/$TEout.$name.agglomorated.max
	
	#awk '{print $1":"$2":"$3"\t"($3-$2)"\t"$5":"$4":"$6"\t"$7"\t"$10}' $name/$TEout.$name.coverage | sort -k3 | awk '{print $3"\t"$5"\t"$2}' | awk -v r=$rpkm '{if (I==$1 || count==0) {SumC=SumC+$2; count=count+1; SumD=SumD+$3;if ($3>maxVal) maxVal=$3;I=$1} else {v=r*SumC/SumD; print I"\t"count"\t"maxVal"\t"v; SumC=$2; SumD=$3; I=$1;count=1; maxVal=0}}' I="LINE:CR1:CR1_Mam" SumD=0 SumC=0 count=0 maxVal=0 > $name/$TEout.$name.agglomorated
    
    #awk '{print $1":"$2":"$3"\t"($3-$2)"\t"$5":"$4":"$6"\t"$7"\t"$10}' RetroTransposons.Female.PGC.SETDB1cKO.RNA-seq.coverage | sort -k3 | awk '{print $3"\t"$5"\t"$2}' | awk -v r=$rpkm '{if (I==$1 || count==0) {SumC=SumC+$2; count=count+1; SumD=SumD+$3;if ($3>maxVal) maxVal=$3;I=$1} else {v=r*SumC/SumD; print I"\t"count"\t"maxVal"\t"v; SumC=$2; SumD=$3; I=$1;count=1; maxVal=0}}' I="LINE:CR1:CR1_Mam" SumD=0 SumC=0 count=0 maxVal=0
    
    awk '{print $1":"$2":"$3"\t"($3-$2)"\t"$5":"$4":"$6"\t"$7"\t"$8"\t"$9"\t"$10}' $name/$TEout.$name.coverage | awk '{print $3"\t"$5"\t"$2"\t"$7}' | sort -k1,1 | awk -v r=$rpkm '{if (I==$1 || count==0) {SumC=SumC+$2; count=count+1; SumD=SumD+$3;I=$1;if ($4>maxVal) maxVal=$4;} else {v=r*SumC/SumD; print I"\t"count"\t"v"\t"maxVal; SumC=$2; SumD=$3; I=$1; count=1; maxVal=0}}' I="" SumD=0 SumC=0 count=0 maxVal=0 > $name/$TEout.$name.agglomorated
	
else
	#echo "Strand specific TE stats not yet implemented"
	#exit 1

	$J -jar -Xmx15G $RegionsCoverageFromWigCalculator -w $wig_neg -r $TEdata -o ./$name/ -s $chrs -n $name.neg
	$J -jar -Xmx15G $RegionsCoverageFromWigCalculator -w $wig_pos -r $TEdata -o ./$name/ -s $chrs -n $name.pos

	gunzip $name/$TEout.$name.neg.coverage.gz
	gunzip $name/$TEout.$name.pos.coverage.gz
	awk '{print $1":"$2":"$3"\t"$5":"$4":"$6"\t"$7"\t"$8"\t"$9}' $name/$TEout.$name.pos.coverage | sort -k1 > x
	awk '{print $1":"$2":"$3"\t"$8"\t"$9}' $name/$TEout.$name.neg.coverage | sort -k1 > x1
	join x x1 | sed 's/ /\t/g' | sed 's/:/\t/' | sed 's/:/\t/' | awk '{print $1":"$2":"$3"\t"($3-$2)"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' > $name/$TEout.$name.pos.neg.coverage

	awk '{if($4=="-") {print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$8} else {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}' $name/$TEout.$name.pos.neg.coverage | sort -k3 | awk '{print $3"\t"$5"\t"$2}' | awk '{if (I==$1) {SumC=SumC+$2; count=count+1; SumD=SumD+$3;I=$1} else {print I"\t"SumC"\t"SumD"\t"count; SumC=$2; SumD=$3; I=$1;count=1;}}' I="LINE:CR1:CR1_Mam" SumD=0 SumC=0 count=0 > $name/$TEout.$name.agglomorated
	
	#awk '{if($4=="-") {print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$8} else {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}' $name/$TEout.$name.pos.neg.coverage | sort -k3 | awk '{print $3"\t"$5"\t"$2}' | awk -v r=$rpkm '{if (I==$1 || count==0) {SumC=SumC+$2; count=count+1; SumD=SumD+$3;if ($3>maxVal) maxVal=$3;I=$1} else {v=r*SumC/SumD; print I"\t"count"\t"maxVal"\t"v; SumC=$2; SumD=$3; I=$1;count=1; maxVal=0}}' I="LINE:CR1:CR1_Mam" SumD=0 SumC=0 count=0 maxVal=0 > $name/$TEout.$name.agglomorated

fi

echo "TE stats complete"
