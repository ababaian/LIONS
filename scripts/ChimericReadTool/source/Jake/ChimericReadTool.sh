#!/bin/bash
set -e
set -o pipefail

PWD=`pwd`

usage()
{
echo "Script to find chimeric reads, generate stats for reads and exons"
echo ""
echo "USAGE:"
echo "ChimericReadTool.sh <BAM file> <WIG coverage file> <Species>"
echo ""
}

if [ $# -ne 3 ]; then
	usage
	exit 1
fi

bam=$1
fwig=$2
species=$3
bai=$bam.bai

if [ ! -r $bam ]; then
	usage
	echo "ERROR: Unable to read BAM file"
	exit 1
fi

if [ ! -r $bai ]; then
	#usage
	echo "ERROR: Cannot find BAM index file ($bai)"
	echo "       Run samtools index on BAM file or find associated index file"
	exit 1
fi

if [ ! -r $fwig ]; then
	usage
	echo "ERROR: Unable to read WIG file"
	exit 1
fi

#input=hg18s_results
#chrs=/projects/03/genereg/projects/SOLEXA/chr_info/hg18.chrom.sizes	
#fwig=/projects/mbilenky/mlc/Jake/hs0988/wig/hs0988.q10.F516.wig.gz

SCRIPT_BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

python3=/gsc/software/linux-x86_64/python-3.2.2/bin/python3

reschrs=$(sh $SCRIPT_BASE/RNAgetRes.sh $species)
res=$(echo $reschrs | cut -f 1 -d ',')
repeats=$(echo $reschrs | cut -f 5 -d ',')
exons=$res/$species"_exons"
chrs=$(echo $reschrs | cut -f 2 -d ',')

# Example files
#exons=/projects/mbilenky/resources/mm9rs_0313/mm9rs_0313_exons
#repeats=/projects/mbilenky/resources/RepeatMasker/ForChimericSearch_mm9
#chrs=/projects/03/genereg/projects/SOLEXA/chr_info/mm9.chrom.sizes

#fwig=/projects/mbilenky/mlc/Jake/PGC_analysis/wig/rpkm_PGC_analysis.q10.F516.wig.gz
#bam=/projects/mbilenky/mlc/Jake/Lorincz_PGC/bams/RNA-seq.PGC.bam
input=tmp_chimeric_data

# Loads BAM file and finds chimeric reads and creates stats on exon/repeat interactions
# Defines ER/DR/etc reads and aggrogates results
#$python3 $SCRIPT_BASE/chimericReadSearch.py $exons $repeats $bam temp.bed > $input

# Put in the exon coords
less $input | cut -f 13,14,15 > $PWD/tmp_repeat_coords

# Put in the repeat coords
less $input | cut -f 13,16,17 >> $PWD/tmp_repeat_coords

# Put in upstream coords of the repeat
less $input | cut -f 13,16,17,18 | awk '{ 
if ($4==1) 
	print $1"\t"($2-50)"\t"($2); 
else 
	print $1"\t"($3)"\t"($3+50); 
}' >> $PWD/tmp_repeat_coords

sort $PWD/tmp_repeat_coords | uniq > $PWD/tmp_uniq_repeat_coords

#J=/gsc/software/linux-x86_64/jre1.7.0_03/bin/java
J=/gsc/software/linux-x86_64/jre-1.6.0_16/bin/java
JAVA_BASE=/projects/03/genereg/projects/SOLEXA/lib/

# Calculates max coverage of exons, repeat and upstream of repeats using Java tool
$J -jar -Xmx5G $JAVA_BASE/RegionsCoverageFromWigCalculator.jar -w $fwig -r $PWD/tmp_uniq_repeat_coords -s $chrs -o $PWD -n results

# All code below is taking coverage data and integrating back into "Excel" spreadsheet

less $input | awk '{ 
if ($18==1) 
	print $0"\tchr"$13":"$14"-"$15"\tchr"$13":"$16"-"$17"\tchr"$13":"($16-50)"-"($16); 
else 
	print $0"\tchr"$13":"$14"-"$15"\tchr"$13":"$16"-"$17"\tchr"$13":"($17)"-"($17+50); 
}' > $PWD/tmp_input_with_key

less tmp_uniq_repeat_coords.results.coverage | awk ' { print $1":"$2"-"$3"\t"$4"\t"$5 } ' | sort -k 1 > $PWD/tmp_sorted_repeat_data

# If adding new columns, the column indices used in the sort/join command need to increased
# at the moment: 21 is exon coordinates
#                22 is repeat coordinates
#                23 is 50bp upstream repeat coordinates

# Need to sort before join command
sort -k 21 $PWD/tmp_input_with_key > $PWD/tmp_input_sorted1
join -1 21 -2 1 $PWD/tmp_input_sorted1 $PWD/tmp_sorted_repeat_data | sed -e 's/ /\t/g' > $PWD/tmp_results_1

sort -k 22 $PWD/tmp_results_1 > $PWD/tmp_input_sorted2
join -1 22 -2 1 $PWD/tmp_input_sorted2 $PWD/tmp_sorted_repeat_data | sed -e 's/ /\t/g' > $PWD/tmp_results_2

sort -k 23 $PWD/tmp_results_2 > $PWD/tmp_input_sorted3
join -1 23 -2 1 $PWD/tmp_input_sorted3 $PWD/tmp_sorted_repeat_data | sed -e 's/ /\t/g' > $PWD/tmp_final

# Add header to final results

echo -en "transcriptID\texonRankInTranscript\trepeatName\tcoordinates\tER_Interaction\tIsExonic\tExonsOverlappingWithRepeat\t" > $PWD/final_results
echo -en "ER\tDR\tDE\tDD\tTotal\tChromosome\tEStart\tEEnd\tRStart\tREnd\tEStrand\tRStrand\tPrevExon\t" >> $PWD/final_results
echo -e "ExonRPKM\tExonMaxCoverage\tRepeatRPKM\tRepeatMaxCoverage\tUpstreamRepeatRPKM\tUpstreamRepeatMaxCoverage" >> $PWD/final_results

# Cut out the first three columns (just coordinate keys used for the sort/join commands above)
cut -f 4- $PWD/tmp_final  >> $PWD/final_results



