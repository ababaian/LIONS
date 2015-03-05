#!/bin/bash
# Chimeric Read Tool v1.0 
# Mod by Artem
# Usage:
# ChimericReadTool.sh <input.bam> <Coverage.wig> <Reference_set>
# output goes to current working directory
# =====================================================================

set -e
set -o pipefail

PWD=`pwd`

usage()
	{
	echo "Script to find chimeric reads, generate stats for reads and exons"
	echo ""
	echo "USAGE:"
	echo "ChimericReadTool.sh <BAM file> <WIG coverage file> <Reference Set>"
	echo ""
	echo "     Human genome 19 with GENCODE v14 = hg19gc_v14"
	echo "     Note: BAM file should have associated bam.bai index file"
	echo ""
	}

if [ $# -ne 3 ]; then
	usage
	exit 1
fi

# INPUT Parameters ----------------------------------------------------
	bam=$1
	name=$(echo $1 | sed 's,^[^ ]*/,,g' | sed 's/\.bam//g' -)
	# For assembly-based method; takes the library name for accessing res
	fwig=$2
	species=$3
	bai=$bam.bai # assumes index is present

	#Check input is readable
	if [ ! -r $bam ]; then
		usage
		echo "ERROR: Unable to read BAM file"
		exit 1
	fi

	if [ ! -r $bai ]; then
		echo "ERROR: Cannot find BAM index file ($bai)"
		echo "       Run samtools index on BAM file or find associated index file"
		exit 1
	fi

	if [ ! -r $fwig ]; then
		usage
		echo "ERROR: Unable to read WIG file"
		exit 1
	fi

# INITIALIZATION ------------------------------------------------------

	SCRIPT_BASE='/home/ababaian/software/ChimericReadTool'
	#python3=/gsc/software/linux-x86_64/python-3.2.2/bin/python3	
	python3='/home/ababaian/software/python3/bin/python3'


	#input=hg18s_results
	#chrs=/projects/03/genereg/projects/SOLEXA/chr_info/hg18.chrom.sizes	
	#fwig=/projects/mbilenky/mlc/Jake/hs0988/wig/hs0988.q10.F516.wig.gz

	#SCRIPT_BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Aquire Resoureces
	reschrs=$(sh $SCRIPT_BASE/RNAgetRes.sh $species)
	res=$(eval echo $reschrs | cut -f 1 -d ',')
	repeats=$(echo $reschrs | cut -f 5 -d ',')
	exons=$res/$species"_exons"
	chrs=$(echo $reschrs | cut -f 2 -d ',')

	echo "================================="
	echo "The Resources used are:"
	echo "     Reference set: " $species
	echo "     name: " $name
	echo "     res: " $res
	echo "     repeats: " $repeats
	echo "     exons: " $exons
	echo "     chrs: " $chrs
	echo "================================="
	echo ""

# Example files
	#exons=/projects/mbilenky/resources/mm9rs_0313/mm9rs_0313_exons
	#repeats=/projects/mbilenky/resources/RepeatMasker/ForChimericSearch_mm9
	#chrs=/projects/03/genereg/projects/SOLEXA/chr_info/mm9.chrom.sizes
	#fwig=/projects/mbilenky/mlc/Jake/PGC_analysis/wig/rpkm_PGC_analysis.q10.F516.wig.gz
	#bam=/projects/mbilenky/mlc/Jake/Lorincz_PGC/bams/RNA-seq.PGC.bam


# READ ANALYIS---------------------------------------------------------

# Temporary file for analysis
	input="$PWD/chimericReadSearch.out"

# Loads BAM file and finds chimeric reads and creates stats on exon/repeat interactions
# Defines ER/DR/etc reads and aggrogates results

	echo "Running Chimeric Read Search python script"
	echo ""
	echo "$python3 $SCRIPT_BASE/chimericReadSearch.py $exons $repeats $bam tmp.bed > $input"
	echo ""
	
	$python3 $SCRIPT_BASE/chimericReadSearch.py $exons $repeats $bam tmp.bed > $input
	echo "read search complete."	
	echo "==============================="
	echo ""

# OUTPUT from chimericReadSearch.py
	# field = colNum

# Transcript
	GeneID='1'
	ExonRank='2'

# Repeat Name:Class:family
	RepeatName='3'
	RepeatCoord='4'

# Intersection 
	ER_intersect='5'
	RepeatExonic='6'  # Is repeat exonic
	Repeat_GeneID='7' # trancsripts overlapping repeat

# Chimeric Reads
	readER='8'
	readDR='9'
	readDE='10'
	readDD='11'
	readTot='12'

# Coordinates/Strand
	chrom='13'
	eStart='14'
	eEnd='15'
	rStart='16'
	rEnd='17'

	eStrand='18'
	rStrand='19'

# Other Information
	RepeatRank='20'
	UpExonStart='21'
	UpExonEnd='22'

# Threading Information
	UpTHREAD='23'
	DownThread='24'

# Parsed Coordinates (from file)
	EXON='25'
	UPEXON='26'
	REPEAT='27'
	UPSTREAM='28'

# Coverage Calculations ---------------------------------------------------

# Regions to input for Coverage Calculation
	# Put in the exon coords
	less $input | cut -f 13,14,15 > $PWD/tmp_repeat_coords

	# Put in the repeat coords
	less $input | cut -f 13,16,17 >> $PWD/tmp_repeat_coords

	# Put in upstream Exon coords
	less $input | cut -f 13,21,22 >> $PWD/tmp_repeat_coords

	# Put in upstream coords of the repeat
	less $input | cut -f 13,16,17,18 | awk '{ 
if ($4==1) 
	print $1"\t"($2-50)"\t"($2); 
else 
	print $1"\t"($3)"\t"($3+50); 
}' >> $PWD/tmp_repeat_coords

	# Sort/Unique Regions in Genome
	sort $PWD/tmp_repeat_coords | uniq > $PWD/tmp_uniq_repeat_coords


# Coverage Calculations

	# JAVA
	#J=/gsc/software/linux-x86_64/jre1.7.0_03/bin/java
	J=/gsc/software/linux-x86_64/jre-1.6.0_16/bin/java
	JAVA_BASE=/projects/03/genereg/projects/SOLEXA/lib/

	# Run Regions Coverage Calculator
	# Exons / Repeats / Upstream Repeats

	$J -jar -Xmx5G $JAVA_BASE/RegionsCoverageFromWigCalculator.jar -w $fwig -r $PWD/tmp_uniq_repeat_coords -s $chrs -o $PWD -n results


# Integrate data into spreadsheet ---------------------------------------------

# From Chimeric Read Search output
# parse coordinates for Exon/UpExon/Repeat/50bp

less $input | awk '{
if ($18==1) # +strand
	print $0"\tchr"$13":"$14"-"$15"\tchr"$13":"$21"-"$22"\tchr"$13":"$16"-"$17"\tchr"$13":"($16-50)"-"($16); 
else # -strand
	print $0"\tchr"$13":"$14"-"$15"\tchr"$13":"$21"-"$22"\tchr"$13":"$16"-"$17"\tchr"$13":"($17)"-"($17+50); 
}' > $PWD/tmp_input_with_key


# From Region Coverage Calculator
# parse/sort coordinates and coverage calculations
	# Col 1: Region Coordinates
	# Col 2: Region RPKM
	# Col 3: Region Max

	less tmp_uniq_repeat_coords.results.coverage | awk ' { print $1":"$2"-"$3"\t"$4"\t"$5 } ' | sort -k 1 > $PWD/tmp_sorted_repeat_data


# For each of Exon / Repeat / Upstream
	# Sort by element coordinates
	# Join <Element_line> <Element_coverage>
	# (note: join requires sorting)
	# Iterate through Exon/Repeat/Upstream and append on coverages
	
# Exon Elements
sort -k $EXON $PWD/tmp_input_with_key > $PWD/tmp_input_sorted0
join -1 $EXON -2 1 $PWD/tmp_input_sorted0 $PWD/tmp_sorted_repeat_data | sed -e 's/ /\t/g' > $PWD/tmp_results_0

# Upstream Exon Elements
sort -k $UPEXON $PWD/tmp_results_0 > $PWD/tmp_input_sorted1
join -1 $UPEXON -2 1 $PWD/tmp_input_sorted1 $PWD/tmp_sorted_repeat_data | sed -e 's/ /\t/g' > $PWD/tmp_results_1

# Repeat Elements
sort -k $REPEAT $PWD/tmp_results_1 > $PWD/tmp_input_sorted2
join -1 $REPEAT -2 1 $PWD/tmp_input_sorted2 $PWD/tmp_sorted_repeat_data | sed -e 's/ /\t/g' > $PWD/tmp_results_2

# Upstream Elements
sort -k $UPSTREAM $PWD/tmp_results_2 > $PWD/tmp_input_sorted3
join -1 $UPSTREAM -2 1 $PWD/tmp_input_sorted3 $PWD/tmp_sorted_repeat_data | sed -e 's/ /\t/g' > $PWD/tmp_final

# Add header to final results


# Parse Final Results -----------------------------------------------

# Add headers
	echo -en "transcriptID\texonRankInTranscript\trepeatName\tcoordinates\tER_Interaction\tIsExonic\tExonsOverlappingWithRepeat\t" > $PWD/final_results

	echo -en "ER\tDR\tDE\tDD\tTotal\tChromosome\tEStart\tEEnd\tRStart\tREnd\tEStrand\tRStrand\tRepeatRank\tUpExonStart\tUpExonEnd\tUpThread\tDownThread\t" >> $PWD/final_results

	echo -e "ExonRPKM\tExonMax\tUpExonRPKM\tUpExonMax\tRepeatRPKM\tRepeatMaxCoverage\tUpstreamRepeatRPKM\tUpstreamRepeatMaxCoverage" >> $PWD/final_results

# Cut out the first three columns
# (Coord Keys for sort/join)
	cut -f 5- $PWD/tmp_final  >> $PWD/final_results


# Clean up
	#rm tmp*
	#rm $input 

# End of script :D
