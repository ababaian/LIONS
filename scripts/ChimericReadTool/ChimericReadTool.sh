#!/bin/bash
# Chimeric Read Tool v1.0 
# Mod by Artem
# Usage:
# ChimericReadTool.sh <input.bam> <Coverage.wig> <Reference_set>
# output goes to current working directory
# =====================================================================

#set -e
#set -o pipefail

input=$1
PWD=$2
JUMBLE=$3

fwig="$PWD/wig/K562_j$JUMBLE.q10.F772.wig.gz"

# INITIALIZATION ------------------------------------------------------

	SCRIPT_BASE='/home/ababaian/software/ChimericReadTool'
	#python3=/gsc/software/linux-x86_64/python-3.2.2/bin/python3	
	python3='/home/ababaian/software/python3/bin/python3'


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

# Exons in the Gene Model
	ExonInGene='25'

# Parsed Coordinates (from file)
	EXON='26'
	UPEXON='27'
	REPEAT='28'
	UPSTREAM='29'

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
	#JAVA_BASE=/projects/03/genereg/projects/SOLEXA/lib/
	JAVA_BASE=/home/ababaian/software/RNAseqPipeline/bin/

	# Run Regions Coverage Calculator
	# Exons / Repeats / Upstream Repeats

	$J -jar -Xmx5G $JAVA_BASE/RegionsCoverageFromWigCalculator.jar -w $fwig -r $PWD/tmp_uniq_repeat_coords -s /home/ababaian/resources/chimeric/hg19r.chr.size -o $PWD -n results


# Integrate data into spreadsheet ---------------------------------------------

# From Chimeric Read Search output
# parse coordinates for Exon/UpExon/Repeat/50bp
# added as additional columns on the right

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

	echo -en "ER\tDR\tDE\tDD\tTotal\tChromosome\tEStart\tEEnd\tRStart\tREnd\tEStrand\tRStrand\tRepeatRank\tUpExonStart\tUpExonEnd\tUpThread\tDownThread\tExonInGene\t" >> $PWD/final_results

	echo -e "ExonRPKM\tExonMax\tUpExonRPKM\tUpExonMax\tRepeatRPKM\tRepeatMaxCoverage\tUpstreamRepeatRPKM\tUpstreamRepeatMaxCoverage" >> $PWD/final_results

# Cut out the first three columns
# (Coord Keys for sort/join)
	cut -f 5- $PWD/tmp_final  >> $PWD/final_results


# Clean up
	rm tmp*
	#rm $input 

# End of script :D
