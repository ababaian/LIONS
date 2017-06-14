#! /bin/bash
# chimIntersect.sh
# Usage: sh chimIntersect.sh <Chimeric_output> <Assembly_GTF> <Gene_UCSC>
# MODIFIED FOR USE WITH LIONS PIPELINE#
#
# Intersects Chimeric_output transcripts to RefSeq_UCSC output
# to return overlapping gene names.
#
# <LIBRARY.lcsv> standard output from ChimericReadTool.sh 
# <Assembly.GTF> is the matched cufflinks assembly GTF file
# <refGene.bed> is the refGene table downloaded from UCSC table

# Extended --------
#
# From refGene parse the table into a bed file
# to define 'genic' and 'intergenic' regions.
#     <1:Chromosome> <2:Start> <3:End> <4:name> <5:score> <6:strand>
#
#
# From Chimeric Output take transcriptID and search Assembly.gtf
#  to find transcript information.
#  Parse transcript information in ordered Bed File
#     <1:Chromosome> <2:Start> <3:End> <4:name> <5:score> <6:strand>
#  Order is equivelant to interaction order in Chimeric_output.
#
# [ More efficient to look this information up in generating
# chimeric output (ChimericReadTool.sh) and then use it here]
#
# Intersect refGene.bed & chimeric.bed files
# Pulling out Name(s) from refGene for each chimeric.bed transcript
# or returning 'Intergenic' for non-overlapping transcripts.
#

# Refseq 1:		----------------------
# Refseq 2:         -------------------
# Assembly 1:         ----------
# Assembly 2: ----

# Here transcript Assembly_1 overlaps Refseq_1,2
# Assembly 2 is 'intergenic'
# NOTE: exon structure not considered, simply genic/intergenic
# -----------------

# INITIALIZE ========================================================

	# Library Name Input
	libName=$1

	# Chimeric Read Tool output raw .lcsv file
	chimFile="$libName.lcsv"

	# Cufflinks output
	assemblyFile="assembly/transcripts.gtf"

	# Reference gene set (BED file)
	refFile="$RESOURCES/annotation/refGene.bed"

	# Output 
	OUTPUT="$libName.pc.lcsv"

	# Check input
	
	if [ -n "$refFile" ]
		then
			echo 'Running chimIntersect.sh script ---------------'
			echo "    Chimeric file: $1"
			echo "    Assembly file: $assemblyFile"
			echo "    Gene Ann File: $refFile"
			echo " "
		else
			echo "Missing parameters!"
			echo "Usage: sh chimIntersect.sh <Chimeric_output> <Assembly_GTF> <Gene_UCSC>"
	
			exit 1
	fi

# CHIMERA TRANSCRIPT LIST ---------------------------------
# Parse transcriptID from chimeric output for lookup
# TranscriptID = column 1b
# [Gene:Transcript]

	cut -f1 $chimFile > chim1.tmp
	cut -f2 -d':' chim1.tmp | sed 1d - > chimID.tmp 


# ASSEMBLY BED FILE ----------------------------------------
# Chimeric File parse
# Grab Transcript Coordinates and drop Exon data from GTF
	grep -P "\ttranscript\t" $assemblyFile |
	sed 's/[";]/ /g' - | sed 's/ /\t/g' - |
	#cut -f1,4,5,7,16,21 > ass1.tmp
	awk 'BEGIN {OFS="\t"}; {print $1,$4,$5,$12,$14,$7}' - |
	sort -u - > assBed.tmp

# Column Order
	# chr
	# txStart
	# txEnd
	# name (TranscriptID from Assembly)
	# score (FPKM)
	# strand

# Drop to R and perform look up on ChimID into AssemblyBed
# and return the unique AssemblyBed entries from ChimID
	ln -s $refFile refGene.bed
	$lBIN/Rscript $SCRIPTS/ChimericReadTool/chimIntLookup.R
	rm refGene.bed
	rm *.int.refGene
	# Final output is named 'appendCol.tmp'

# Merge Data
paste $chimFile appendCol.tmp > $OUTPUT

# Clean up
rm *.tmp
rm *.bed

# End of script *<:P
