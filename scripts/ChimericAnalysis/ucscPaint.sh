#!/bin/bash
# ucscPaint.sh
# 
# USAGE:
#	bash ucscPaint.sh <input.list>
#	ran from westLion.sh
#
#
# CONTROL PANEL -----------------------------------------------------
#

# CORE SCRIPT -------------------------------------------------------
#

# Create a wig file
#
for wigFile in $(ls | grep wig.gz)
do
	echo ------------------------------------

	#sh wigChr.sh D null $X tmp.wig.gz
	
	CHR='null'

	gzip -dc $wifFile | awk "/chrom=$CHR/{p=0}/chrom=[^$CHR]/{p=1}p" - | gzip > $tmp.wig.gz

	OUTPUT=$( echo $wigFile | sed 's/wig.gz/bw/g' )

	wigToBigWig tmp.wig.gz ~/resources/index/hg19r/hg19r.chr.size $OUTPUT

	echo -----------------------------------
done

# Junction File

# Parse Bed File
	# Remove Header
	# Sort by coordinates
	# If Score (col5) > 1000, make it 1000 (max)
	# Make it tab-delimited
	tail -n +2 $INPUT | sort -k1,1 -k2,2n - |\
	awk -F $'\t' '{ if ( $5 <= 1000 )
	print $0;
	else
	$5='1000'
	print $0
	}' |\
	sed 's/ /\t/g' - > tmp.bed


# Convert to Big Bed Ouput
	bedToBigBed tmp.bed $CHR $OUTPUT

# Clean up
	rm tmp.bed

# Transcripts File
#
