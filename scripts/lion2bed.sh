#!/bin/sh
#
# lion2bed.sh < input.lion file> < output.bed file >
# 
# Conversion of `.lion` output to `.bed` format
# Note: `.bed` is a bit limited in the data it contains,
# this is useful for an overview of the data in a genome
# browser.

INPUT=$1     # .lion file input
OUTPUT=$2    # .bed file output

HEADER='T'   # Include Header? <T/F>

## BED6 Format
#
# 1. chrom
# 2. chromStart
# 3. chromEnd
# 4. name
# 5. score
# 6. strand
#

## lions Format (extraction)
#
# 1. chrom
# 2. Exon/TE Start
# 3. Exon/TE End
# 4. RepeatID
# 5. Total Read Fragment Count
# 6. Exon Strand
#

echo "Converting the lions file: $INPUT"
echo "  to the bed6-format file: $OUTPUT"
echo ""

# BED header
if [ "$HEADER" = "T" ]
then
	echo "track name=\"$INPUT\" description=\"LION file in BED6\" visibility=2 colorByStrand=\"255,0,0 0,0,255\"" > $OUTPUT
else
	touch $OUTPUT
fi

# Coordinates
# Cut coord field; exclude header
cut -f4 $INPUT    | tail -n +2  | \
sed 's/[:-]/\t/g' - > bed1.tmp

# Cut repeatName and Total Reads
# append columns to output
cut -f3,12 $INPUT | tail -n +2  | \
paste bed1.tmp -    > bed2.tmp

# Cut exon strand, convert format
# and append
cut -f18 $INPUT   | tail -n +2  | \
sed 's/\-1/-/g' -  | sed 's/0/\./g' - | \
sed 's/1/+/g'  -  | paste bed2.tmp - > bed3.tmp


mv bed3.tmp $OUTPUT
rm *.tmp

# end
