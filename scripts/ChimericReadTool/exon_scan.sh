#!/bin/bash
# Exon Scan
#
# Input an <geneSet_exon> file
# returns <geneSet_exon2> file
#
# with an additional column (column 9) with
# the maximum number of exons per transcriptid
#
# This has to be one of the ugliest scripts
# ever written!


# Parameters
	INPUT=$1

LINE_X='Line_to_Modify____index_X'

#Initilaized Paramters
LINE_Y='Trailing_line_____index_X-1'
LY_transcript='Trailing_transcriptID'

# Hold-over
HOLDER='FALSE'
FIRSTLINE='TRUE'

# Initialized Output
OUTPUT="$INPUT"_2
echo "$OUTPUT is being generated"


# Column Count
	# If there are 8 column; then biotypes are used
	# If there are 7 columns; then biotypes not used; assign as NA
COLUMN8=$(sed -n 1p $INPUT | cut -f8 -)

	if [ "$COLUMN8" == '' ]
	then
		#echo "No Biotypes assigned; use NA"

		MULTI='NA 2'
		LONE='NA 1'
	else
		MULTI='2'
		LONE='1'
	fi

sed -n 1p $INPUT > tmp.header
echo $MULTI | sed 's/ /\t/g' - | paste tmp.header - > $OUTPUT

rm tmp.header

while read p; do
	if [ "$HOLDER" == 'TRUE' ]
	then
		RankAssign=$MULTI
	fi
	HOLDER='FALSE'

	# Leading Line of the file 
	
	LINE_X=$(echo $p)
	LX_transcript=$(echo $p | cut -f2 -d' ' -)
	LX_rank=$(echo $p | cut -f7 -d' ' - )

	# If transcriptID is the same	
	if [ "$LX_transcript" == "$LY_transcript" ]
	then
		# Transcript has greater then 1 exon
		RankAssign=$MULTI
		HOLDER='TRUE'
	fi

	if [ "$FIRSTLINE" == 'TRUE' ]
	then
		#Do nothing
		FIRSTLINE='FALSE'
	else
		echo $LINE_Y $RankAssign | sed 's/ /\t/g' - >> $OUTPUT
	fi

	# Move Index Line // Reset
	LINE_Y=$LINE_X
	LY_transcript=$LX_transcript
	LY_rank=$LX_rank
	RankAssign=$LONE

done < $INPUT

# End of Script :D
