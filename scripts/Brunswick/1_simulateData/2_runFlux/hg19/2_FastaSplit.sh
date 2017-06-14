#!/bin/bash
# Fasta Split
# Take a multiple entry Fasta file and split it into its components
#
#
# This is required for FluxSim to run

# Go to Scripts folder for general work
#cd ~/scripts/

# ===================================================================================
# INPUT
# Fasta file
	FASTA=$1
	
# INITILIZE
	touch headline.tmp #HEADER LINES
	touch endline.tmp #END LINES
# ===================================================================================
# SCRIPT
# Line Numbers of Interest

	# Find header lines, return line numbers
	grep -n "^>" $FASTA | cut -f1 -d":" > headline.tmp
		
	# End Lines of each Fasta (HEADERLINE - 1)
	awk '{ $2 = $1 - 1; print $2 }' headline.tmp > endline.tmp
	
	# Last Line of file
	wc -l $FASTA | cut -f1 -d" " >> endline.tmp

# FASTA coords
# two column, gives [start, end] of each fasta file
	# Create coordinate file
	sed 1d endline.tmp | paste -d',' headline.tmp - > coords.tmp
	
	# Cleanup
	rm headline.tmp endline.tmp

# Make Individual FASTA files
for COORDS in $(cat coords.tmp)
do
	# Extract Header string
	HEADLINE=$(echo $COORDS | cut -f1 -d',' -)
	HEADER=$(sed -n "$HEADLINE"p $FASTA | sed 's/>//g' -)
	
	# Initialize Fasta File and populate it
	FILE="$HEADER".fa
	touch "$HEADER".fa
	sed -n "$COORDS"p $FASTA > $FILE
	
done

# ===================================================================================
rm *.tmp
