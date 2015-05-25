#!/bin/bash
# wigChr.sh
#----------------------------------------------------------
# Extract/Manipulate a wiggle file
# Usage:
#
# sh wigChr.sh <Delete/eXtract> <chr> <input.wig.gz> <output.wig.gz>
#
	# Delete: removes the chromsome <chr> from input.wig.gz
	# Extract: takes out the chromosome <chr> from input.wig.gz
# ---------------------------------------------------------
# PARAMETERS

# Run Delete or Extraction
	DX=$1

# Input Wiggle File
	INPUT=$3

# Output Wiggle File
	OUTPUT=$4

# Chromsome
	CHR=$2

# --------------------------------------------------------
# SCRIPT

# Remove Nulls (outdated)
	# 'gzip -dc $INPUT |
	# pcregrep -M "^fixed.*chrom=$CHR.*(\n[0-9]+)*" >> tmp

if [ $DX == 'X' ]
	then
		# Extract Chromsome
		gzip -dc $INPUT | gzip -dc $INPUT | awk "/chrom=$CHR/{p=1}/chrom=[^$CHR]/{p=0}p" - | gzip > $OUTPUT
	elif [ $DX == 'D' ]
	then
		# Delete Chromsome
		gzip -dc $INPUT | gzip -dc $INPUT | awk "/chrom=$CHR/{p=0}/chrom=[^$CHR]/{p=1}p" - | gzip > $OUTPUT
	else
		echo 'Please select [D]elete or e[X]tract'
	fi


