#!/bin/bash
# Usage: .lions.sh <parameter.ctrl (opt.)>
#
# ===================================================================
# LIONS analysis pipeline
# ===================================================================
#
# Analyze an input .bam RNAseq file for transcripts which initiate in
# transposable elements and create an annotation file. Compare TE
# files between biological groups.
#
# Details can be found in README
#

# INITIALIZATION ===================================================
# Start-up script which checks all requisites are operational for
# LIONS to run. Also initializes the project space  
# *** WRITE LAST ***

# Read parameter file (imports run parameters)
if [ -z $1 ]
then
	echo " No parameter input file specified. Using default file:"
	echo "      ./LIONS/parameter.ctrl"

	PARAMETER="parameter.ctrl"
else
	echo " Custom parameter file set. Using file:"
	echo "     $1"
fi
	# Run parameter script
	.$PARAMETER # works in bash only

# Run Initialization Script
	bash $SCRIPTS/initialize.sh

# EAST LION =========================================================

#sh runTH2.sh <Input Bam> <output_name>
