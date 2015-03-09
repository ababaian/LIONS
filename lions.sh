#!/bin/bash
# Usage: .lions.sh <parameter.ctrl (opt.)>
set -e
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
echo ''
echo ''
echo '==============================================================='
echo '========================= L I O N S ==========================='
echo '==============================================================='
echo '''                             _   _
                           _/ \|/ \_
                          /\\/   \//\
                          \|/<\ />\|/   *RAWR*
                          /\   _   /\  /
                          \|/\ Y /\|/
                           \/|v-v|\/
                            \/\_/\/
'''
echo ''

# INITIALIZATION ===================================================
# Start-up script which checks all requisites are operational for
# LIONS to run. Also initializes the project space  
# *** WRITE LAST ***

# Read parameter file (imports run parameters)
if [ -z $1 ]
then
	echo " No parameter input file specified. Importing default file:"
	echo "      ./LIONS/parameter.ctrl"
	export PARAMETER="parameter.ctrl"
	echo ''
else
	echo " Import parameter file."
	echo "    Project Parameters: ./$1"
	export PARAMETER=$1
	echo ''
fi
	# Run parameter script
	source $PARAMETER # works in bash only

# Run Initialization Script
echo ' running initializeLIONS.sh'

	bash $SCRIPTS/initializeLIONS.sh

echo ' initialization completed successfully.'
echo ''


# EAST LION =========================================================
echo ''
echo '                     E A S T       L I O N                     '
echo ''
echo ' ./LIONS/scripts/eastLion.sh' 
echo '==============================================================='
echo '  Align reads to genome and perform chimeric analysis'
echo ''
cd $pDIR #./LIONS/projects/<projectName>

# Loop through each library in input file
iterN=$(wc -l $INPUT_LIST | cut -f1 -d' ' -)

for nLib in $(seq $iterN)
do
	# Extract row of entries from input list
	rowN=$(sed -n "$nLib"p $INPUT_LIST)
	# Library Name
	libName=$(echo $rowN | cut -f1 -d' ')
	
	# ****** Add CONTROL PROTOCOL FOR CLUSTER/LOCAL RUN

	echo " Iteration $nLib: $libName ------------------------------------------"
	echo "      run: eastLion.sh $libName"
	bash $SCRIPTS/eastLion.sh $libName
done


