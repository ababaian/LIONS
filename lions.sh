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
	echo "      ./LIONS/controls/parameter.ctrl"
	export PARAMETER="controls/parameter.ctrl"
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

	bash $SCRIPTS/Initialize/initializeLIONS.sh

echo ' initialization completed successfully.'
echo ''


# EAST LION =========================================================
echo ''
echo '                     E A S T       L I O N                     '
echo ''
echo ' ./LIONS/scripts/eastLion.sh ' 
echo '==============================================================='
echo '  Align reads to genome and perform TE-initiation analysis'
echo ''
cd $pDIR #./LIONS/projects/<projectName>

# Initialize Summit Log file
touch summitLog_$RUNID

# Loop through each library in input file
iterN=$(wc -l $INPUT_LIST | cut -f1 -d' ' -)

for nLib in $(seq $iterN)
do
	# Extract row of entries from input list
	rowN=$(sed -n "$nLib"p $INPUT_LIST)
        
	# test to ignore empty lines in $INPUT_LIST
        if [[ ! -z "$rowN" ]]
        then 
		# Library Name
		libName=$(echo $rowN | cut -f1 -d' ')

		echo " Iteration $nLib: $libName ------------------------------------------"
		echo "      run: $QSUB eastLion.sh $libName"

		if [ ! -e $pDIR/$libName/$libName.lions ]
		then
		# Lions output for this library doesn't exist
		# so let's make it.

			if [ $CLUSTER == '1' ]
			then # Cluster QSUB
				$QSUB $SCRIPTS/eastLion.sh $libName

			else # Local (no) QSUB
				$SCRIPTS/eastLion.sh $libName
			fi

		elif [ $SORTBYPASS = '0' ]
		then
		# Lions output already exists but
		# East Lion bypass is set to false, re-calculate lions file
		# so let's make it.

			if [ $CLUSTER == '1' ]
			then # Cluster QSUB
				$QSUB $SCRIPTS/eastLion.sh $libName

			else # Local (no) QSUB
				$SCRIPTS/eastLion.sh $libName
			fi

		else
		# East Lion file already exists and bypass is true (1)
		# Skip the east lion

			echo "   East Lions has previously been completed. "
			lionSuccess='1'
			echo $libName $lionSuccess $(date) >> $pDIR/summitLog_$RUNID
		fi

		echo " ... run complete -------------------------------------------"
		echo ''
		echo ''
	fi
done

# Check that all libraries have completed
	#iterN is the number of libraries
summitN=$(wc -l $pDIR/summitLog_$RUNID | cut -f1 -d' ' )

while [  $summitN -lt $iterN ]  # Not all EAST LION iterations have completed
do 

	# Verbose
	echo " $summitN / $iterN East Lion scripts completed. Waiting..."
	date

	# Wait 10 minutes
	sleep 600s # Actual

	# Recalculate summitN
	summitN=$(wc -l $pDIR/summitLog_$RUNID | cut -f1 -d' ')

done

# All runs completed

# Check if they are each succesful
for Log in $(cut -f2 $pDIR/summitLog_$RUNID)
do
	if [ $Log = '0' ];
	then
		echo " ERROR 15: One of the East Lion Libraries didn't finish."
		exit 15
	fi
done

# Clear summit log
rm $pDIR/summitLog_$RUNID

echo ''
echo ' All EAST LION scripts have completed. '


# WEST LION =========================================================
echo ''
echo '                     W E S T       L I O N                     '
echo ''
echo ' ./LIONS/scripts/westLion.sh ' 
echo '==============================================================='
echo '  Group and analyze lions files of TE-initiation events'
echo ''
cd $pDIR #./LIONS/projects/<projectName>

# Run West Lions Script
	echo ' run: westLion.sh'
	bash $SCRIPTS/westLion.sh

echo ''
echo ' WEST LION scripts have completed. '

