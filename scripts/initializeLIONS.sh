#!/bin/bash
# initializeLIONS.sh
set -e # quit script if errors occur
#
# 1) Check neccesary bin/scripts exist and are accesible
# 2) Create project workspace
# 3) (optional/ not implemented) run toy data to ensure pipeline works

echo ""
echo " ./LIONS/scripts/initilizeLIONS.sh ..."
echo " =============================================================="
echo "     Project Name: $PROJECT "
echo "     Run Identification Number: $RUNID"
echo "     Library List: $INPUT_LIST"
echo "     Genome: $INDEX "
echo "     System: $SYSTEM"
echo "             cores: $THREADS"
echo "             qsub: $QSUB (if applicable)"
echo "     LIONS base dir: $BASE"
echo ""

#cd $BASE #Go to base folder
 
# FILE CHECK FUNCTION----------------------------------------------

# Function: Checks $FILE exists + permissions, returns that it does
# or returns a does not exist error 2 and exits

FCHECK_rs='if [ -s $FILE -a -r $FILE ]; then echo "     $FILE found."; else echo "     $FILE not found (empty or non-readable)."; echo " ===== ERROR 2: MISSING REQUISITE FILE ===== "; exit 2; fi'

FCHECK_x='if [ -s $FILE -a -x $FILE ]; then echo "     $FILE found."; else echo "     $FILE not found (empty or non-executable)."; echo " ===== ERROR 2: MISSING REQUISITE FILE ===== "; exit 2; fi'

FILE='' # File to check

# Resource Check
echo ' ------ Run LIONS self-check procedures ------ '
echo ''

# SCRIPT CHECK --------------------------------------------
	cd $BASE/scripts #goto Script folder

	FILE='initializeScripts.sh'
		eval $FCHECK_x

	echo ' ... checking scripts'

	bash initializeScripts.sh # Run init script

	echo ' ... script check completed successfully!'
	echo ''

# BIN(aries) CHECK ------------------------------------------
cd $BASE/bin # Go to Binary folder

	echo ' ... checking binaries'

# Check for system specific initializeBin.sh scripts else use default
if [ -s initializeBin_$SYSTEM.sh ]
then
	echo '     Found System-specific intitializeBin file to use'
	INITBIN="initializeBin_$SYSTEM.sh"
	FILE=$INITBIN
		eval $FCHECK_x
else

	if [ -e initializeBin.sh ] # initializeBin.sh exists
	then
		echo '     Using .LIONS/bin binary initialization file'
		INITBIN="initializeBin.sh" # ./BASE/bin/
		FILE="$INITBIN"
	
	else # DNE: copy initializeBin.sh from scripts folder
		echo '     Using default binary initilization file'
		cp $BASE/scripts/initializeBin.sh $lBIN/initializeBin.sh

		INITBIN='initializeBin.sh'
		FILE="$INITBIN"

		eval $FCHECK_x
	fi
fi

	echo ' attempting to run initializeBin.sh'
	source $INITBIN # Run initializeBin.sh
	rm $INITBIN

	echo ' ... binary check completed successfully!'
	echo ''

# Resource CHECK ------------------------------------------
cd $BASE

# Check for resource files 
	echo " ... checking resource: $INDEX "
	echo '     Check genome files, repeat files and annotations are in order'
	echo ''

	source $SCRIPTS/initializeRes.sh $INDEX

	echo ' ... resource check completed successfully!'
	echo ''

# WORKSPACE ----------------------------------------------
echo " ---------- Set-up Project Workspace ---------- "
cd $BASE

# Make project folder 
# ./LIONS/projects/<Project>
echo " Initializing $PROJECT Directory: $pDIR"
	mkdir -p $pDIR # ./LIONS/projects/<projectName>
	echo $INPUT_LIST $PARAMETER
	mkdir $pDIR/run$RUNID
	cp $INPUT_LIST $pDIR/run$RUNID/input.list
	cp $PARAMETER $pDIR/run$RUNID/parameter.ctrl

# Make and Initialize Library Folders
	iterN=$(wc -l $INPUT_LIST | cut -f1 -d' ' -)

for nLib in $(seq $iterN)
do
	# Extract row of entries from input list
	rowN=$(sed -n "$nLib"p $INPUT_LIST)

	# Library Name
	libName=$(echo $rowN | cut -f1 -d' ')

	# Bam File Path
	bamPath=$(echo $rowN | cut -f2 -d' ')
	
	# Make library directory in project folder
	mkdir -p $pDIR/$libName

	# Link the bam file to the project folder
	ln -fs $bamPath $pDIR/$libName/input.bam
	
done

# Initialize parameter folder  ******** NOT IMPLEMENTED YET *********
	# Copy scripts + parameters of this run
	# into one place for data-keeping 

# Data Demo ------------------------------------------
# NOT IMPLEMENTED ***
##


# End of Script = D
