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
	echo ' Found System-specific intitializeBin file to use'
	INITBIN="initializeBin_$SYSTEM.sh"
	FILE=$INITBIN
		eval $FCHECK_x
else
	echo ' Using default binary initialization file'


	if [ -e initializeBin.sh ]
		INITBIN="initializeBin.sh" # ./BASE/bin/
		FILE=$INITBIN

	else # copy initilizeBin from scripts folder
		echo ' Using default binary initilization file'
		cp $BASE/scripts/inializeBin.sh

		INITBIN='initializeBin.sh'
		eval $FCHECK_x
fi



	source $INITBIN # Run initializeBin.sh

	echo ' ... binary check completed successfully!'
	echo ''

# Resource CHECK ------------------------------------------
cd $BASE


# Check for resource files 
	echo " ... checking resource: $INDEX "

	source $SCRIPTS/initializeRes.sh $INDEX

	echo ' ... resource check completed successfully!'
	echo ''

# WORKSPACE ----------------------------------------------

echo " ---------- Set-up Project Workspace ---------- "

# Make project folder 
# ./LIONS/projects/<Project>
echo " Initializing $PROJECT Directory"
	mkdir -p $pDIR

# Make Library Folders
for LINE in $(cat $INPUT_LIST)
do
	echo $LINE
	#mkdir -p $pDIR/$LIBRARY
	
done


# Data Demo ------------------------------------------
# NOT IMPLEMENTED ***
##
