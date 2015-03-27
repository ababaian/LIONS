#!/bin/bash
# westLion.sh
# 


# CONTROL PANEL -----------------------------------------------------

echo "     ... westLion.sh running"
echo "     Library: $libName"
echo "     Ouput Directory: $outDir"
echo "     Working Directory: $WORK"
echo ''

# ===================================================================
# CORE SCRIPT========================================================
# ===================================================================
cd $pDIR #./LIONS/projects/<projectName>

# Master Flow Control
## **** To be implemented *****


# PROJECT LIONS FILE ------------------------------------------------

# Initialize a project-wide .lions file: <pLionFile>
	# If it exists, append RUNID
if [ ! -e $PROJECT.lions ]
then
	#No <Project>.lions file exists; initialize it
	pLionFile="$PROJECT.lions"

	echo " $PROJECT.lions does not exist, creating it."

	# Initialize file
	touch $pLionFile

else
	# <Project>.lions exists; initialize new file with RUNID
	pLionFile="$PROJECT.$RUNID.lions"
	
	echo " Previous project file exists."
	echo " Initializing new file $pLionFile"

	touch $pLionFile
fi


# Populate Master LIONS file for this project/run
# Loop through each library in input file
iterN=$(wc -l $INPUT_LIST | cut -f1 -d' ' -)

for nLib in $(seq $iterN)
do
	# Extract row of entries from input list
	rowN=$(sed -n "$nLib"p $INPUT_LIST)
	# Library Name
	libName=$(echo $rowN | cut -f1 -d' ')

	# Most Recently generated lions file
	recentLion=$(ls -t $libName/*.lions | sed -n 1p)
	
	# If first loop; initialize master output file
	if [ $nLib == '1' ]
	then
		#From the output file generate a header file for the output
		head -n1 $recentLion > $pLionFile
	fi

	# Append filtered list to master list
	echo "     append $libName.lions to $pLionFile"
	sed '1d' $recentLion >> $pLionFile

done

	echo ''

# ChimGroup ---------------------------------------------------------
#



# End of Script *<=D
