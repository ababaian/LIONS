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
echo "     Call Settings: $CALLSETTINGS"
echo ""

#cd $BASE #Go to base folder
 
# FILE CHECK FUNCTION----------------------------------------------

# Function: Checks $FILE exists + permissions, returns that it does
# or returns a does not exist error 2 and exits

FCHECK_rs='if [ -s $FILE -a -r $FILE ]; then echo "     $FILE found."; else echo "     $FILE not found (empty or non-readable. Check permissions.)."; echo " ===== ERROR 2: MISSING REQUISITE FILE ===== "; exit 2; fi'

FCHECK_x='if [ -s $FILE -a -x $FILE ]; then echo "     $FILE found."; else echo "     $FILE not found (empty or non-executable. Check permissions.)."; echo " ===== ERROR 2: MISSING REQUISITE FILE ===== "; exit 2; fi'

FILE='' # File to check

# create links function
lionlinkf () {
  ln -fs $1 $2 || ln -f $1 $2 || cp -f $1 $2
}

# Resource Check
echo ' ------ Run LIONS self-check procedures ------ '
echo ''

# SCRIPT CHECK --------------------------------------------
	cd $BASE/scripts #goto Script folder

	FILE='Initialize/initializeScripts.sh'
		eval $FCHECK_x

	echo ' ... checking scripts'

	bash Initialize/initializeScripts.sh # Run init script

	echo ' ... script check completed successfully!'
	echo ''

# BIN(aries) CHECK ------------------------------------------
mkdir -p $BASE/bin
cd $BASE/bin # Go to Binary folder

	echo ' ... checking binaries'

# Check for system specific initializeBin.sh scripts else use default
if [ -s initializeBin_$SYSTEM.sh ]
then
	echo '     Found System-specific intitializeBin file to use'
	INITBIN="initializeBin_$SYSTEM.sh"
	FILE=$INITBIN
		eval $FCHECK_x

	echo ' attempting to run initializeBin.sh'
	source $INITBIN # Run initializeBin.sh

else

	if [ -e initializeBin.sh ] # initializeBin.sh exists
	then
		echo '     Using .LIONS/bin binary initialization file'
		INITBIN="initializeBin.sh" # ./BASE/bin/
		FILE="$INITBIN"
	
	else # DNE: copy initializeBin.sh from scripts folder
		echo '     Using default binary initilization file'
		cp $BASE/scripts/Initialize/initializeBin.sh $lBIN/initializeBin.sh

		INITBIN='initializeBin.sh'
		FILE="$INITBIN"

		eval $FCHECK_x
	fi

	echo ' attempting to run initializeBin.sh'
	source $INITBIN # Run initializeBin.sh
	rm $INITBIN
fi

	echo ' ... binary check completed successfully!'
	echo ''

# Resource CHECK ------------------------------------------
cd $BASE

# Check for resource files 
	echo " ... checking resource: $INDEX "
	echo '     Check genome files, repeat files and annotations are in order'
	echo ''

	source $SCRIPTS/Initialize/initializeRes.sh $INDEX

	echo ' ... resource check completed successfully!'
	echo ''

# WORKSPACE ----------------------------------------------
echo " ---------- Set-up Project Workspace ---------- "
cd $BASE
mkdir -p $BASE/projects

# Make project folder 
# ./LIONS/projects/<Project>
echo " Initializing $PROJECT Directory: $pDIR"
	mkdir -p $pDIR # ./LIONS/projects/<projectName>
	echo $INPUT_LIST $PARAMETER
	mkdir -p $pDIR/logs/run$RUNID
	awk -f $SCRIPTS/Initialize/input.list.awk  $INPUT_LIST > $pDIR/input.list
	export INPUT_LIST="$pDIR/input.list" # <libName> <libPath> <group> csv file
	cp $INPUT_LIST $pDIR/logs/run$RUNID/input.list
	cp $PARAMETER $pDIR/logs/run$RUNID/parameter.ctrl

# Make and Initialize Library Folders
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

		# Bam File Path (or FastQ files)
		bamPath=$(echo $rowN | cut -f2 -d' ') 

		# Make library directory in project folder
		mkdir -p $pDIR/$libName

		# Check input files and link them to the project folder
		if [[ $bamPath == *'.bam' ]]
		then
		# Input file is bam (standard)

			if [ -s $bamPath ] #File is not empty
			then
			
				# Link the bam file to the project folder
				lionlinkf $bamPath $pDIR/$libName/input.bam
			else
				# File is empty
				echo ' ERROR 7A: Input Bam File is not readable / empty'
				echo " file: $bamPath"
				echo '  a) The input file (.bam or .fq_1 & .fq_2) isnt found'
				echo '  b) If youre using FASTQ; the input name for the read pairs'
				echo '    should be suffixed with _0'
				exit 7
			fi

		else
		# Input file should then be fastqs (comma seperated)
			fq1=$(echo $bamPath | cut -f1 -d',' - )
			fq2=$(echo $bamPath | cut -f2 -d',' - )

			# check file type (uncompressed or compressed) in input.list matches available files
			if [[ ${fq1: -3} == ".gz" ]] && [ -s ${fq1%.gz} ]
			then
				fq1=${fq1%.gz}
			elif [[ ${fq1: -3} != ".gz" ]] && [ -s ${fq1}.gz ]
			then
				fq1=${fq1}.gz
			fi 

			if [[ ${fq2: -3} == ".gz" ]] && [ -s ${fq2%.gz} ]
			then
				fq2=${fq2%.gz}
			elif [[ ${fq2: -3} != ".gz" ]] && [ -s ${fq2}.gz ]
			then
				fq2=${fq2}.gz
			fi 
 
			if [ -s $fq1 ] && [ -s $fq2 ] # Both files exist/ aren't empty
			then

				if [[ ${fq1: -3} == ".gz" ]] 
				then
					# Link both fastq files as input.fq_1 and input.fq_2
					lionlinkf $fq1 $pDIR/$libName/temp.1.fq.gz
				else
					# gzip both uncompressed fastq files as temp.1.fq.gz and temp.2.fq.gz
					gzip -c $fq1 > $pDIR/$libName/temp.1.fq.gz
				fi

				if [[ ${fq2: -3} == ".gz" ]] 
				then
					# Link both fastq files as input.fq_1 and input.fq_2
					lionlinkf $fq2 $pDIR/$libName/temp.2.fq.gz
				else
					# gzip both uncompressed fastq files as temp.1.fq.gz and temp.2.fq.gz
					gzip -c $fq2 > $pDIR/$libName/temp.2.fq.gz
				fi

		
			else
				echo ' ERROR 7B: Input File Not Accesible (Fastq)'
				echo " files: $fq1 ; $fq2"
				echo '  a) The non-bam input file (.fq_1 & .fq_2) isnt found'
				echo '  b) If youre using FASTQ; make sure youre listing two'
				echo '     files in the input.list file seperated by a comma'
				exit 7
			fi
		fi	
        fi
done

# Initialize parameter folder  ******** NOT IMPLEMENTED YET *********
	# Copy scripts + parameters of this run
	# into one place for data-keeping 

# Data Demo ------------------------------------------
# NOT IMPLEMENTED ***
##


# End of Script = D
