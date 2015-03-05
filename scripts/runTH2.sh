#!/bin/bash
# Tophat alignment pipe (BAM input)
# Moves fastq files to TMP that are input

# USAGE:
#	sh runTH2.sh <.../Input Bam> <output_name>
#	output name = <output_name>.bam
#	Input bam = </home/ababaian/data/input.bam>
#
# CONTROL PANEL --------------------------------------
# Input Name (bam)
	INPUT=$1

# Output name (no .bam suffix)
	OUTPUT=$2 # .bam aligned output
	OUTDIR="$BASE/projects/$PROJECT/$OUTPUT" # directory for final output

# Imported Parameters from <parameter.ctrl>
	# Tophat 2
	 #THREADS [-p]
	 #INREAD Inner Read distance [-r]

	# Directories
 	 #RESOURCES '.../resources/hg19r/'

	# Genome Index
	 #INDEX 'hg19r'

# CLUSTER/LOCAL ALTERNATIVE PROTOCOLS
if [ $SYSTEM == 'gsc' ]
then #Cluster
	
	# Working directories
	WORK=$TMP # work on temporary space

	# BT2 Genome Index (copy to work space)
	cp -R $RESOURCES/genome/* $WORK

	# Bam input file (cp)
	cp $INPUT $WORK

else # Local

	# Working directories
	WORK=$OUTDIR # work in output space

	# BT2 Genome Index (link to work space)
	ln -s $RESOURCES/genome/* $WORK/

	# Bam input file (link)
	ln -s $INPUT $WORK
fi 


cd $WORK # Go to work-space

# SCRIPT --------------------------------------------------
echo " runTH2.sh initialized ... "

# Declare the input/output files
	echo "The bam input file is: $INPUT"
	echo "The output filename will be $OUTPUT.bam"

# Working Directory
	echo " listing working folder contents: "
	ls -lh

# Sort the input bam file
	echo " Sorting input bam file with samtools ..."
	samtools sort -n $INPUT $TMP/temp_sort
	rm $INPUT

# Convert to fastq files
	# Produces two files for paired-end reads
	# temp.1.fq temp.2.fq
	echo " Converting sorted bam file to fastq file"
	bam2fastx -Q -q -A -P -N -o $TMP/temp.fq $TMP/temp_sort.bam
	
	# Clean Up
	rm temp_sort.bam

# Run Tophat2
	
	echo " Running tophat2 ..."
	echo "  cmd: $TOPHAT -o $PWD $INDEX $WORK/temp.1.fq $WORK/temp.2.fq"
	$TOPHAT -o $PWD $INDEX $WORK/temp.1.fq $WORK/temp.2.fq
	
	echo ' ... tophat2 completed.'

# Cleanup -------------------------------------------------
	rm temp* # Clear temporary files
	rm $INDEX* # Clear bowtie index files
	
	echo 'After clearup'
	ls -lh

# Append to single output file
	samtools cat -o unsorted.bam accepted_hits.bam unmapped.bam
	samtools sort unsorted.bam $OUTPUT
	samtools index $OUTPUT.bam

	# run a command; if bam file output is greater then 100 Mb
	# then remove fastq files and copy everything over
	# return a message the the alignment worked
	# else copy back the fastq files and output that the alignment
	# iddn't work
	if [ -s $OUTPUT.bam ] 
	then
		echo 'Alignment likley worked!'
		# Cleanup
		rm accepted_hits.bam
		rm unmapped.bam
		rm unsorted.bam
		rm *.fq
		
		# Move all remaining output to home folder
		mv $TMP $WORK/$OUTDIR
	else
		echo 'Alignment probably didnt work'
		# copy files to output
		cd /tmp
		mv $TMP $WORK/$OUTDIR
	fi

# Done script :D
echo 'Done... maybe'
