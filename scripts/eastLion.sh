#!/bin/bash
# eastLion.sh
# 
# USAGE:
#	.eastLion.sh <libName>
#	ran from lions.sh
#
# Wrapper for processing .bam file into .fastq
# and running Bowtie2 and Tophat.
# Runs RNAseqPipeline and generates a per-library
# complete lions file and filters to TE-initiated transcripts
# only.
# Output: <library.lions> file
#

# CONTROL PANEL -----------------------------------------------------
	# Run parameters are imported from the ./LIONS/parameter.ctrl

# Library Name
	libName=$1

# Input Bam File
	INPUT="input.bam"

# Output Bam File with no .bam suffix
	OUTPUT="$libName"

# Output name (no .bam suffix)
	outDir="$BASE/projects/$PROJECT/$libName" # directory for final output

# Imported Parameters from <parameter.ctrl>
	# Tophat 2
	 #THREADS [-p]
	 #INREAD Inner Read distance [-r]

	# Directories
 	 #RESOURCES='./LIONS/resources/hg19r/'

	# Genome Index
	 #INDEX='hg19r'

# CLUSTER/LOCAL ALTERNATIVE PROTOCOLS -------------------------------
if [ $SYSTEM == 'gsc' ]
then #Cluster
# Copy over files for analysis on cluster node

	# Working directories
	WORK=$TMP # work on temporary space

	# BT2 Genome Index (copy to work space)
	cp -R $RESOURCES/genome/* $WORK

	# Bam input file (cp)
	cp $INPUT $WORK

else # Local
# Create symbolic links to the output directory and work from there
	# Working directories
	WORK=$outDir # work in output space

	# BT2 Genome Index (link to work space)
	ln -s $RESOURCES/genome/* $WORK

	# Bam input file (link)
	#ln -s $INPUT $WORK/input.bam
fi 

# ===================================================================
# CORE SCRIPT========================================================
# ===================================================================

echo "     ... eastLion.sh running"
echo "     Library: $libName"
echo "     Ouput Directory: $outDir"
echo "     Working Directory: $WORK"
echo ''
cd $WORK # go to working directory


# ALIGNMENT ---------------------------------------------------------
# If an aligned file already exists, don't recalculate it
if [ -s $libName.bam ]
then
	echo " $libName.bam is already generated."
	echo " ... skipping alignment"
	echo ''
	rm $INDEX*

else # Generate Alignment

# Declare the input/output files
	echo "  No previous alignment detected"
	echo "  Aligning reads to the genome"
	echo "     Bam input: $INPUT"
	echo "     Bam output: $OUTPUT.bam"
	echo "     Genome: $INDEX"

# Working Directory ** Comment out **
	#echo " listing working folder contents: "
	#ls -lh

# Sort the input bam file
	echo " Sorting input bam file with samtools ..."
	$lBIN/samtools sort -n $INPUT $WORK/temp_sort
	#rm $INPUT # cleanup input

# Convert to fastq files
	# Produces two files for paired-end reads
	# temp.1.fq temp.2.fq
	echo " Converting sorted bam file to fastq file"
	$lBIN/bam2fastx -Q -q -A -P -N -o $WORK/temp.fq $WORK/temp_sort.bam
	
	# Clean Up
	#rm temp_sort.bam

# Run Tophat2
	echo " Running tophat2 ..."

	echo "  cmd: $lbin/tophat2 $ctrlTH2 -o $PWD $INDEX $WORK/temp.1.fq $WORK/temp.2.fq"
	$lBIN/tophat2 $ctrlTH2 -o $PWD $INDEX $WORK/temp.1.fq $WORK/temp.2.fq
	
	echo ' ... tophat2 completed.'
	echo ''

# -------- CLEANUP
	rm temp* # Clear temporary files
	rm $INDEX* # Clear bowtie index files
	
	#echo 'Directory after clearup'
	#ls -lh

# Append to single output file
	$lBIN/samtools cat -o unsorted.bam accepted_hits.bam unmapped.bam
	$lBIN/samtools sort unsorted.bam $OUTPUT
	$lBIN/samtools index $OUTPUT.bam

	# run a command; if bam file output is greater then 1 Mb
	# then remove fastq files and copy everything over
	# return a message the the alignment worked
	# else copy back the fastq files and output that the alignment
	# didn't work

	fileSize=$(du -k $OUTPUT.bam | cut -f1) # kb size of $OUTPUT.bam
	minSize=20 # Minimum size of file to be accepted

	if [ $fileSize -ge $minSize ] # bam file must be 1 Mb 
	then
		echo 'Alignment likley completed successfully!'
		# Cleanup
		rm accepted_hits.bam
		rm unmapped.bam
		rm unsorted.bam
		
	else
		echo 'Alignment probably didnt work'
		echo ' ============= ERROR 10: Alignment Not Generated ============='

		# Cluster migration upon failure
		if [ $SYSTEM == 'gsc' ]
		then
			# copy files to output
			cd $BASE
			mv $WORK $outDir
		fi
		# Else the working directory is the output directory

		exit 10 # Exit with error 10
	fi

# Post Alignment Organization
	# Create an 'alignment' folder and move all the tophat2 files into
	# this folder
	mkdir .tmp
	mv * ./.tmp
	mv .tmp alignment

	# Move input/output bam files back to main LIB directory
	mv alignment/input.bam ./input.bam
	ln -s ./alignment/$OUTPUT.bam ./$OUTPUT.bam
	ln -s ./alignment/$OUTPUT.bam.bai ./$OUTPUT.bam.bai

fi 

# ASSEMBLY ----------------------------------------------------------
# Flow control
if [ -s assembly/transcripts.gtf ]
then
# Assembly already ran; skip	
	echo "  transcripts.gtf Assembly file already generated."
	echo "  ... skipping cufflinks assembly"
	echo ''
else
# Run Cufflinks
 
# Make an ouput folder for Cufflinks Assembly
	mkdir -p assembly

# Declare the input/output files
	echo "  No assembly detected"
	echo "     Bam input: $OUTPUT.bam"
	echo "     Label: $libName"
	echo "     Genome: $INDEX"

if [ $deNovo == '1']
then
	echo '     De Novo Assembly'
else
 	echo '     Guided Assembly'
fi
	echo "     Controls:     $ctrlCL2"

# Run Cufflinks
	# Parameters set and imported from parameter.ctrl
	echo ' Running Cufflinks ...'
	echo "   cmd: $lBIN/cufflinks -o ./assembly -L $libName $ctrlCL2 $OUTPUT.bam"

	$lBIN/cufflinks -o ./assembly -L $libName $ctrlCL2 $OUTPUT.bam

	echo " ... cufflinks completed."


fi

# RESOURCE GENERATION -----------------------------------------------
# Flow Control
if [ -s resources/assembly_exons ]
then
# Resources already generated
	echo '  Resources already generated.'
	echo '  ... skipping.'
	echo ''
else
# Resources have not been generated yet
	# Initialize Resources folder for library
	mkdir -p resources
	cd resources

	# link transcripts.gtf to resources folder
	ln -fs ../assembly/transcripts.gtf transcripts.gtf

	echo "  Building resources for $libName"

	# Run buildResourceGTF.sh
	$BASE/scripts/RNAseqPipeline/resourceGeneration/buildResourceGTF.sh transcripts.gtf assembly
	cd ..
fi


# RNASEQPIPELINE ----------------------------------------------------

# For LIONS; REF='assembly'

	# RNAseq Analysis
	echo " RNAseq Pipeline Analysis"
	echo "     cmd: RNAseqPipeline.sh $REF $libName R $OUTPUT.bam"
	
	bash $SCRIPTS/RNAseqPipeline/RNAseqPipeline.sh $REF $libName R $OUTPUT.bam

	echo ""

	# Chimeric Analysis
	echo " Chimeric Reads Analysis"
	echo "      ChimericReadTool.sh $FPATH $PWD/$name/wig/$name.$QUALITY.wig.gz $REF"
	
	bash ChimericReadTool.sh $FPATH $PWD/$libName/wig/$libName.$QUALITY.wig.gz $REF
	
	echo ""

# CHIMERICREADTOOL --------------------------------------------------


# CLEAN-UP ----------------------------------------------------------

# Cluster migration to output workspace
	# Cluster migration to final output
		if [ $SYSTEM == 'gsc' ]
		then
			# copy files to output
			cd $BASE
			mv $WORK $outDir
		fi

# Done script :D
