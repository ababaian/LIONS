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
# Output: <library.lion> file
#

# CONTROL PANEL -----------------------------------------------------
	# Run parameters are imported from the ./LIONS/parameter.ctrl

# Library Name
	export libName=$1

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


# MASTER FLOW CONTROL ===============================================
if [ ! -e $outDir/$libName.lion ] # no primary LIONS file exists
then
# RUN LIONS PIPELINE FROM THE START

# CLUSTER/LOCAL ALTERNATIVE PROTOCOLS -------------------------------
if [ $SYSTEM == 'gsc' ]
then #Cluster
# Copy over files for analysis on cluster node

	# Working directories
	export WORK=$TMP # work on temporary space

	# BT2 Genome Index (copy to work space)
	cp -R $RESOURCES/genome/* $WORK

	# Check if there is an input fastq file or bam file
	if [ ! -s $outDir/temp.1.fq.gz ] && [ ! -s $outDir/temp.2.fq.gz ]
	then
		inputType='bam'
		# Bam input file (cp)
		cp $outDir/$INPUT $WORK

		# Create a temporary symbolic link to the input
		ln -s $outDir/$INPUT $outDir/$INPUT.tmp
	else
		inputType='fastq'
		# Fastq input file
		cp $outDir/temp.1.fq.gz $WORK
		cp $outDir/temp.2.fq.gz $WORK
	fi
	


	# Running on temporary space
	echo " Running on temporary space on cluster"
	echo " ls -alh"
	ls -alh $WORK

else # Local
# Create symbolic links to the output directory and work from there
	# Working directories
	export WORK=$outDir # work in output space

	# BT2 Genome Index (link to work space)
	ln -sf $RESOURCES/genome/* $WORK
	# Bam input file (link)
	#ln -s $INPUT $WORK/input.bam
fi 

# ===================================================================
# CORE SCRIPT========================================================
# ===================================================================

echo "     ... eastLion.sh running"
echo "     Library: $libName"
echo "     Output Directory: $outDir"
echo "     Working Directory: $WORK"
echo "     Alignment Bypass: $ALIGNBYPASS"
echo ''
cd $WORK # go to working directory


# ALIGNMENT ---------------------------------------------------------
# Use Tophat2 to generate an alignment with the parameters defined
# in paramters.ctrl file
# LIONS defaults are tuned for aligning reptative sequences; feel free to
# change it up and try new runs

# Alignment Flow Control
# If an aligned file already exists, don't recalculate it

if [ -s alignment/$libName.bam ]
then
	echo " $libName.bam is already generated."
	echo " ... skipping alignment"
	echo ''
	rm $INDEX*
	
	# Check if index is generated
	if [ ! -s $libName.bam.bai ]
	then
		$lBIN/samtools index $OUPUT.bam
		mv $OUTPUT.bam.bai alignment/$OUTPUT.bam.bai
		ln -s alignment/$OUTPUT.bam.bai $OUTPUT.bam.bai
	fi	

else # Generate Alignment

# Alignment Bypass Flow Control Starts
if [ $ALIGNBYPASS == '1' ] && [ -s $INPUT ] # Bypass is True
then
	# Create a symbolic link between the input bam file
	# and the final output.bam file
	echo "  No previous alignment detected"
	echo "  Aligment Bypass is set to true"
	echo "  A new alignment won't be calculated,"
	echo "  The input alignment will be used instead"

	ln -s $(readlink -f $INPUT) $PWD/$OUTPUT.bam
	$lBIN/samtools index $OUTPUT.bam

	# Clean up index files
	rm $INDEX*

else # Alignment Bypass is False, calculate Alignment
# Declare the input/output files
	echo "  No previous alignment detected"
	echo "  Aligning reads to the genome"
	echo "     Bam (or fq) input type: $inputType"
	echo "     Bam output: $OUTPUT.bam"
	echo "     Genome: $INDEX"

if [ ! -s temp.1.fq.gz ] && [ ! -s temp.2.fq.gz ]
then
# if the input is a bam file generate fastq files
# else fq files exist and will be used

# Sort the input bam file for generating fastq files for input
	echo " Sorting input bam file with samtools ..."
	$lBIN/samtools sort -n $INPUT $WORK/temp_sort
	#rm $INPUT # cleanup input

# Convert bam to fastq files
	# Produces two files for paired-end reads
	# temp.1.fq temp.2.fq
	echo " Converting sorted bam file to fastq file"
	$lBIN/bam2fastx -Q -q -A -P -N -o $WORK/temp.fq $WORK/temp_sort.bam
	gzip temp.1.fq
	gzip temp.2.fq
	
	# Clean Up from this point
	rm temp_sort.bam
fi

# Run Tophat2 Alignment
	echo " Running tophat2 ..."

	echo "  cmd: $lBIN/tophat2 $ctrlTH2 -o $PWD $INDEX $WORK/temp.1.fq.gz $WORK/temp.2.fq.gz"
	$lBIN/tophat2 $ctrlTH2 -o $PWD $INDEX $WORK/temp.1.fq.gz $WORK/temp.2.fq.gz

# Alignment Sanity Check I
	# run a command; if bam file output is greater then 1 Mb
	# then remove fastq files and copy everything over
	# return a message the the alignment worked
	# else copy back the fastq files and output that the alignment
	# didn't work

	fileSize=$(du -k accepted_hits.bam | cut -f1) # kb size of $OUTPUT.bam
	minSize=20 # Minimum size of file to be accepted

	if [ $fileSize -ge $minSize ] # bam file must be 1 Mb 
	then
		echo 'Alignment likely completed successfully!'
	else
		echo "Alignment probably didn't work"
		echo ' ============= ERROR 10: Alignment Not Generated ============='

	lionsSuccess='0'
	echo $libName $lionSuccess $(date) >> $pDIR/summitLog_$RUNID

	exit 10 # Exit with error 10
	fi

	echo ' ... tophat2 completed.'
	echo ''

# -------- Post Alignment Cleanup
	rm temp* # Clear temporary files
	rm $INDEX* # Clear bowtie index files
	
	#echo 'Directory after clearup'
	#ls -lh

# Append tophat2 output bam files to a to single bam file
# (iff the bam file is less than 20G)
	##$lBIN/samtools cat -o unsorted.bam accepted_hits.bam unmapped.bam
	##$lBIN/samtools sort unsorted.bam $OUTPUT

fileSize=$(du -m $(readlink -f accepted_hits.bam ) | cut -f1) # Mb size

if [ $fileSize -ge 20000 ] # Bam > 20gb
	then
	# Rename accepted hits to output.bam
	mv accepted_hits.bam $OUTPUT.bam
	mv unmapped.bam $OUTPUT.unmapped.bam
	$lBIN/samtools index $OUTPUT.bam
	else
	# Merge mapped/unmapped bam files and rename to output.bam
	$lBIN/samtools merge -r $OUTPUT.bam accepted_hits.bam unmapped.bam
	$lBIN/samtools index $OUTPUT.bam
fi

fi # Alignment Bypass Flow Control Ends


# Calculate Flagstat for bam file (Read statistics)
	$lBIN/samtools flagstat $OUTPUT.bam > $OUTPUT.flagstat

# Alignment Sanity Check II

	fileSize=$(du -k $(readlink -f $OUTPUT.bam) | cut -f1) # kb size of $OUTPUT.bam
	minSize=20 # Minimum size of file to be accepted

	if [ $fileSize -ge $minSize ] # bam file must be 1 Mb 
	then
		echo 'Bam processing succesful'
		# Cleanup
		rm -f accepted_hits.bam
		rm -f unmapped.bam
		rm -f unsorted.bam
		
	else
		echo "Alignment probably didn't work"
		echo ' ============= ERROR 10: Alignment Not Generated ============='

		# Cluster migration upon failure
		if [ $SYSTEM == 'gsc' ]
		then
			# copy files to output
			cd $BASE
			mv $WORK $outDir
		fi
		# Else the working directory is the output directory

		lionsSuccess='0'
		echo $libName $lionSuccess $(date) >> $pDIR/summitLog_$RUNID

		exit 10 # Exit with error 10
	fi

# Post Alignment Organization of project folder
	# Create an 'alignment' folder and move all the tophat2 files into
	# this folder
	mkdir .tmp
	mv * ./.tmp
	mv .tmp alignment

	# Move input/output bam files back to main LIB directory
	mv alignment/input.bam ./input.bam
	ln -s ./alignment/$OUTPUT.bam ./$OUTPUT.bam
	ln -s ./alignment/$OUTPUT.bam.bai ./$OUTPUT.bam.bai

fi # End all Alignment flow

# ASSEMBLY ----------------------------------------------------------
# Use Cufflinks to generate an assembly of transcripts from which
# TE-Assembly interactions may be inferred
# Previously we used RefSeq instead of assembly (and it's still possible)
# but found the results to be poorer.

# Assembly Flow control
if [ -s assembly/transcripts.gtf ]
then
# Assembly already ran; skip	
	echo "  transcripts.gtf Assembly file already generated."
	echo "  ... skipping cufflinks assembly"
	echo ''
else


# Make an ouput folder for Cufflinks Assembly
	mkdir -p assembly

# Declare the input/output files
	echo "  No assembly detected"
	echo "     Bam input: $OUTPUT.bam"
	echo "     Label: $libName"
	echo "     Genome: $INDEX"

if [ $deNovo == '1' ]
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

	$lBIN/cufflinks -q -o ./assembly -L $libName $ctrlCL2 $OUTPUT.bam

	echo " ... cufflinks completed."


fi

# RESOURCE GENERATION -----------------------------------------------
# Create resource files from the assembly. Most of these are not neccesary
# but were/are part of the RNAseqPipeline software from the GSC.
# Thus 

# Resource Flow Control
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

	# From FlagStat File; create a file of 'mappedReads' and put it into resources
	sed -n 3p alignment/$libName.flagstat | cut -f1 -d' ' - > resources/mappedReads
	
fi # End resource flow control

# RNASEQPIPELINE ----------------------------------------------------
# Generates RPKM values for the transcripts
# Creates a wig file ** Edit this for UCSC compatability **

# Expression Flow Control
if [ -s expression/wig/$libName.$QUALITY.wig.gz ]
then
	# RNA seq already ran
	echo " RNAseq Pipeline already performed"
	echo " .... skipping "

else
# Note: For LIONS; REF='assembly'

	# RNAseq Analysis
	echo " RNAseq Pipeline Analysis"
	echo "     cmd: RNAseqPipeline.sh $REF $libName R $OUTPUT.bam"
	
	bash $SCRIPTS/RNAseqPipeline/RNAseqPipeline.sh $REF $libName R $OUTPUT.bam

	echo ""

	# Wig Sanity Check
	echo expression/wig/$libName.$QUALITY.wig.gz

	if [ -s expression/wig/$libName.$QUALITY.wig.gz ] # was the wig file generated
	then
		echo "Wig file generated successfully."
		# Continue
	else
		echo "Wig file not generated. Errors are afoot."
		echo ' ============= ERROR 12: Wig Not Generated ============='
		lionSuccess='0'
		echo $libName $lionSuccess $(date) >> $pDIR/summitLog_$RUNID
		exit 12
		
	fi
	echo ""

fi # End expression flow control

# CHIMERICREADTOOL --------------------------------------------------
# Core protocol for detecting TE-Transcript interactions
# This run script sets up the file architecture 
# and runs the python script to generate the core values for the .lion
# files

# Chimera Read Tool Flow Control
if [ -s $libName.pc.lcsv ];
then
	# RNA seq already ran
	echo " Chimeric Read Tool already ran"
	echo " .... skipping "
	echo ''

else

	# Chimeric Analysis
	echo " Chimeric Reads Analysis"
	echo "      ChimericReadTool.sh $libName.bam"

	#bash $SCRIPTS/ChimericReadTool/ChimericReadTool.sh $OUTPUT.bam $pDIR/$libName/expression/wig/$libName.$QUALITY.wig.gz $REF

	bash $SCRIPTS/ChimericReadTool/ChimericReadTool.sh $OUTPUT.bam expression/wig/$libName.$QUALITY.wig.gz $REF
	# Output is <libName>.lcsv raw file
	# Total TE-Transcript interaction table

	# ChimIntersect
	# Intersection with Defined Protein Coding Genes
	echo "  ChimIntersect"
	echo "       chimIntersect.sh $libName"
	bash $SCRIPTS/ChimericReadTool/chimIntersect.sh $libName	

	echo " ... Chimeric analysis completed."
	echo ""
fi # End chimera read tool flow


# CHIMSORT ----------------------------------------------------------
# Sort the pcRaw output from ChimericReadTool/chimIntersect into
# a set of TE which are likely to be the initiation events
# the sort script is customizable and parameters are set in the
# parameter.ctrl file

# Read Mapped Read number
	mappedReads=$(cat resources/mappedReads)

	# Run ChimSort Script
	echo ' ChimSort'
	echo "     Library: $libName"
	echo "     Raw Lions File: $libName.pc.lcsv" 
	echo "     Output: $libName.lion"
	echo ''
	echo "     --- classification parameters --- "

	if [ $CALLSETTINGS == 'custom' ]
	then 
		# A custom setting was defined and will be printed out
		echo "     Mapped Reads: $mappedReads"
		echo "     # Read Support = $crtReads"
		echo "     ThreadRatio = $crtThread"
		echo "     DownStream Threads = $crtDownThread"
		echo "     RPKM = $crtRPKM"
		echo "     TE Contribution = $crtContribution"
		echo "     Upstream Coverage = $crtUpstreamCover"
		echo "     Upstream Exon Expression = $crtUpstreamExonCover"
	else
		# A pre-defined setting was used
		echo "     Predefined sort settings: $CALLSETTINGS"
	fi

	echo ''

	# use ANN-based classifier or custom classifier
	if [ $CALLSETTINGS = 'transcriptomeANN' ]
	then
		# Chimeric Filtering use ANN model
		echo " Run ChimANNSort"
		echo "     Rscript chimANNSort.R $libName.pc.lcsv $libName.ann.lion $ANNMODEL"

		$lBIN/Rscript $SCRIPTS/ChimericReadTool/chimAnnSort.R $libName.pc.lcsv $libName.ann.lion $ANNMODEL

	else	
		# Chimeric Filtering using thresholds
		echo "  Run ChimSort"
		echo "     Rscript chimSort.R $libName.pc.lcsv $libName.lion $mappedReads $CRT"
	
		$lBIN/Rscript $SCRIPTS/ChimericReadTool/chimSort.R $libName.pc.lcsv $libName.lion $mappedReads $CRT

		# Sanity Check
		if [ -s $libName.lion ] # was the lions file generated
		then
			echo "Lions file generated successfully."
			lionSuccess='1'
		else
			echo "Lions file not generated. Something is amiss."
			lionSuccess='0'
		fi

	fi
	echo ""

# CLEAN-UP ----------------------------------------------------------

# Cluster migration to output workspace
	# Cluster migration to final output
		if [ $SYSTEM == 'gsc' ]
		then
			# discard input bam file (which was copied to tmp)
			rm $WORK/$INPUT
			# copy files to output
			cd $BASE
			mv $WORK/* $outDir
		fi

# End LIONS pipeline

# MASTER FLOW CONTROL ===============================================
else # LIONS file already exists

	# Recalculate chimSort (Bypass Protocol)
	if [ $SORTBYPASS = '1' ] 
	then
	# East Lion already complete; don't recalculate sort
		echo " East Lions has already been completed. "
		echo "    $libName.lion exists"
		echo "    not re-calculating chimSort (SORTBYPASS = 1)"

	else
	# East Lion already completed; re-calcualte sort though
		echo " East Lions has already been completed. "
		echo "    $libName.lion exists"
		echo "    SORTBYPASS = 0, "
		echo "    will re-recalcualte chimSort ..."

		# Move to project/library directory
		WORK=$outDir # work in output space
		cd $WORK # go to working directory

		# Read Mapped Read number
		mappedReads=$(cat resources/mappedReads)

		# Run ChimSort Script
		echo ' Re-calculate ChimSort'
		echo "     Library: $libName"
		echo "     Raw Lions File: $libName.pc.lcsv" 
		echo "     Output: $libName.lion"
		echo "     Run ID: $RUNID"
		echo ''
		echo "     --- Classification parameters --- "

		if [ $CALLSETTINGS == 'custom' ]
		then 
			# A custom setting was defined and will be printed out
			echo "     Mapped Reads: $mappedReads"
			echo "     # Read Support = $crtReads"
			echo "     ThreadRatio = $crtThread"
			echo "     DownStream Threads = $crtDownThread"
			echo "     RPKM = $crtRPKM"
			echo "     TE Contribution = $crtContribution"
			echo "     Upstream Coverage = $crtUpstreamCover"
			echo "     Upstream Exon Expression = $crtUpstreamExonCover"
		else
			# A pre-defined setting was used
			echo "     Predefined sort settings: $CALLSETTINGS"
		fi

		# use ANN-based classifier or custom classifier
		if [ $CALLSETTINGS = 'transcriptomeANN' ]
		then
			# Chimeric Filtering use ANN model
			echo " Run ChimANNSort"
			echo "     Rscript chimANNSort.R $libName.pc.lcsv $libName.$RUNID.ann.lion $ANNMODEL"

			$lBIN/Rscript $SCRIPTS/ChimericReadTool/chimAnnSort.R $libName.pc.lcsv $libName.$RUNID.ann.lion $ANNMODEL

		else	
			# Chimeric Filtering using thresholds
			echo "  Run ChimSort"
			echo "     Rscript chimSort.R $libName.pc.lcsv $libName.$RUNID.lion $mappedReads $CRT"
	
			$lBIN/Rscript $SCRIPTS/ChimericReadTool/chimSort.R $libName.pc.lcsv $libName.$RUNID.lion $mappedReads $CRT

			# Sanity Check
			if [ -s $libName.$RUNID.lion ] # was the lions file generated
			then
				echo "Lions file generated successfully."
				lionSuccess='1'
			else
				echo "Lions file not generated. Something is amiss."
				lionSuccess='0'
			fi

		fi

		echo ""

	fi # End Bypass flow control

fi # END Master Flow control

# Summit Log
	# Add counter to summitLog_$RUNID
	# <LibName> <Successful? 1/0> <date>
	echo $libName $lionSuccess $(date) >> $pDIR/summitLog_$RUNID

# Done script :D
