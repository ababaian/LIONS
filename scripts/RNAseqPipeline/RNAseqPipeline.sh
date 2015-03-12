#!/bin/bash
# RNAseqPipeline Script
# Modified for use by Artem Babaian
# :D

set -e

# USAGE -------------------------------------------------------------------------------
# if parameters not ran correctly
usage()
	{
	echo "RNAseqPipeline.sh "
	echo "    Calculate RPKM and coverage data plus TE RPKM data and WIG file"
	echo ""
	echo "USAGE (4 input parameters):"
	echo "     RNAseqPipeline.sh <1: reference_set> <2:Project_name> <3:paiRed/Stranded (R/S)> <4: input_BAM>"
	echo ""
	echo "Prior to running - set-up a reference set (also called species in older versions)"
	echo "     Human genome 19 with GENCODE v14 = hg19gc_v14"
	}

	if [ $# -ne 4 ]; then
		usage
		exit 1
	fi

	# Some example usages (for reference only)
	# sh RNArpkmMaster.sh hg19gc_v14 test_run R HS0988.bam 
	# sh RNArpkmMaster.sh mm9v65 test_run R Lorincz_PGC/bams/RNA-seq.PGC.bam 

# -------------------------------------------------------------------
# INITIALIZAITON ----------------------------------------------------
# -------------------------------------------------------------------

# Script and Binary Directories for RNAseqPipeline

# JAVA access
	export J="$lBIN/java"

# RNAseqPipeline Shell Base
	export SHELL_BASE="$SCRIPTS/RNAseqPipeline"

# Java (Jar) Binaries Base
	export JAVA_BASE="$SHELL_BASE/bin"

# Samtools
	export SAMTOOLS="$lBIN/samtools"

# Dependent Scripts
	rnaMasterScript=$SHELL_BASE/RNAseqMaster.sh
	RegionsCoverageFromWigCalculator=$SHELL_BASE/bin/RegionsCoverageFromWigCalculator_May6_2011.jar


# INPUT PARAMETERS
	species=$1 # Reference Set
	name=$2 # Library Name
	sr=$3 # Stranded/paiRed read type
	bamFile=$4 # BAM-file

# Read length in Bam File
	length=$($SAMTOOLS view $bamFile | awk '($1!~/@/) {s=s+1; if(s==1){print length($10)} else {exit}}')

# Resource for RNAseqPipeline

	# Reference Exon Annotation Folder
	res="$pDIR/$name/resources"

	# Chromosome Sizes File
	chrs="$RESOURCES/genome/$INDEX.chr.size"

	# BWA Conversion File
	chrfile="$RESOURCES/genome/$INDEX.bwa.names"

	echo "=======  Run Parameters  ========" 
	echo " bamFile: $bamFile"
	echo " species: $species"
	echo " name: $name"
	echo " read length: $length bp"
	echo " resource folder: $res" # **********
	echo " chr sizes: $chrs" # ********
	echo " chr name conv.: $chrfile" # *******
	echo ""

#---------------------------------------------------------------------
# EXECUTE SCRIPTS ----------------------------------------------------

# Calculate Coverage and RPKM
	# RNAseqMaster.sh <BAM> <Project_name> <output_folder>
	# <Reference> <S/R> <Quality Threshold> <Running Mask>
	#
	# Where; Running Mask = COVERAGE, RPKM, LEAKAGE,PROFILE,
	# REPORT (1:run, 0:norun)
	
	echo " cmd: "
	echo " sh $SHELL_BASE/RNAseqMaster.sh $bamFile $name `pwd` $species R 0 1,1,0,0,1 $res $chrs $chrfile"
	echo ''
	sh $SHELL_BASE/RNAseqMaster.sh $bamFile $name `pwd` $species $sr 10 1,1,0,0,1 $res $chrs $chrfile

	echo "RNAseqMaster complete"

# Check the files generated are in fact there
	if [ "$sr" == "R" ]; then
		wig=`readlink -f $name/wig/$name.q*.F*.wig.gz`
		mrwig=`readlink -f $name/wig_multiread/$name.q*.F*.wig.gz`
		echo "wig=$wig"
	else
		wig_pos=`readlink -f $name/wig/$name.q*.F*.pos.wig.gz`
		wig_neg=`readlink -f $name/wig/$name.q*.F*.neg.wig.gz`
		mrwig_pos=`readlink -f $name/wig_multiread/$name.q*.F*.pos.wig.gz`
		mrwig_neg=`readlink -f $name/wig_multiread/$name.q*.F*.neg.wig.gz`
		echo "wig_pos=$wig_pos"
		echo "wig_neg=$wig_neg"
	fi


# Calculate RPKM
	# Exonic Reads: protein coding, no MT, no ribo, exclude top 0.005 exons
	exonicreads=`grep excluded $name/coverage/$name.rpkm.log | head -n 1 | cut -f 2`
	rpkm=`echo "(1000000000.0/$length)/$exonicreads" | bc -l`

echo "########## Running TE_stats.sh ##########"

exonicreads=`grep excluded $name/coverage/$name.rpkm.log | head -n 1 | cut -f 2`
rpkm=`echo "(1000000000.0/$length)/$exonicreads" | bc -l`
if [ "$sr" == "R" ]; then
    sh $SHELL_BASE/TE_stats.sh $name $rpkm $species $sr $mrwig
else
    sh $SHELL_BASE/TE_stats.sh $name $rpkm $species $sr $mrwig_pos $mrwig_neg
fi

# run WIG RPKM code

echo "########## Running WIG_RPKM.sh ##########"

if [ "$sr" == "R" ]; then
	rpkm_wig=$(echo `dirname $wig`/rpkm_`basename $wig`)
	sh $SHELL_BASE/WIG_RPKM.sh $wig $rpkm $rpkm_wig
else
	rpkm_wig_pos=$(echo `dirname $wig_pos`/rpkm_`basename $wig_pos`)
	rpkm_wig_neg=$(echo `dirname $wig_neg`/rpkm_`basename $wig_neg`)
	sh $SHELL_BASE/WIG_RPKM.sh $wig $rpkm_pos $rpkm_wig_pos
	sh $SHELL_BASE/WIG_RPKM.sh $wig $rpkm_neg $rpkm_wig_neg
fi

# Script done
# :D
