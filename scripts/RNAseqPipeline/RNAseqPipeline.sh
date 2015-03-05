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

# -------------------------------------------------------------------------------------
# INITIALIZAITON ----------------------------------------------------------------------

# Which folder is the script ran out of?
	#SHELL_BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
	SHELL_BASE='/home/ababaian/software/RNAseqPipeline'
# Dependent Scripts
	rnaMasterScript=$SHELL_BASE/RNAseqMaster.sh
	#RegionsCoverageFromWigCalculator=/projects/03/genereg/projects/SOLEXA/lib/RegionsCoverageFromWigCalculator_May6_2011.jar
	RegionsCoverageFromWigCalculator=/home/ababaian/software/RNAseqPipeline/bin/RegionsCoverageFromWigCalculator_May6_2011.jar

# BINARIES
# Java
	#J=java #local install
	#J=/gsc/software/linux-x86_64/jre-1.6.0_16/bin/java # Use older version of Java for use on Apollo
	J=/gsc/software/linux-x86_64/jre1.7.0_03/bin/java

# Samtools
	SAMTOOLS=samtools
	#SAMTOOLS=/gsc/software/linux-x86_64/samtools-0.1.13/samtools

# INPUT PARAMETERS
	species=$1 # Reference Set
	name=$2 # Project Name
	sr=$3 # Stranded/paiRed read type

	inputA=$4 # BAM-file
	inputB=$5 #Obsolete

# Check that the input file is BAM
# 
if [ $# -eq 4 ]; then
	# Check if BAM or FASTQ
	if [ ${inputA: -4} == ".bam" ]; then
		bamFile=`readlink -f $inputA`
	else
		echo "ERROR: Expected bam. Got $inputA" 1>&2
		exit 1
	fi
fi

# Read length in Bam File
	length=$($SAMTOOLS view $bamFile | awk '($1!~/@/) {s=s+1; if(s==1){print length($10)} else {exit}}')

	echo "=======  Run Parameters  ========" 
	echo "bamFile:$bamFile"
	echo "species:$species"
	echo "name:$name"
	echo "read length:$length"
	echo ""

# Resources
	# hg19gc_v14 : hg19 and Gencode v14	
	reschrs=$(sh $SHELL_BASE/RNAgetRes.sh $species)

	# Reference Annotation
	res=$(eval echo $reschrs | cut -f 1 -d ',')

	# Chromosome Sizes
	chrs=$(echo $reschrs | cut -f 2 -d ',')

	echo $res
	echo $chrs

#---------------------------------------------------------------------
# EXECUTE SCRIPTS ----------------------------------------------------
echo ""
echo ""
echo "########## Running RNAseqMaster.sh ##########"

# Calculate Coverage and RPKM
	# RNAseqMaster.sh <BAM> <Project_name> <output_folder>
	# <Reference> <S/R> <Quality Threshold> <Running Mask>
	#
	# Where; Running Mask = COVERAGE, RPKM, LEAKAGE,PROFILE,
	# REPORT (1:run, 0:norun)
	
	echo "sh $SHELL_BASE/RNAseqMaster.sh $bamFile $name `pwd` $species R 0 1,1,0,0,1"
	sh $SHELL_BASE/RNAseqMaster.sh $bamFile $name `pwd` $species $sr 10 1,1,0,0,1 $J

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
