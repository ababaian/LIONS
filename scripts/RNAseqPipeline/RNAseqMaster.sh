#!/bin/sh
# RNAseqMaster.sh
# Script used to run RNAseq Pipeline
	set -e

# USAGE ---------------------------------------------------------------------
usage()
	{
	echo "Master script to run number of tools and generate an extensive QC report for the strand specific RNA-seq library"
	echo ""
	echo "USAGE (five input parameters):"
	echo "    RNAseqMaster.sh <1: bam file (long path)> <2: name> <3: folder with all output (will be this PATH/name)> <4: species (hg19v66/...)> <5: strand specific(S)/regular(R)> <6: quality threshold> <7: running mask COVERAGE,RPKM,LEAKAGE,PROFILE,REPORT (1==run, 0==don't run))> [<8: java>]"
	echo ""
	echo "OUTPUT: coverage files and coverage distributions; note that for strand specific RNA-seq coverages are calculated for proper strand and then cat together"
	}

if [ $# -lt 7 ]; then
        usage
        exit 1
fi

# INITIALIZE ----------------------------------------------------------
date
echo "Running - $0"
#svn info https://svn01.bcgsc.ca/svn/Solexa_Shell/RNAseqMaster.sh | awk '/Revision/ || /Last Changed Date/'

# Binaries for Files
# JAVA
	if [ -z "$8" ]; then
		J=/gsc/software/linux-x86_64/jre1.7.0_03/bin/java
	else
		J=$8
	fi
	echo "Using java: $J"
	#JAVA_BASE=/projects/03/genereg/projects/SOLEXA/lib/
	JAVA_BASE=/home/ababaian/software/RNAseqPipeline/bin/

#SHELL
	SHELL_BASE='/home/ababaian/software/RNAseqPipeline'
	#SHELL_BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
#SAMTOOLS
	#SAMTOOLS=samtools
	SAMTOOLS=/gsc/software/linux-x86_64/samtools-0.1.13/samtools
	#SAMTOOLS=/home/pubseq/BioSw/samtools/samtools-0.1.16/samtools
	echo "Using samtools: $SAMTOOLS"

# SCRIPT -----------------------------------------------------------------
echo "********** $2 **********";

# input is BAM file directory
	bam=$1 # total path to the input bam file
	name=$2 # name of the library
	dir=$3 # output folder
	species=$4 # species (hg19v66, hg18, etc)
	sr=$5 # strand specific = S, not any other letter, perhaps R (regular)
	#length=$6 # read length (75bp, 100bp,...)
	QC=$6 # quality threshold
		if [ $QC -lt 0 ]; then
			$QC=0
		fi
	
	wtd=$7
		RUN_COVERAGE=$(echo $wtd | cut -f1 -d",");
		RUN_RPKM=$(echo $wtd | cut -f2 -d",");
		RUN_LEAKAGE=$(echo $wtd | cut -f3 -d",");
		RUN_PROFILE=$(echo $wtd | cut -f4 -d",");
		RUN_REPORT=$(echo $wtd | cut -f5 -d",");

# Samtools [-F] flag: reads to ignore in wig
	# 516 =  Read unmapped, read fails QC [Default]
	# 772 =  Read unmapped, read fails QC, not primary alignment
	# 1796 = Read unmapped, read fails QC, not primary alignment, PCR duplicate
	flag='772'

# Directories
	Bdir=$dir"/"$name # folder for individual library
	mkdir -p $Bdir
	cd $Bdir

	if [ -s $name".bam" ]; then
		echo "BAM file is soft-linked"
	else
		ln -s $bam $name".bam" # generating the soft link
	fi

# Calculate Read Length of BAM file
	length=$($SAMTOOLS view $name".bam" | awk '($1!~/@/) {s=s+1; if(s==1){print length($10)} else {exit}}')
	echo "Read length: $length bp"

# defining chr file
	reschrs=$(sh $SHELL_BASE/RNAgetRes.sh $species)
	chrfile=$(echo $reschrs | cut -f 3 -d ',')

# BAM2WIG.jar
# generate wig file from bam
	
	# Make 
	Wdir=$Bdir"/wig"
	mkdir -p $Wdir

	MRWdir=$Bdir"/wig_multiread"
	mkdir -p $MRWdir

# check if file does not exist, then generate

if [ "$sr" == "S" ]; then
	if [ -s $Wdir"/"$name".q"$QC".F$flag.pos.wig.gz" ] && [ -s $Wdir"/"$name".q"$QC".F$flag.neg.wig.gz" ]; then
		echo "Wig files (pos and neg strand) exist"
	else	
		date > $Wdir/$name.wig.log
		echo "$J -jar -Xmx10G $JAVA_BASE/BAM2WIG.jar -bamFile $name".bam" -out $Wdir -q $QC -F $flag -s -samtools $SAMTOOLS -chr $chrfile >> $Wdir/$name.wig.log"
		$J -jar -Xmx10G $JAVA_BASE/BAM2WIG.jar -bamFile $name".bam" -out $MRWdir -q 0 -F $flag -s -samtools $SAMTOOLS -chr $chrfile >> $Wdir/$name.wig.log
		$J -jar -Xmx10G $JAVA_BASE/BAM2WIG.jar -bamFile $name".bam" -out $Wdir -q $QC -F $flag -s -samtools $SAMTOOLS -chr $chrfile >> $Wdir/$name.wig.log
	fi
else
	if [ -s $Wdir"/"$name".q"$QC".F$flag.wig.gz" ]; then
		echo "Wig file exist"
	else	
		date > $Wdir/$name.wig.log
		echo "$J -jar -Xmx10G $JAVA_BASE/BAM2WIG.jar -bamFile $name".bam" -out $Wdir -q $QC -F $flag -samtools $SAMTOOLS -chr $chrfile >> $Wdir/$name.wig.log"
		$J -jar -Xmx10G $JAVA_BASE/BAM2WIG.jar -bamFile $name".bam" -out $MRWdir -q 0 -F $flag -samtools $SAMTOOLS -chr $chrfile >> $Wdir/$name.wig.log
		$J -jar -Xmx10G $JAVA_BASE/BAM2WIG.jar -bamFile $name".bam" -out $Wdir -q $QC -F $flag -samtools $SAMTOOLS -chr $chrfile >> $Wdir/$name.wig.log
	fi
fi	

#echo "Forced stop..."
#exit 1

# Generate Coverage
	Cdir=$Bdir"/coverage";
	mkdir -p $Cdir

	if [ "$RUN_COVERAGE" == "1" ]; then
		echo "Calculating Genomic Coverage..."
		
		echo "$SHELL_BASE/RNAseqCoverageCalculator.sh $name $Wdir $Cdir $species $sr $J > $Cdir/$name.coverage.log"

		$SHELL_BASE/RNAseqCoverageCalculator.sh $name $Wdir $Cdir $species $sr $J > $Cdir/$name.coverage.log
		echo "Coverage Complete"
		#mkdir -p $Cdir"/covDist"
		#mv -f $Cdir"/*.covDist" $Cdir"/covDist/."
	else
		echo "Skipping COVERAGE..."
	fi

# RPKM
cd $Cdir
	if [ "$RUN_RPKM" == "1" ]; then
		echo "RPKM..."
		$SHELL_BASE/RPKM.sh $name $length G A $species $sr
	else
		echo "Skipping RPKM..."
	fi

	cd $Bdir

# Leakage Claculations (if S reads)
if [ "$sr" == "S" ]; then
	Ldir=$Bdir"/leakage/"
	mkdir -p $Ldir
	if [ "$RUN_LEAKAGE" == "1" ]; then
		echo "Leakage..."
		$SHELL_BASE/TCReadsAnalyzer.sh $name".bam" $Ldir $name $species $J > $Ldir/$name.leakage.log
	else
		echo "Skipping LEAKAGE..."
	fi	
else
	echo "Not a strand specific library (no LEAKAGE calculations)"
fi

# cDNA profiles
Pdir=$Bdir"/cDNAProfile"
mkdir -p $Pdir
if [ "$RUN_PROFILE" == "1" ]; then
	echo "cDNA profiles..."
	if [ "$sr" == "S" ]; then
		$SHELL_BASE/cDNAProfileAnalyzer.sh $Wdir"/"$name".*pos.wig*" $Wdir"/"$name".*neg.wig*" $Pdir $name $species $J
	else
		$SHELL_BASE/cDNAProfileAnalyzer.sh $Wdir"/"$name"*wig*" $Wdir"/"$name"*wig*" $Pdir $name $species $J
	fi
else
	echo "Skipping PROFILE..."
fi	

################################
#                              #
# Generating the final report  #
#                              #
################################

if [ "$RUN_REPORT" == "1" ]; then
	$SHELL_BASE/QCreport.sh $name $dir $sr
fi

# Script over
echo " RNAseqMaster.sh Ran to completion "
