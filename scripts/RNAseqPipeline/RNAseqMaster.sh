#!/bin/sh
# RNAseqMaster.sh < Bam_input > <Libray Name> <Output dir> <Resource Name> <S/R> <QC Flag> <Run Profile> <chrfile BWA_conversion file>
# Script used to run RNAseq Pipeline
	set -e

# USAGE ---------------------------------------------------------------------
usage()
	{
	echo "Master script to run number of tools and generate an extensive"
	echo "QC report for the strand specific RNA-seq library"
	echo ""
	echo "USAGE (input parameters):"
	echo " bash RNAseqMaster.sh"
	echo "  <1: bam file (long path)>"
	echo "  <2: name>"
	echo "  <3: folder with all output (will be this PATH/name)>"
	echo "  <4: species (hg19v66/...)>"
	echo "  <5: strand specific(S)/regular(R)>"
	echo "  <6: quality threshold>"
	echo "  <7: running mask COVERAGE,RPKM,LEAKAGE,PROFILE,REPORT (1==run, 0==don't run))>"
	echo "  <8: BWA UCSC Conversion File>"
	echo ""
	echo "OUTPUT: coverage files and coverage distributions; "
	echo " note that for strand specific RNA-seq coverages are calculated"
	echo " for proper strand and then cat together"
	}

# INITIALIZE ----------------------------------------------------------
date
echo "Running - $0"
echo ''



# SCRIPT -----------------------------------------------------------------

# input is BAM file directory
	bam=$1 # input bam file
	name=$2 # name of the library
	dir=$3 # output folder
	species=$4 # species (hg19v66, hg18, etc)
	sr=$5 # strand specific = S, not any other letter, perhaps R (regular)
	#length=$6 # read length (75bp, 100bp,...)
	QC=$6 # quality threshold
		if [ $QC -lt 0 ]; then
			$QC=0
		fi
	
	wtd=$7 # Run Parameters
		RUN_COVERAGE=$(echo $wtd | cut -f1 -d",");
		RUN_RPKM=$(echo $wtd | cut -f2 -d",");
		RUN_LEAKAGE=$(echo $wtd | cut -f3 -d",");
		RUN_PROFILE=$(echo $wtd | cut -f4 -d",");
		RUN_REPORT=$(echo $wtd | cut -f5 -d",");

	# Resourcse
	res=$8
	chrSize=$9
	chrfile=$10

# Ouput 
# Samtools [-F] flag: reads to ignore in wig
	# 516 =  Read unmapped, read fails QC [Default]
	# 772 =  Read unmapped, read fails QC, not primary alignment
	# 1796 = Read unmapped, read fails QC, not primary alignment, PCR duplicate
	#flag='772'
	flag=$(echo $QUALITY | cut -f2 -d'F' - )
	QC=$(echo $QUALITY | cut -f1 -d'.' - | cut -f2 -d'q' - )

# Directory Initialization
	Bdir=$dir"/"expression # folder for individual library
	mkdir -p $Bdir
	cd $Bdir
	ln -fs ../alignment/$bam ./$bam


# Calculate Read Length of BAM file
	length=$($lBIN/samtools view $bam | awk '($1!~/@/) {s=s+1; if(s==1){print length($10)} else {exit}}')

# BAM2WIG.jar
# generate wig file from bam
	# Make 
	Wdir=$Bdir"/wig"
	mkdir -p $Wdir

	MRWdir=$Bdir"/wig_multiread"
	mkdir -p $MRWdir

# check if file does not exist, then generate

if [ "$sr" = "S" ]
then

	if [ -s $Wdir"/"$name".q"$QC".F$flag.pos.wig.gz" ] && [ -s $Wdir"/"$name".q"$QC".F$flag.neg.wig.gz" ]; then
		echo "Wig files (pos and neg strand) exist"
	else	
		date > $Wdir/$name.wig.log
		echo "Generating Wig Files"
		echo "$J -jar -Xmx6G $JAVA_BASE/BAM2WIG.jar -bamFile $bam -out $Wdir -q $QC -F $flag -s -samtools $SAMTOOLS -chr $chrfile >> $Wdir/$name.wig.log"
		echo ''
		$J -jar -Xmx6G $JAVA_BASE/BAM2WIG.jar -bamFile $bam -out $MRWdir -q 0 -F $flag -s -samtools $SAMTOOLS -chr $chrfile >> $Wdir/$name.wig.log
		$J -jar -Xmx6G $JAVA_BASE/BAM2WIG.jar -bamFile $bam -out $Wdir -q $QC -F $flag -s -samtools $SAMTOOLS -chr $chrfile >> $Wdir/$name.wig.log
	fi

else

	if [ -s $Wdir"/"$name".q"$QC".F$flag.wig.gz" ]; then
		echo "Wig file exist"
	else	
		date > $Wdir/$name.wig.log
		echo "$J -jar -Xmx6G $JAVA_BASE/BAM2WIG.jar -bamFile $bam -out $Wdir -q $QC -F $flag -samtools $lBIN/samtools -chr $chrfile >> $Wdir/$name.wig.log"
		echo ''
		$J -jar -Xmx6G $JAVA_BASE/BAM2WIG.jar -bamFile $bam -out $MRWdir -q 0 -F $flag -samtools $SAMTOOLS -chr $chrfile >> $Wdir/$name.wig.log
		$J -jar -Xmx6G $JAVA_BASE/BAM2WIG.jar -bamFile $bam -out $Wdir -q $QC -F $flag -samtools $SAMTOOLS -chr $chrfile >> $Wdir/$name.wig.log
	fi

fi	

echo ''

# Generate Coverage
	Cdir=$Bdir"/coverage";
	mkdir -p $Cdir

	if [ "$RUN_COVERAGE" = "1" ]; then
		echo "Calculating Genomic Coverage..."
		
		echo "$SHELL_BASE/RNAseqCoverageCalculator.sh $name $Wdir $Cdir $species $sr $Cdir/$name.coverage.log $res $chrSize"
		echo ''

		bash $SHELL_BASE/RNAseqCoverageCalculator.sh $name $Wdir $Cdir $species $sr $Cdir/$name.coverage.log $res $chrSize

		echo "Coverage Complete"
		#mkdir -p $Cdir"/covDist"
		#mv -f $Cdir"/*.covDist" $Cdir"/covDist/."
	else
		echo "Skipping COVERAGE..."
	fi

echo ''

# RPKM
#cd $Cdir
#	if [ "$RUN_RPKM" = "1" ]; then
#		echo " RPKM..."
#		$SHELL_BASE/RPKM.sh $name $length G A $species $sr
#	else
#		echo "Skipping RPKM..."
#	fi
#
#	cd $Bdir
#
#echo ''

# Leakage Claculations (if S reads)
if [ "$sr" = "S" ]; then
	Ldir=$Bdir"/leakage/"
	mkdir -p $Ldir
	if [ "$RUN_LEAKAGE" = "1" ]; then
		echo " Leakage..."
		$SHELL_BASE/TCReadsAnalyzer.sh $name".bam" $Ldir $name $species $J > $Ldir/$name.leakage.log
	else
		echo "Skipping LEAKAGE..."
	fi	
else
	echo "Not a strand specific library (no LEAKAGE calculations)"
fi

echo ''

# cDNA profiles
Pdir=$Bdir"/cDNAProfile"
mkdir -p $Pdir
if [ "$RUN_PROFILE" = "1" ]; then
	echo " cDNA profiles..."
	if [ "$sr" = "S" ]; then
		$SHELL_BASE/cDNAProfileAnalyzer.sh $Wdir"/"$name".*pos.wig*" $Wdir"/"$name".*neg.wig*" $Pdir $name $species $J
	else
		$SHELL_BASE/cDNAProfileAnalyzer.sh $Wdir"/"$name"*wig*" $Wdir"/"$name"*wig*" $Pdir $name $species $J
	fi
else
	echo "Skipping PROFILE..."
fi	

echo ''

################################
#                              #
# Generating the final report  #
#                              #
################################

if [ "$RUN_REPORT" = "1" ]; then
	#$SHELL_BASE/QCreport.sh $name $dir $sr
	# not used
fi

# Script over
echo " RNAseqMaster.sh Ran to completion "
