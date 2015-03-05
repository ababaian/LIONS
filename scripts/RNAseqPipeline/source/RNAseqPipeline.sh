#!/bin/bash
set -e
#set -o pipefail # Commented due to incorrect return values from SAMtools

usage()
{
echo "Script to calculate RPKM and coverage data plus TE RPKM data and WIG file"
echo ""
echo "USAGE (five input parameters):"
echo "/home/jlever/bin/RNAseqPipeline.sh <1: species> <2: name> <3: strand (R/S)> <4/5: input files (1 BAM/1 FASTQ/2 FASTQs)"
echo ""
}

if [ $# -ne 4 -a $# -ne 5 ]; then
	usage
	exit 1
fi

# Some example usages (for reference only)
# sh /home/jlever/bin/RNArpkmMaster.sh hg19 test_run R colon_data/HS0988.bam 
# sh /home/jlever/bin/RNArpkmMaster.sh mm9v65 test_run R Lorincz_PGC/bams/RNA-seq.PGC.bam 
# sh /home/jlever/bin/RNArpkmMaster.sh mm9v65 test_run R Lorincz_PGC/FASTQ/1.fastq Lorincz_PGC/FASTQ/2.fastq

SHELL_BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
rnaMasterScript=$SHELL_BASE/RNAseqMaster.sh
#J=/gsc/software/linux-x86_64/jre-1.6.0_16/bin/java # Use older version of Java for use on Apollo
J=/gsc/software/linux-x86_64/jre1.7.0_03/bin/java
RegionsCoverageFromWigCalculator=/projects/03/genereg/projects/SOLEXA/lib/RegionsCoverageFromWigCalculator_May6_2011.jar
SAMTOOLS=/gsc/software/linux-x86_64/samtools-0.1.13/samtools

species=$1
name=$2
sr=$3

inputA=$4
inputB=$5

if [ $# -eq 4 ]; then
	# Check if BAM or FASTQ
	if [ ${inputA: -4} == ".bam" ]; then
		bamFile=`readlink -f $inputA`
	elif [ ${inputA: -5} == ".fastq" ]; then
		echo "ERROR: Cannot handle single FASTQ file yet" 1>&2
		exit 1
	else
		echo "ERROR: Expected bam or fastq file. Got $inputA" 1>&2
		exit 1
	fi
elif [ $# -eq 5 ]; then
	# Check if both are FASTQ
	if [ ${inputA: -6} == ".fastq" -a ${inputB: -6} == ".fastq" ]; then
		echo "ERROR: TopHatting not implemented yet" 1>&2
		exit 1

		mkdir -p $name
		mkdir -p $name/bam
		outfolder=`readlink -f $name/bam`

		PATH=/home/mkarimi/bin:$PATH

		/projects/03/genereg/projects/mES/tophat/tophat-1.3.1.Linux_x86_64/tophat -r 200 -o $outfolder /projects/mbilenky/resources/bowtie_mm9/mm9 `readlink -f $inputA` `readlink -f $input B`

	else
		echo "ERROR: Both files must be FASTQ files (for TopHatting)" 1>&2
		exit 1
	fi
else
	echo "ERROR: Unexpected number of arguments" 1>&2
	exit 1
fi

length=$($SAMTOOLS view $bamFile | awk '($1!~/@/) {s=s+1; if(s==1){print length($10)} else {exit}}')

echo "bamFile:$bamFile"
echo "species:$species"
echo "name:$name"
echo "length:$length"

# resources
reschrs=$(sh $SHELL_BASE/RNAgetRes.sh $species)
res=$(echo $reschrs | cut -f 1 -d ',')
chrs=$(echo $reschrs | cut -f 2 -d ',')

echo $res
echo $chrs

echo "########## Running RNAseqMaster.sh ##########"

# Execute coverage and RPKM
echo "sh $SHELL_BASE/RNAseqMaster.sh $bamFile $name `pwd` $species R 0 1,1,0,0,0"
#exit 1
sh $SHELL_BASE/RNAseqMaster.sh $bamFile $name `pwd` $species $sr 10 1,1,0,0,0 $J

# check everything is there as needed

echo "RNAseqMaster complete"

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


# extract rpkm from RPKM file
# Total number of exonic reads for RPKM (protein coding; no MT, no ribo proteins, top expressed 0.005 exons excluded)

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


