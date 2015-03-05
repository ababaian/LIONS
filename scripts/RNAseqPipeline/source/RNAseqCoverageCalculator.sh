#!/bin/sh
set -e
set -o pipefail

usage()
{
echo "Script to calculate 'exonic' coverage for strand-specific RNA-seq"
echo ""
echo "USAGE (five input parameters):"
echo "/home/jlever/bin/RNAseqCoverageCalculator.sh <1: library name> <2: PATH to IN directory with the wig files> <3: PATH to OUT directory where coverage files will be stored> <4: species> <5: strand specific(S)/regular(R)> [<6: java>]"
echo ""
echo "OUTPUT: coverage files and coverage distributions; note that for strand specific RNA-seq coverages are calculated for proper strand and then cat together"
}

if [ $# -lt 5 ]; then
        usage
        exit 1
fi

date
echo "Running - $0"
svn info https://svn01.bcgsc.ca/svn/Solexa_Shell/RNAseqCoverageCalculator.sh | awk '/Revision/ || /Last Changed Date/'

name=$1
wigd=$2
out=$3
species=$4
# strand specific or regular S/R
sr=$5

if [ -z "$6" ]; then
	J=/gsc/software/linux-x86_64/jre1.7.0_03/bin/java
else
	J=$6
fi

echo "Using java: $J"
JAVA_BASE=/projects/03/genereg/projects/SOLEXA/lib/
SHELL_BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


# resources
reschrs=$(sh $SHELL_BASE/RNAgetRes.sh $species)
res=$(echo $reschrs | cut -f 1 -d ',')
chrs=$(echo $reschrs | cut -f 2 -d ',')

mkdir -p $out

if [ "$sr" == "S" ]; then
	# strand specific RNA-seq
	prefs=(pos neg)
	for pref in ${prefs[*]}; do
		
		fwig=$wigd"/"$name.q*.F*.$pref.wig.gz
		
		if [ -s $fwig ]; then
			echo "$fwig exists"
		else
			echo "$fwig does not exists. Exiting..."
			exit 1
		fi
			
		$J -jar -Xmx5G $JAVA_BASE/RegionsCoverageFromWigCalculator.jar -w $fwig -r $res"/"$species"_genes_for_introns_"$pref -s $chrs -o $out -n $name
		$J -jar -Xmx5G $JAVA_BASE/RegionsCoverageFromWigCalculator.jar -w $fwig -r $res"/"$species"_exons_for_genes_"$pref -s $chrs -o $out -n $name
		$J -jar -Xmx5G $JAVA_BASE/RegionsCoverageFromWigCalculator.jar -w $fwig -r $res"/"$species"_genes_for_coverage_"$pref -s $chrs -o $out -n $name
		$J -jar -Xmx5G $JAVA_BASE/RegionsCoverageFromWigCalculator.jar -w $fwig -r $res"/"$species"_exons_for_coverage_"$pref -s $chrs -o $out -n $name
	
	done
	# combining all together
	cd $out
	cat $species"_genes_for_introns_pos."$name".coverage" $species"_genes_for_introns_neg."$name".coverage" > $species"_genes_for_introns."$name".coverage"
	cat $species"_exons_for_genes_pos."$name".coverage" $species"_exons_for_genes_neg."$name".coverage" > $species"_exons_for_genes."$name".coverage"
	cat $species"_genes_for_coverage_pos."$name".coverage" $species"_genes_for_coverage_neg."$name".coverage" > $species"_genes_for_coverage."$name".coverage"		
	cat $species"_exons_for_coverage_pos."$name".coverage" $species"_exons_for_coverage_neg."$name".coverage" > $species"_exons_for_coverage."$name".coverage"
	rm -rf *pos* *neg*
else
	# Regular RNA-seq
	fwig=$wigd"/"$name.q*.F*.wig.gz

	if [ -s $fwig ]; then
		echo "$fwig exists"
	else
		echo "$fwig does not exists. Exiting..."
		exit 1
	fi
	
	
	echo "$J -jar -Xmx5G $JAVA_BASE/RegionsCoverageFromWigCalculator.jar -w $fwig -r $res\"/\"$species\"_genes_for_introns\" -s $chrs -o $out -n $name"
	$J -jar -Xmx5G $JAVA_BASE/RegionsCoverageFromWigCalculator.jar -w $fwig -r $res"/"$species"_genes_for_introns" -s $chrs -o $out -n $name
	echo "$J -jar -Xmx5G $JAVA_BASE/RegionsCoverageFromWigCalculator.jar -w $fwig -r $res\"/\"$species\"_exons_for_genes\" -s $chrs -o $out -n $name"
	$J -jar -Xmx5G $JAVA_BASE/RegionsCoverageFromWigCalculator.jar -w $fwig -r $res"/"$species"_exons_for_genes" -s $chrs -o $out -n $name
	echo "$J -jar -Xmx5G $JAVA_BASE/RegionsCoverageFromWigCalculator.jar -w $fwig -r $res\"/\"$species\"_genes_for_coverage\" -s $chrs -o $out -n $name"
	$J -jar -Xmx5G $JAVA_BASE/RegionsCoverageFromWigCalculator.jar -w $fwig -r $res"/"$species"_genes_for_coverage" -s $chrs -o $out -n $name
	echo "$J -jar -Xmx5G $JAVA_BASE/RegionsCoverageFromWigCalculator.jar -w $fwig -r $res\"/\"$species\"_exons_for_coverage\" -s $chrs -o $out -n $name"
	$J -jar -Xmx5G $JAVA_BASE/RegionsCoverageFromWigCalculator.jar -w $fwig -r $res"/"$species"_exons_for_coverage" -s $chrs -o $out -n $name

fi

#covDist=$out"/covDist"
#mkdir -p $covDist
#mv *.covDist $covDist/.



