#!/bin/bash
# RNAseqCoverageCalculator.sh
# Modified for LIONS pipeline

set -e

usage()
{
echo "Script to calculate 'exonic' coverage for strand-specific RNA-seq"
echo ""
echo "USAGE (five input parameters):"
echo "/home/jlever/bin/RNAseqCoverageCalculator.sh <1: library name> <2: PATH to IN directory with the wig files> <3: PATH to OUT directory where coverage files will be stored> <4: species> <5: strand specific(S)/regular(R)> [<6: java>]"
echo ""
echo "OUTPUT: coverage files and coverage distributions; note that for strand specific RNA-seq coverages are calculated for proper strand and then cat together"
}


date
echo "Running RNAseqCoverageCalculator - $0 ------------"

# Input Parameters
# RNAseqCoverageCalculator.sh $name $Wdir $Cdir $species $sr $Cdir/$name.coverage.log $res $chrSize

	echo "Run Parameters"
	echo "     Name: "$1
	echo "     Wigd: "$2
	echo "     output: "$3 # Directory or file?
	echo "     Reference Set: "$4
	echo "     S/R: "$5
	echo "     res: "$7
	echo "     chrSize: "$8
	echo ""
	echo ""

	name=$1
	wigd=$2
	out=$3
	species=$4
	sr=$5	# strand specific or regular S/R
	res=$7
	chrSize=$8

	mkdir -p $out

# SCRIPT ------------------------------------------------------------


if [ "$sr" == "S" ]; then
# Strand specific RNA-seq
	prefs=(pos neg)
	for pref in ${prefs[*]}; do
		
		fwig=$wigd"/"$name.q*.F*.$pref.wig.gz
		
		if [ -s $fwig ]; then
			echo "$fwig exists"
			echo ""
		else
			echo "$fwig does not exists. Exiting..."
			exit 1
		fi
			
		$J -jar -Xmx5G $JAVA_BASE/RegionsCoverageFromWigCalculator.jar -w $fwig -r $res"/"$species"_genes_for_introns_"$pref -s $chrSize -o $out -n $name
		$J -jar -Xmx5G $JAVA_BASE/RegionsCoverageFromWigCalculator.jar -w $fwig -r $res"/"$species"_exons_for_genes_"$pref -s $chrSize -o $out -n $name
		$J -jar -Xmx5G $JAVA_BASE/RegionsCoverageFromWigCalculator.jar -w $fwig -r $res"/"$species"_genes_for_coverage_"$pref -s $chrSize -o $out -n $name
		$J -jar -Xmx5G $JAVA_BASE/RegionsCoverageFromWigCalculator.jar -w $fwig -r $res"/"$species"_exons_for_coverage_"$pref -s $chrSize -o $out -n $name
	
	done

	# combining all together
	cd $out
	cat $species"_genes_for_introns_pos."$name".coverage" $species"_genes_for_introns_neg."$name".coverage" > $species"_genes_for_introns."$name".coverage"
	cat $species"_exons_for_genes_pos."$name".coverage" $species"_exons_for_genes_neg."$name".coverage" > $species"_exons_for_genes."$name".coverage"
	cat $species"_genes_for_coverage_pos."$name".coverage" $species"_genes_for_coverage_neg."$name".coverage" > $species"_genes_for_coverage."$name".coverage"		
	cat $species"_exons_for_coverage_pos."$name".coverage" $species"_exons_for_coverage_neg."$name".coverage" > $species"_exons_for_coverage."$name".coverage"
	rm -rf *pos* *neg*


else #SR = R
# Regular RNAseq ----------------------------------------------------

	fwig=$wigd"/"$name.q*.F*.wig.gz

	if [ -s $fwig ]; then
		echo " Check that wig file was generated"
		echo " $fwig exists"
		echo ""
	else
		echo " $fwig does not exists. Exiting..."
		exit 1
	fi
	
	echo "Running RegionsCoverageFromWigCalculator"
	echo ""

	echo "$J -jar -Xmx5G $JAVA_BASE/RegionsCoverageFromWigCalculator.jar -w $fwig -r $res"/"$species"_genes_for_introns" -s $chrSize -o $out -n $name"
	$J -jar -Xmx5G $JAVA_BASE/RegionsCoverageFromWigCalculator.jar -w $fwig -r $res"/"$species"_genes_for_introns" -s $chrSize -o $out -n $name
	echo "$J -jar -Xmx5G $JAVA_BASE/RegionsCoverageFromWigCalculator.jar -w $fwig -r $res"/"$species"_exons_for_genes" -s $chrSize -o $out -n $name"
	$J -jar -Xmx5G $JAVA_BASE/RegionsCoverageFromWigCalculator.jar -w $fwig -r $res"/"$species"_exons_for_genes" -s $chrSize -o $out -n $name
	echo "$J -jar -Xmx5G $JAVA_BASE/RegionsCoverageFromWigCalculator.jar -w $fwig -r $res"/"$species"_genes_for_coverage" -s $chrSize -o $out -n $name"
	$J -jar -Xmx5G $JAVA_BASE/RegionsCoverageFromWigCalculator.jar -w $fwig -r $res"/"$species"_genes_for_coverage" -s $chrSize -o $out -n $name
	echo "$J -jar -Xmx5G $JAVA_BASE/RegionsCoverageFromWigCalculator.jar -w $fwig -r $res"/"$species"_exons_for_coverage" -s $chrSize -o $out -n $name"
	$J -jar -Xmx5G $JAVA_BASE/RegionsCoverageFromWigCalculator.jar -w $fwig -r $res"/"$species"_exons_for_coverage" -s $chrSize -o $out -n $name

fi

#covDist=$out"/covDist"
#mkdir -p $covDist
#mv *.covDist $covDist/.

