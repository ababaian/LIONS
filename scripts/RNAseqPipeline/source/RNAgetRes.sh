#!/bin/sh
set -e
set -o pipefail


	#usage

species=$1


if [ "$species" == "mm9v61" ]; then
	res=/projects/03/genereg/projects/SOLEXA/resources/Ensembl_v61/mm9v61/
	chrs=/projects/03/genereg/projects/SOLEXA/chr_info/mm9.chrom.sizes
	chrfile=/projects/03/genereg/projects/SOLEXA/resources/UCSC_chr/mm9.bwa2ucsc.names
	bwtindex=/projects/mbilenky/resources/bowtie_mm9/
	TEdata=/projects/mbilenky/mlc/RepeatMasker/RetroTransposons
	chimeric=/projects/mbilenky/resources/RepeatMasker/ForChimericSearch_mm9
elif [ "$species" == "mm9v65" ]; then
	res=/projects/03/genereg/projects/SOLEXA/resources/Ensembl_v65/mm9v65/
	chrs=/projects/03/genereg/projects/SOLEXA/chr_info/mm9.chrom.sizes
	chrfile=/projects/03/genereg/projects/SOLEXA/resources/UCSC_chr/mm9.bwa2ucsc.names
	bwtindex=/projects/mbilenky/resources/bowtie_mm9/
	TEdata=/projects/mbilenky/mlc/RepeatMasker/RetroTransposons
	chimeric=/projects/mbilenky/resources/RepeatMasker/ForChimericSearch_mm9
elif [ "$species" == "mm9refSeq" ]; then
	res=/projects/03/genereg/projects/SOLEXA/resources/UCSC_mm9/refGene_20120726/
	chrs=/projects/03/genereg/projects/SOLEXA/chr_info/mm9.chrom.sizes
	chrfile=/projects/03/genereg/projects/SOLEXA/resources/UCSC_chr/mm9.bwa2ucsc.names
	bwtindex=/projects/mbilenky/resources/bowtie_mm9/
	TEdata=/projects/mbilenky/mlc/RepeatMasker/RetroTransposons
	chimeric=/projects/mbilenky/resources/RepeatMasker/ForChimericSearch_mm9
elif [ "$species" == "mm9v67" ]; then
 	res=/projects/03/genereg/projects/SOLEXA/resources/Ensembl_v67/mm9v67/
	chrs=/projects/03/genereg/projects/SOLEXA/chr_info/mm9.chrom.sizes
	chrfile=/projects/03/genereg/projects/SOLEXA/resources/UCSC_chr/mm9.bwa2ucsc.names
	bwtindex=/projects/mbilenky/resources/bowtie_mm9/
	TEdata=/projects/mbilenky/mlc/RepeatMasker/RetroTransposons
	chimeric=/projects/mbilenky/resources/RepeatMasker/ForChimericSearch_mm9
	
elif [ "$species" == "hg18" ]; then
 	res=/projects/03/genereg/projects/SOLEXA/resources/Ensembl_v54/hg18/
	chrs=/projects/03/genereg/projects/SOLEXA/chr_info/hg18.chrom.sizes
	chrfile=/projects/03/genereg/projects/SOLEXA/resources/UCSC_chr/hg18.bwa2ucsc.names
	bwtindex=/projects/mbilenky/resources/bowtie_hg18/
	TEdata=/projects/mbilenky/resources/RepeatMasker/SINES_LINES_LTRS_hg18
	chimeric=/projects/mbilenky/resources/RepeatMasker/ForChimericSearch_hg18
elif [ "$species" == "hg18v54" ] ; then
	res=/projects/03/genereg/projects/SOLEXA/resources/Ensembl_v54/hg18v54/
	chrs=/projects/03/genereg/projects/SOLEXA/chr_info/hg18.chrom.sizes	
	chrfile=/projects/03/genereg/projects/SOLEXA/resources/UCSC_chr/hg18.bwa2ucsc.names
	bwtindex=/projects/mbilenky/resources/bowtie_hg18/
	TEdata=/projects/mbilenky/resources/RepeatMasker/SINES_LINES_LTRS_hg18
	chimeric=/projects/mbilenky/resources/RepeatMasker/ForChimericSearch_hg18
	
elif [ "$species" == "hg19v59" ]; then
 	res=/projects/03/genereg/projects/SOLEXA/resources/Ensembl_v59/hg19v59/
	chrs=/projects/03/genereg/projects/SOLEXA/chr_info/hg19.chrom.sizes
	chrfile=/projects/03/genereg/projects/SOLEXA/resources/UCSC_chr/hg19.bwa2ucsc.names
	bwtindex=/projects/mbilenky/resources/bowtie_hg19/
	TEdata=/projects/mbilenky/resources/RepeatMasker/SINES_LINES_LTRS_hg19
	chimeric=/projects/mbilenky/resources/RepeatMasker/ForChimericSearch_hg19
elif [ "$species" == "hg19v61" ]; then
 	res=/projects/03/genereg/projects/SOLEXA/resources/Ensembl_v61/hg19v61/
	chrs=/projects/03/genereg/projects/SOLEXA/chr_info/hg19.chrom.sizes
	chrfile=/projects/03/genereg/projects/SOLEXA/resources/UCSC_chr/hg19.bwa2ucsc.names
	bwtindex=/projects/mbilenky/resources/bowtie_hg19/
	TEdata=/projects/mbilenky/resources/RepeatMasker/SINES_LINES_LTRS_hg19
	chimeric=/projects/mbilenky/resources/RepeatMasker/ForChimericSearch_hg19
elif [ "$species" == "hg19v66" ]; then
 	res=/projects/03/genereg/projects/SOLEXA/resources/Ensembl_v66/hg19v66/
	chrs=/projects/03/genereg/projects/SOLEXA/chr_info/hg19.chrom.sizes
	chrfile=/projects/03/genereg/projects/SOLEXA/resources/UCSC_chr/hg19.bwa2ucsc.names
	bwtindex=/projects/mbilenky/resources/bowtie_hg19/
	TEdata=/projects/mbilenky/resources/RepeatMasker/SINES_LINES_LTRS_hg19
	chimeric=/projects/mbilenky/resources/RepeatMasker/ForChimericSearch_hg19
elif [ "$species" == "hg19v65" ]; then
 	res=/projects/03/genereg/projects/SOLEXA/resources/Ensembl_v65/hg19v65/
	chrs=/projects/03/genereg/projects/SOLEXA/chr_info/hg19.chrom.sizes
	chrfile=/projects/03/genereg/projects/SOLEXA/resources/UCSC_chr/hg19.bwa2ucsc.names
	bwtindex=/projects/mbilenky/resources/bowtie_hg19/
	TEdata=/projects/mbilenky/resources/RepeatMasker/SINES_LINES_LTRS_hg19
	chimeric=/projects/mbilenky/resources/RepeatMasker/ForChimericSearch_hg19
elif [ "$species" == "hg19v69" ]; then
 	res=/projects/03/genereg/projects/SOLEXA/resources/Ensembl_v69/hg19v69/
	chrs=/projects/03/genereg/projects/SOLEXA/chr_info/hg19.chrom.sizes
	chrfile=/projects/03/genereg/projects/SOLEXA/resources/UCSC_chr/hg19.bwa2ucsc.names
	bwtindex=/projects/mbilenky/resources/bowtie_hg19/
	TEdata=/projects/mbilenky/resources/RepeatMasker/SINES_LINES_LTRS_hg19
	chimeric=/projects/mbilenky/resources/RepeatMasker/ForChimericSearch_hg19
	
elif [ "$species" == "mm9rs_0313" ]; then
	res=/projects/mbilenky/resources/mm9rs_0313/
	chrs=/projects/03/genereg/projects/SOLEXA/chr_info/mm9.chrom.sizes
	chrfile=/projects/03/genereg/projects/SOLEXA/resources/UCSC_chr/mm9.bwa2ucsc.names
	bwtindex=/projects/mbilenky/resources/bowtie_mm9/
	TEdata=/projects/mbilenky/mlc/RepeatMasker/RetroTransposons
	chimeric=/projects/mbilenky/resources/RepeatMasker/ForChimericSearch_mm9
elif [ "$species" == "hg18rs_0313" ]; then
	res=/projects/mbilenky/resources/hg18rs_0313/
	chrs=/projects/03/genereg/projects/SOLEXA/chr_info/hg18.chrom.sizes
	chrfile=/projects/03/genereg/projects/SOLEXA/resources/UCSC_chr/hg18.bwa2ucsc.names
	bwtindex=/projects/mbilenky/resources/bowtie_hg18/
	TEdata=/projects/mbilenky/resources/RepeatMasker/SINES_LINES_LTRS_hg18
	chimeric=/projects/mbilenky/resources/RepeatMasker/ForChimericSearch_hg18
elif [ "$species" == "hg19rs_0313" ]; then
	res=/projects/mbilenky/resources/hg19rs_0313/
	chrs=/projects/03/genereg/projects/SOLEXA/chr_info/hg19.chrom.sizes
	chrfile=/projects/03/genereg/projects/SOLEXA/resources/UCSC_chr/hg19.bwa2ucsc.names
	bwtindex=/projects/mbilenky/resources/bowtie_hg19/
	TEdata=/projects/mbilenky/resources/RepeatMasker/SINES_LINES_LTRS_hg19
	chimeric=/projects/mbilenky/resources/RepeatMasker/ForChimericSearch_hg19
    
elif [ "$species" == "hg18gc_v3" ]; then
	res=/projects/mbilenky/resources/hg18gc_v3/
	chrs=/projects/03/genereg/projects/SOLEXA/chr_info/hg18.chrom.sizes
	chrfile=/projects/03/genereg/projects/SOLEXA/resources/UCSC_chr/hg18.bwa2ucsc.names
	bwtindex=/projects/mbilenky/resources/bowtie_hg18/
	TEdata=/projects/mbilenky/resources/RepeatMasker/SINES_LINES_LTRS_hg18
	chimeric=/projects/mbilenky/resources/RepeatMasker/ForChimericSearch_hg18
elif [ "$species" == "hg19gc_v14" ]; then
	res=/projects/mbilenky/resources/hg19gc_v14/
	chrs=/projects/03/genereg/projects/SOLEXA/chr_info/hg19.chrom.sizes
	chrfile=/projects/03/genereg/projects/SOLEXA/resources/UCSC_chr/hg19.bwa2ucsc.names
	bwtindex=/projects/mbilenky/resources/bowtie_hg19/
	TEdata=/projects/mbilenky/resources/RepeatMasker/SINES_LINES_LTRS_hg19
	chimeric=/projects/mbilenky/resources/RepeatMasker/ForChimericSearch_hg19
	
else
	echo "RNAgetRes.sh ERROR: No matching resources for species: $species" 1>&2
	exit 1
fi

echo $res,$chrs,$chrfile,$TEdata,$chimeric
