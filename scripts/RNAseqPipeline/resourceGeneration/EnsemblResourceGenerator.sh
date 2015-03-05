#!/bin/sh
set -e

# Resources
# Download from biomart
#----------------------
# exon file: geneId transcriptId chr exon_start exon_end strand rank_in_transcript gene_biotype
# gene file: geneId chr gene_start gene_end strand gene_biotype MGI_symbol/HGNC_symbol description
# transcript file: geneId transcriptId chr transcript_start transcript_end strand gene_biotype transcript_biotype MGI_symbol/HGNC_symbol description
#
date
echo "Running - $0"
#svn info https://svn01.bcgsc.ca/svn/Solexa_Shell/EnsemblResourceGenerator.sh | awk '/Revision/ || /Last Changed Date/'

usage()
{
echo "Script to generate resources for coverage calculations EnsemblResourceGenerator.sh v1.0 (27/03/2012)"
echo ""
echo "USAGE:"
echo "/home/mbilenky/bin/resourceGenerator.sh <prefix [mm9v61/hg18v54/hg19v61/hg19v66..]>"
echo ""
}

if [ "$#" -ne 1 ] ; then
        usage
        exit 1
fi

prefix=$1

e=$prefix"_exons";

if [ -s $e ]; then
	echo "$e exists"
else
	echo "$e does not exists. You have to be in a folder with this file! Download the file from Ensembl. Exiting..."
	exit 1
fi

# test number of fields in the file and consistency

less $e | awk '!/Biotype/ && ($3!~/_/) {if($3=="MT"){$3="M"}; print $0}' | sed 's/ /\t/g' | sort -k1,1 > x; mv -f x $e

eft=$prefix"_exons_for_transcripts"
less $e  | awk '{print $3"\t"$4"\t"$5"\t"$1"_"$2"\t"$6}' | sort -k4,4 > $eft

eftp=$prefix"_exons_for_transcripts_pos";
eftn=$prefix"_exons_for_transcripts_neg";

efg=$prefix"_exons_for_genes";
efgp=$prefix"_exons_for_genes_pos";
efgn=$prefix"_exons_for_genes_neg";

efc=$prefix"_exons_for_coverage";
efcp=$prefix"_exons_for_coverage_pos";
efcn=$prefix"_exons_for_coverage_neg";

eb=$prefix"_exons_boundaries";
ebp=$prefix"_exons_boundaries_pos";
ebn=$prefix"_exons_boundaries_neg";
ebt=$prefix"_exons_boundaries_t";

multi=$prefix"_multi_exon_transcripts";
single=$prefix"_single_exon_transcripts";

echo "$efg - is being generated..."
less $eft | sed 's/_/\t/' | cut -f1-4,6 | sort -k4,4 -k2,2n | awk '{if(id==$4 && $2<=end) {if($3>end){end=$3}} else {if(id!=null){print chr"\t"start"\t"end"\t"id"\t"strand}; id=$4;chr=$1;start=$2;end=$3;strand=$5}} END{print chr"\t"start"\t"end"\t"id"\t"strand}' | sort -k4,4 | uniq  > $efg

echo "Transcripts ..."
less $eft | awk '{if($5==1) print $0}' > $eftp
less $eft | awk '{if($5==-1) print $0}' > $eftn

echo "Genes ..."
less $efg | awk '{if($5==1) print $0}' > $efgp
less $efg | awk '{if($5==-1) print $0}' > $efgn
less $efg | awk '{if($5==1) print $0}' | sort -k4,4 -k2,2n | awk '{s=s+1; print $1"\t"$2"\t"$3"\t"$5"\t"$4"_"s}' > $efgp".profile"
less $efg | awk '{if($5==-1) print $0}'| sort -k4,4 -k2,2nr | awk '{s=s+1; print $1"\t"$2"\t"$3"\t"$5"\t"$4"_"s}' > $efgn".profile"

echo "Coverage ..."
less $efg | cut -f1-3 | sort -k1,1 -k2,2n | uniq | awk '{if(first==0) {chr=$1;start=$2;end=$3;first=1}; if(chr==$1 && $2<=end) {if($3>end) {end=$3}} else {print chr"\t"start"\t"end; chr=$1; start=$2; end=$3}} END{print chr"\t"start"\t"end}' > $efc
less $efg | awk '{if($5==1) print $0}' | cut -f1-3 | sort -k1,1 -k2,2n | uniq | awk '{if(first==0) {chr=$1;start=$2;end=$3;first=1}; if(chr==$1 && $2<=end) {if($3>end) {end=$3}} else {print chr"\t"start"\t"end; chr=$1; start=$2; end=$3}} END{print chr"\t"start"\t"end}' > $efcp
less $efg | awk '{if($5==-1) print $0}' | cut -f1-3 | sort -k1,1 -k2,2n | uniq | awk '{if(first==0) {chr=$1;start=$2;end=$3;first=1}; if(chr==$1 && $2<=end) {if($3>end) {end=$3}} else {print chr"\t"start"\t"end; chr=$1; start=$2; end=$3}} END{print chr"\t"start"\t"end}' > $efcn

echo "Single/multi exonic transcripts..."
less $eft | cut -f4 | sort | uniq -c | awk '{if($1>1) print $2}' | sort > $multi
less $eft | cut -f4 | sort | uniq -c | awk '{if($1==1) print $2}' | sort > $single

echo "Boundaries..."
join -1 4 -2 1 $eft $multi | sed 's/ /\t/g' | sort -k1,1 -k3,3n | awk '{print $1"\t"$2"\t"$3"\t"$3"\t"$5"\n"$1"\t"$2"\t"$4"\t"$4"\t"$5}' | sort -k1,1 -k3,3 | awk '{if($1==id) {if(k>0){print out}; out=$0;k=k+1;} else {id=$1; k=0; }}' | awk '{print $2"\t"$3"\t"$4"\t"$1"\t"$5}' > $ebt
less $ebt | sed 's/_ENS/\tENS/' | cut -f1-4,6 | sort -k1,1 -k2,2n | uniq > $eb
less $eb | awk '{if($5==1) print $0}' > $ebp
less $eb | awk '{if($5==-1) print $0}' > $ebn


g=$prefix"_genes";

if [ -s $g ]; then
	echo "$g exists"
else
	echo "$g does not exists. You have to be in a folder with this file! Download the file from Ensembl. Exiting..."
	exit 1
fi

less $g | sed 's/ /_/g' | awk '!/Biotype/ && ($2!~/_/) {if($2=="MT"){$2="M"}; print $0}' | sed 's/ /\t/g' | awk '{if(NF==8) {print $0} else if(NF==7) {out=$1;for(i=2;i<=6;i++){out=out"\t"$i};out=out"\tNA\t"$7; print out} else {print $0"\tNA\tNA"}}' > x; mv -f x $g;

echo "Genes for introns..."
less $g | cut -f1-5 | awk '{print $2"\t"$3"\t"$4"\t"$1"\t"$5}' | sort -k1,1 | uniq > $g"_for_introns"
less $g | cut -f1-5 | awk '{if($5==1){print $2"\t"$3"\t"$4"\t"$1"\t"$5}}' | sort -k1,1 | uniq > $g"_for_introns_pos"
less $g | cut -f1-5 | awk '{if($5==-1){print $2"\t"$3"\t"$4"\t"$1"\t"$5}}' | sort -k1,1 | uniq > $g"_for_introns_neg"

echo "Genes for coverage..."
less $g | cut -f2-4 | sort -k1,1 -k2,2n | uniq | awk '{if(first==0) {chr=$1;start=$2;end=$3;first=1}; if(chr==$1 && $2<=end) {if($3>end) {end=$3}} else {print chr"\t"start"\t"end; chr=$1; start=$2; end=$3}} END{print chr"\t"start"\t"end}' > $g"_for_coverage"
less $g | awk '{if($5==1) print $0}' | cut -f2-4 | sort -k1,1 -k2,2n | uniq | awk '{if(first==0) {chr=$1;start=$2;end=$3;first=1}; if(chr==$1 && $2<=end) {if($3>end) {end=$3}} else {print chr"\t"start"\t"end; chr=$1; start=$2; end=$3}} END{print chr"\t"start"\t"end}' > $g"_for_coverage_pos"
less $g | awk '{if($5==-1) print $0}' | cut -f2-4 | sort -k1,1 -k2,2n | uniq | awk '{if(first==0) {chr=$1;start=$2;end=$3;first=1}; if(chr==$1 && $2<=end) {if($3>end) {end=$3}} else {print chr"\t"start"\t"end; chr=$1; start=$2; end=$3}} END{print chr"\t"start"\t"end}' > $g"_for_coverage_neg"

echo "Gene types.."
less $g | awk '/ribosomal_protein/' | cut -f1 | sort -k1,1 | uniq > $prefix"_genes.rb.EnsID"
less $g | awk '($2=="M")' | cut -f1 | sort -k1,1 | uniq > $prefix"_genes.mt.EnsID"
less $g | awk '($2!="M") && !/ribosomal_protein/ && ($6=="protein_coding")' | cut -f1 | sort -k1,1 | uniq > $prefix"_genes.pc.EnsID"
less $g | awk '($2!="M") && !/ribosomal_protein/ && ($6!="protein_coding")' | cut -f1 | sort -k1,1 | uniq > $prefix"_genes.nc.EnsID"

###############
# Exon files for different types
types=(mt nc pc rb);
for type in ${types[*]}; do
echo "$type"
join $e $prefix"_genes."$type".EnsID" | sed 's/ /\t/g' | cut -f3-5 | sort -k1,1 -k2,2n | uniq | awk '{if(first==0) {chr=$1;start=$2;end=$3;first=1}; if(chr==$1 && $2<=end) {if($3>end) {end=$3}} else {print chr"\t"start"\t"end; chr=$1; start=$2; end=$3}} END{print chr"\t"start"\t"end}' | sed 's/\t/:/g' | sort -k1,1 >  $prefix"_exons_for_coverage."$type
join $e $prefix"_genes."$type".EnsID" | sed 's/ /\t/g' | awk '{if($6==1) print $0}' | cut -f3-5 | sort -k1,1 -k2,2n | uniq | awk '{if(first==0) {chr=$1;start=$2;end=$3;first=1}; if(chr==$1 && $2<=end) {if($3>end) {end=$3}} else {print chr"\t"start"\t"end; chr=$1; start=$2; end=$3}} END{print chr"\t"start"\t"end}' | sed 's/\t/:/g' > tmp
join $e $prefix"_genes."$type".EnsID" | sed 's/ /\t/g' | awk '{if($6==-1) print $0}' | cut -f3-5 | sort -k1,1 -k2,2n | uniq | awk '{if(first==0) {chr=$1;start=$2;end=$3;first=1}; if(chr==$1 && $2<=end) {if($3>end) {end=$3}} else {print chr"\t"start"\t"end; chr=$1; start=$2; end=$3}} END{print chr"\t"start"\t"end}' | sed 's/\t/:/g' >> tmp
sort -k1,1 tmp > $prefix"_exons_for_coverage.ss."$type; rm -rf tmp
done
###############

# for 5-3 bias we need to do this at the transcript level for protein coding genes only; will keep track of number of exons in the gene

t=$prefix"_transcripts";

if [ -s $t ]; then
	echo "$t exists"
else
	echo "$t does not exists. You have to be in a folder with this file! Download the file from Ensembl. Exiting..."
	exit 1
fi

less $t | sed 's/ /_/g' | awk '!/Biotype/ && ($3!~/_/) {if($3=="MT"){$3="M"}; print $0}' | sed 's/ /\t/g' | awk '{if(NF==10) {print $0} else if(NF==9) {out=$1;for(i=2;i<=8;i++){out=out"\t"$i};out=out"\tNA\t"$9; print out} else {print $0"\tNA\tNA"}}' > x; mv -f x $t;

less $t | awk '{if($7=="protein_coding" && $8=="protein_coding") print $0}' | cut -f2 | sort -k1,1 | uniq > $t.pc.pc.EnstID
less $e | sort -k2,2 > x

join -1 1 -2 2 $t.pc.pc.EnstID x | awk '{print $3"\t"$4"\t"$5"\t"$6"\t"$1"_"$2"_"$7}' > $prefix"_exons_for_transcripts.pc.pc.profile"
join -1 1 -2 2 $t.pc.pc.EnstID x | awk '($6==1) {print $3"\t"$4"\t"$5"\t"$6"\t"$1"_"$2"_"$7}' > $prefix"_exons_for_transcripts_pos.pc.pc.profile"
join -1 1 -2 2 $t.pc.pc.EnstID x | awk '($6==-1) {print $3"\t"$4"\t"$5"\t"$6"\t"$1"_"$2"_"$7}' > $prefix"_exons_for_transcripts_neg.pc.pc.profile"

rm -rf x


	#less hg19v59_genes | awk '{s=s+1; print s"\t"$2"\t"$3"\t"$4"\t"$1"_"$5"_"$6}' > hg19v59_genes.peaks	
	
	
	#	less mm9v67_transcripts | awk 'BEGIN{l=2000}{if($6==1){print $3"\t"$4-l"\t"$4+(l-1)"\t"$6"\t"$1":"$7"_"$2":"$8} else {print $3"\t"$5-(l-1)"\t"$5+l"\t"$6"\t"$1":"$7"_"$2":"$8}}' > mm9v67_transcripts_TSS_2000
	
