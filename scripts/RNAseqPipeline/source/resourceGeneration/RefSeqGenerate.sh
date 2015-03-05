#!/bin/bash
set -e
set -o pipefail
set -x

db=$1
name=$2

echo "Generating Gencode data from db $db with name $name"

# exon file: geneId transcriptId chr exon_start exon_end strand rank_in_transcript gene_biotype

# gene file: geneId chr gene_start gene_end strand gene_biotype MGI_symbol/HGNC_symbol description

# transcript file: geneId transcriptId chr transcript_start transcript_end strand gene_biotype transcript_biotype MGI_symbol/HGNC_symbol description

# less mm9rs_0313_exons | awk ' {print "chr"$3"\t"$4"\t"$5"\t"$1":"$2} ' > mm9rs_0313_exons.bed


mysql -h ucsc -u ucsc -pucsc -e "use $db -A; select * from refGene" > refGene.table

less refGene.table | awk '!/bin/ {if($4=="+"){strand=1} else {strand=-1}; if (substr($2,0,2)=="NM"){type="protein_coding"} else{type="non_coding"}; chr=gensub("chr","","g", $3); print $13"\t"$2"\t"chr"\t"$5"\t"$6"\t"strand"\t"type"\t"type"\tname\tdescr"}' | sort -k1,1 | uniq > $name"_transcripts"

less $name"_transcripts" | cut -f1,3-7,9- | sort -k1,1 -k2,2 -k3,3n -k4,4n | uniq | awk '{if($1==id){e=$4} else {if(id!=null){print id"\t"$2"\t"s"\t"e"\t"$5"\t"$6"\t"$7"\t"$8}; id=$1;s=$3;e=$4;}} END{print id"\t"$2"\t"s"\t"e"\t"$5"\t"$6"\t"$7"\t"$8}' > $name"_genes"


less refGene.table | cut -f2-6,10,11,13 | sort -k8,8 | uniq | awk '{n=split($6,s,",");n=split($7,e,",");if (substr($1,0,2)=="NM"){type="protein_coding"} else{type="non_coding"}; chr=gensub("chr","","g", $2); if($3=="+"){strand=1;for(i=1;i<n;i++){print $8"\t"$1"\t"chr"\t"s[i]"\t"e[i]"\t"strand"\t"i"\t"type}} else {strand=-1;for(i=1;i<n;i++){print $8"\t"$1"\t"chr"\t"s[i]"\t"e[i]"\t"strand"\t"(n-i)"\t"type}}}' > $name"_exons"

sh /home/jlever/bin/EnsemblResourceGenerator.sh $name

