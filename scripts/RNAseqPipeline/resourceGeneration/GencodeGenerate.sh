#!/bin/bash
set -e
set -o pipefail
set -x

db=$1
name=$2

# Example sh GencodeGenerate.sh hg18 hg18gc

echo "Generating Gencode data from db $db with name $name"

# exon file: geneId transcriptId chr exon_start exon_end strand rank_in_transcript gene_biotype

# gene file: geneId chr gene_start gene_end strand gene_biotype MGI_symbol/HGNC_symbol description

# transcript file: geneId transcriptId chr transcript_start transcript_end strand gene_biotype transcript_biotype MGI_symbol/HGNC_symbol description

if [ "$db" == "hg19" ]; then

	# for hg19
	mysql -h ucsc -u ucsc -pucsc -e "use $db -A; SELECT b.*, a.* FROM wgEncodeGencodeCompV14 b, wgEncodeGencodeAttrsV14 a WHERE a.transcriptID=b.name" > gencode.table
	
	less gencode.table | awk '!/bin/ {if($4=="+"){strand=1} else {strand=-1}; chr=gensub("chr","","g", $3); print $13"\t"$2"\t"chr"\t"$5"\t"$6"\t"strand"\t"$23"\t"$23"\tname\tdescr"}' | sort -k1,1 | uniq > $name"_transcripts"

	less $name"_transcripts" | cut -f1,3-7,9- | sort -k1,1 -k2,2 -k3,3n -k4,4n | uniq | awk '{if($1==id){e=$4} else {if(id!=null){print id"\t"$2"\t"s"\t"e"\t"$5"\t"$6"\t"$7"\t"$8}; id=$1;s=$3;e=$4;}} END{print id"\t"$2"\t"s"\t"e"\t"$5"\t"$6"\t"$7"\t"$8}' > $name"_genes"

	less gencode.table | cut -f2-6,10,11,13 | sort -k8,8 | uniq | awk '{n=split($6,s,",");n=split($7,e,",");chr=gensub("chr","","g", $2); if($3=="+"){strand=1;for(i=1;i<n;i++){print $8"\t"$1"\t"chr"\t"s[i]"\t"e[i]"\t"strand"\t"i"\t"$23}} else {strand=-1;for(i=1;i<n;i++){print $8"\t"$1"\t"chr"\t"s[i]"\t"e[i]"\t"strand"\t"(n-i)"\t"$23}}}' > $name"_exons"
	
elif [ "$db" == "hg18" ]; then

	# for hg19
	mysql -h ucsc -u ucsc -pucsc -e "use $db -A; SELECT b.*, a.* FROM wgEncodeGencodeAutoV3 b, wgEncodeGencodeClassesV3 a WHERE a.name=b.name" > gencode.table
	mysql -h ucsc -u ucsc -pucsc -e "use $db -A; SELECT b.*, a.* FROM wgEncodeGencodeManualV3 b, wgEncodeGencodeClassesV3 a WHERE a.name=b.name" >> gencode.table
	mysql -h ucsc -u ucsc -pucsc -e "use $db -A; SELECT b.*, a.* FROM wgEncodeGencodePolyaV3 b, wgEncodeGencodeClassesV3 a WHERE a.name=b.name" >> gencode.table
	
	
	less gencode.table | awk '!/bin/ {if($4=="+"){strand=1} else {strand=-1}; chr=gensub("chr","","g", $3); print $13"\t"$2"\t"chr"\t"$5"\t"$6"\t"strand"\t"$19"\t"$19"\tname\tdescr"}' | sort -k1,1 | uniq > $name"_transcripts"

	less $name"_transcripts" | cut -f1,3-7,9- | sort -k1,1 -k2,2 -k3,3n -k4,4n | uniq | awk '{if($1==id){e=$4} else {if(id!=null){print id"\t"$2"\t"s"\t"e"\t"$5"\t"$6"\t"$7"\t"$8}; id=$1;s=$3;e=$4;}} END{print id"\t"$2"\t"s"\t"e"\t"$5"\t"$6"\t"$7"\t"$8}' > $name"_genes"

	less gencode.table | cut -f2-6,10,11,13 | sort -k8,8 | uniq | awk '{n=split($6,s,",");n=split($7,e,",");chr=gensub("chr","","g", $2); if($3=="+"){strand=1;for(i=1;i<n;i++){print $8"\t"$1"\t"chr"\t"s[i]"\t"e[i]"\t"strand"\t"i"\t"$19}} else {strand=-1;for(i=1;i<n;i++){print $8"\t"$1"\t"chr"\t"s[i]"\t"e[i]"\t"strand"\t"(n-i)"\t"$19}}}' > $name"_exons"

	#exit 1
	
elif [ "$db" == "hg18X" ]; then

	# for hg18
	mysql -h ucsc -u ucsc -pucsc -e "use $db -A; SELECT b.* FROM wgEncodeGencodeAutoV3 b" > gencode.table
	mysql -h ucsc -u ucsc -pucsc -e "use $db -A; SELECT b.* FROM wgEncodeGencodeManualV3 b" >> gencode.table
	
	less gencode.table | cut -f2-6,10,11,13 | sort -k8,8 | uniq | awk '{n=split($6,s,",");n=split($7,e,",");chr=gensub("chr","","g", $2); if($3=="+"){strand=1;for(i=1;i<n;i++){print $8"\t"$1"\t"chr"\t"s[i]"\t"e[i]"\t"strand"\t"i"\t"$23}} else {strand=-1;for(i=1;i<n;i++){print $8"\t"$1"\t"chr"\t"s[i]"\t"e[i]"\t"strand"\t"(n-i)"\t"$23}}}' > $name"_exons"
	
	exit 1
	
else
	echo "Script can only process hg18/hg19"
	exit 1
fi

#mysql -h ucsc -u ucsc -pucsc -e "USE hg19 -A; SELECT b.name, b.chrom, b.txStart, b.txEnd, b.strand, a.geneType, a.geneName, 'description' FROM wgEncodeGencodeBasicV14 b, wgEncodeGencodeAttrsV14 a WHERE a.transcriptID=b.name" | tail -n+2 | awk '{ chr=gensub("chr","","g", $2); strand=1; if ($5=="-") strand=-1; print $1"\t"chr"\t"$3"\t"$4"\t"strand"\t"$6"\t"$7"\t"$8}' > $db"_genes"

#awk '{print $1"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$6"\t"$7"\t"$8}' $db"_genes"  > $db"_transcripts"

#mysql -h ucsc -u ucsc -pucsc -e "USE hg19 -A; SELECT b.name, b.name, b.chrom, b.txStart, b.txEnd, b.strand, a.geneType, a.geneType, a.geneName, 'description' FROM wgEncodeGencodeBasicV14 b, wgEncodeGencodeAttrsV14 a WHERE a.transcriptID=b.name" | tail -n+2 > $db"_transcripts"


#less gencode.table | awk '!/bin/ {if($4=="+"){strand=1} else {strand=-1}; chr=gensub("chr","","g", $3); print $13"\t"$2"\t"chr"\t"$5"\t"$6"\t"strand"\t"$23"\t"$23"\tname\tdescr"}' | sort -k1,1 | uniq > $name"_transcripts"

#less $name"_transcripts" | cut -f1,3-7,9- | sort -k1,1 -k2,2 -k3,3n -k4,4n | uniq | awk '{if($1==id){e=$4} else {if(id!=null){print id"\t"$2"\t"s"\t"e"\t"$5"\t"$6"\t"$7"\t"$8}; id=$1;s=$3;e=$4;}} END{print id"\t"$2"\t"s"\t"e"\t"$5"\t"$6"\t"$7"\t"$8}' > $name"_genes"


#less gencode.table | cut -f2-6,10,11,13 | sort -k8,8 | uniq | awk '{n=split($6,s,",");n=split($7,e,",");chr=gensub("chr","","g", $2); if($3=="+"){strand=1;for(i=1;i<n;i++){print $8"\t"$1"\t"chr"\t"s[i]"\t"e[i]"\t"strand"\t"i"\t"$23}} else {strand=-1;for(i=1;i<n;i++){print $8"\t"$1"\t"chr"\t"s[i]"\t"e[i]"\t"strand"\t"(n-i)"\t"$23}}}' > $name"_exons"

BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
sh $BASE/EnsemblResourceGenerator.sh $name

