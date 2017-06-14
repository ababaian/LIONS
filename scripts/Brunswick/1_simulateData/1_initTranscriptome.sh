#!/bin/bash
#
# subsetTranscripts.sh <isoforms.fpkm_tracking> <top N transcripts> <ref_transcriptome.gtf>

# Controls

# Input File:
# Output from cufflinks quantification
# isoforms.fpkm_tracking
	INPUT=$1

# N Transcripts:
# Number of transcripts to extract
# (selecting the top N expressed transcripts)
	Nscript=$2

# Transcriptome GTF (original input)
# GENCODE v19 for now
	Transcriptome="gencode.v19.annotation.gtf"

# Script

# Extract FPKM values from file
	cut -f10 $INPUT > fpkm_only

# Sort (decreasing) by expression
	sort -g fpkm_only | tac - > sort_fpkm_only


# Find the Nth highest expression value

fpkmCutoff=$(sed -n "$Nscript"p sort_fpkm_only)

echo Threshold FPKM value is: $fpkmCutoff
echo will return all transcripts with greater expression

rm sort_fpkm_only fpkm_only

# In the original file extract
# all transcripts with greater then cutoff
# expression value

awk -v N=$fpkmCutoff '$10 >= N' $INPUT > expressed.genes_$fpkmCutoff

cut -f1 expressed.genes_$fpkmCutoff > exTranID.$fpkmCutoff 

echo successfully extracted $(wc -l expressed.genes_$fpkmCutoff) transcripts


# Unique Transcript Names in File

cut -f9 gencode.v19.annotation.gtf | cut -f2 -d';' |
sed 's/ transcript_id \"//g' - |
sed 's/\"//g' - |
sed 's/##.*/##/g' > tranID.gencode


# Make a list of Expressed Transcripts
	# It's super slow I know... sorry : (

for transcript in $(cat exTranID.$fpkmCutoff)
do
	#grep -n $transcript tranID.gencode >> exTranLines.$fpkmCutoff
	
	# Extract lines matching Transcript
	grep -n $transcript tranID.gencode | cut -f1 -d':' -  > exLINES	
	
	# Append those lines from the Gencode file into an expressed
	# transcriptome gtf file
	for LINE in $(cat exLINES)
	do
		sed "${LINE}q;d" $Transcriptome >> transcriptome.$fpkmCutoff.gtf
	done
		
done

# Clean up
rm exLINES exTranID.$fpkmCutoff expressed.genes_$fpkmCutoff

