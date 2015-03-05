#!/bin/bash
# buildResourceGTF.sh
#
# USAGE:
# sh buildResourceGTF.sh <assembly.gtf> <resource_name>
#
# Input is a Cufflinks assembly; constructs a resource
# for that file for RNAseqAnalysis // ChimericToolAnalysis
#
#

# INPUT ---------------------------------------------------

GTF=$1
NAME=$2

echo "Generating resource data for $GTF with the name $NAME"

# ---------------------------------------------------------

# GTF File Organization
#	1:	Chromosome name
#	2:	Source of transcript 
#	3:	Feature <transcript> <exon>
#	4:	START
#	5:	END
#	6:	score; most abundant isoform=1000
#	7:	Strand_guess
#	8:	frame [not used]
#	9:	Atributes 'name1 "value1"; name2 "value2"
#	-:		gene_id
#	-:		transcript_id
#	-:		FPKM
#	-:		frac
#	-:		conf_lo [Lower confidance interval]
#	-:		conf_hi [Higher confidence interval]
# ---------------------------------------------------------

# Exon File Organization
#	1:	geneId
#	2:	transcriptID
#	3:	chr
#	4:	exon_start
#	5:	exon_end
#	6:	strand
#	7:	rank_in_transcript
#	8:	gene_biotype
# ---------------------------------------------------------

# Gene File Organization
#	1:	geneId
#	2:	chr
#	3:	gene_start
#	4:	gene_end
#	5:	strand
#	6:	gene_biotype
#	7:	MGI_symbol/HGNC_symbol
#	8:	description
# ---------------------------------------------------------

# Transcript File Organization
#	1:	geneId
#	2:	transcriptId
#	3:	chr
#	4:	transcript_start
#	5:	transcript_end
#	6:	strand
#	7:	gene_biotype
#	8:	transcript_biotype
#	9:	MGI_symbol/HGNC_symbol
#	10:	description


# Script Core --------------------------------------------

# Seperate transcript and exon annotation
# into two seperate temporary files

	# Temporary files
	EXONS="$NAME.ex.tmp"
	TRANS="$NAME.tr.tmp"

	# Parse Exons/Transcripts for easy file-building
	grep 'exon' $GTF | sed 's/"//g' - | sed 's/; /\t/g' - |
	sed 's/gene_id //g' - | sed 's/transcript_id //g' - |
	sed 's/exon_number //g' - | sed 's/chr//g' - > $EXONS
	
	
	grep -P '\ttranscript' $GTF | sed 's/"//g' - | sed 's/; /\t/g' - |
	sed 's/gene_id //g' - | sed 's/transcript_id //g' - |
	sed 's/exon_number //g' - | sed 's/chr//g' > $TRANS
	


# Build Exon File -------------------------------

	awk '{print ($9,$10,$1,$4,$5,$7,$11,"NA")}' $EXONS |
	sed 's/+/1/g' - | sed 's/-/-1/g' - |
	sed 's/ /\t/g' - |
	sed 's/\t\.\t/\t0\t/g' - > "$NAME"_exons

# Build Transcript File -------------------------

	awk '{print ($9,$10,$1,$4,$5,$7,$2,"NA","NA","NA")}' $TRANS |
	sed 's/+/1/g' - | sed 's/-/-1/g' - |
	sed 's/ /\t/g' - |
	sed 's/\t\.\t/\t0\t/g' - > "$NAME"_transcripts


# Build Gene File -------------------------

	awk '{print ($9,$1,$4,$5,$7,$2,"NA","NA")}' $TRANS |
	sed 's/+/1/g' - | sed 's/-/-1/g' - |
	sed 's/ /\t/g' - |
	sed 's/\t\.\t/\t0\t/g' - > "$NAME"_genes

# Clean Up!
#rm $EXONS $TRANS



# Assemble Reference Set ----------------------------------

echo 'Precursor files generated'
echo ' Running EnsembleResourceGenerated.sh'

	sh EnsemblResourceGenerator.sh $NAME

# Cleanup
rm *.tmp
# End of script
