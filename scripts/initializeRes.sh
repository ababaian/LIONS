#!/bin/bash
# initializeRes.sh
# -----------------------------------------------
# source initializeRes.sh <resource_name>
# run from $BASE folder
# 
# For this <project> establish variables and symbolic links for 
# all annotations and resources required by LIONS.
#
# A symbolically linked up folder will also be created
# with all the resources as needed by each library in LIONS
# which will allow for parallel computing

# Resource Name:
	# $INDEX

# Resource Directory (from parent script)
	#RESOURCES="$BASE/resources/$INDEX"
	cd $RESOURCES

# Resource list----------------------------------

# RESOURCES/genome
	
	genomeFa="$INDEX.fa" # Genome Fasta
	genomeFai="$INDEX.fa.fai" # Genome Fasta Index

	bt_1="$INDEX.1.bt2" #Bowtie 2 index files for genome
	bt_2="$INDEX.2.bt2"
	bt_3="$INDEX.3.bt2"
	bt_4="$INDEX.4.bt2"
	bt_r1="$INDEX.rev.1.bt2"
	bt_r2="$INDEX.rev.2.bt2"
	bt_chr="$INDEX.chr.size"
	bt_bwa="$INDEX.bwa.names"

# RESOURCES/repeat
	repeatMasker="RepeatMasker.hg19.ucsc" # RepeatMasker UCSC download	
	rmSearch="ForChimericSearch_hg19" # LINE SINE LTR DNA Other elements only

# RESOURCES/annotation

	geneSet="$GENESET"
	# Will be used to make --->
	# refGene.bed: A Bed set of protein coding genes to intersect the assembly with. 


# Functions -------------------------------------

# Test if file exist and can be read
# Returns an error if it doesn't else it returns a 'S'all good' and continues
FCHECK_rs='if [ -s $FILE -a -r $FILE ]; then echo "... $FILE found."; else echo "     $FILE not found (empty or non-readable)."; echo " ===== ERROR 5: MISSING REQUISITE RESOURCE ===== "; echo Check README for more information; exit 2; fi'

# GENOME FILE CHECK -----------------------------
cd genome # go to genome directory

# ----- Check if genome.fa file is present
	FILE="$genomeFa"; eval $FCHECK_rs

# -----  Check if genome.fa.fai is present, if not generate it
if [ ! -r $genomeFai ];
then
	echo ' Genome fasta index file not found. Generating...'
	$SAMTOOLS faidx $genomeFa
	
	# Generate res.chr.size file
	cut -f1,2 $INDEX.fa.fai > $INDEX.chr.size
	
	# Generate INDEX.bwa.names
	cut -f1 $INDEX.fa.fai > tmp
	cut -f1 $INDEX.fa.fai | sed 's/chr//g' - > tmp2
	paste tmp2 tmp > $INDEX.bwa.names
	rm tmp*
	
fi

# -----  Check if bowtie2 index is present, if not generate it
if [ ! -r $bt_1 ];
then
	echo ' Bowtie2 index file not found. Generating...'
	$lBIN/bowtie-build $genomeFa $INDEX
fi

cd ..

# REPEAT FILE CHECK ------------------------------
# ***** Implement UCSC download/ File Conversion *****

# ANNOTATION FILE CHECK- -------------------------
# Reference UCSC Table File
# (download ucsc RefSeq table https://genome.ucsc.edu/cgi-bin/hgTables)
cd annotation

# Check if the Gene Set annotation file defined in parameter.ctrl exists
	FILE="$geneSet"; eval $FCHECK_rs

if [ ! -r refGene.bed ];
then
	echo " Protein coding gene subset not found. Generating..."

	# UCSC RefSeq Standard Column Order
	# to be extracted into the ouput file
	chr=3
	txStart=5
	txEnd=6
	name=13 # name2 (GENEID) used
	score=9 # Exon Count used
	strand=4

	# Remove ncRNA entries (those with GeneName containing 'unk')
	grep -P -v "\tunk\t" $geneSet  > pcSet.tmp

	# Cut works in ascending order only so use 3 cuts
	# to generate Bed file and paste them together

	cut -f$chr,$txStart,$txEnd,$name pcSet.tmp > ref1.tmp
	cut -f$score pcSet.tmp > ref2.tmp
	cut -f$strand pcSet.tmp > ref3.tmp

	paste ref1.tmp ref2.tmp ref3.tmp | sed 1d - > refGene.bed

	# NOTE: this generates a reference set with isoforms
	# or gene redundancy. To simplify it might be a good
	# idea to find unique gene boundaries since
	# that's what I really need here.

	# Clean Up
	rm *.tmp

fi

cd ..

# OTHER FILE CHECK -------------------------------

