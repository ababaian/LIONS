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
	export res=$1 # Resource name/ genome version name
	# also called $INDEX

# Resource Directory (from parent script)
	#RESOURCES="$BASE/resources/$INDEX"
	cd $RESOURCES

# Resource list----------------------------------

# RESOURCES/genome
	
	genomeFa="$res.fa" # Genome Fasta
	genomeFai="$res.fa.fai" # Genome Fasta Index

	bt_1="$res.1.bt2" #Bowtie 2 index files for genome
	bt_2="$res.2.bt2"
	bt_3="$res.3.bt2"
	bt_4="$res.4.bt2"
	bt_r1="$res.rev.1.bt2"
	bt_r2="$res.rev.2.bt2"
	bt_chr="$res.chr.size"

# RESOURCES/repeat
	repeatMasker="RepeatMasker.hg19.ucsc" # RepeatMasker UCSC download	
	rmSearch="ForChimericSearch_hg19" # LINE SINE LTR DNA Other elements only

# RESOURCES/chimera 


# Functions -------------------------------------

# Test if file exist and can be read
FCHECK_rs='if [ -s $FILE -a -r $FILE ]; then echo "... $FILE found."; else echo "     $FILE not found (empty or non-readable)."; echo " ===== ERROR 5: MISSING REQUISITE RESOURCE ===== "; exit 2; fi'

# GENOME FILE CHECK -----------------------------
cd genome # go to genome directory

# ----- Check if genome.fa file is present
	FILE="$genomeFa"; eval $FCHECK_rs

# -----  Check if genome.fa.fai is present, if not generate it
if [ ! -r $genomeFai ];
then
	echo ' Genome fasta index file not found. Generating...'
	$SAMTOOLS faidx $genomeFa
fi

# -----  Check if bowtie2 index is present, if not generate it
if [ ! -r $bt_1 ];
then
	echo ' Bowtie2 index file not found. Generating...'
	$BOWTIE2_BUILD $genomeFa $res
fi

cd ..

# REPEAT FILE CHECK ------------------------------

# OTHER FILE CHECK -------------------------------

