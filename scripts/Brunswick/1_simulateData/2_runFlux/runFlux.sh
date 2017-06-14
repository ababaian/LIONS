#!/bin/bash
#
# runFlux.sh <parameter file .par>
#

# Run flux-simulator
# full run
#
flux-simulator -p $1

# Compress the output .bed file
gzip *.bed 

# Split and compress the fastq output
FASTQ=$(ls *.fastq)

# First Read Pair (Print lines 1-4 from 8)
sed -n '1~8{N;N;N;p}' $FASTQ | gzip -c -  > $FASTQ.1.gz

# Second Read Pair (Print lines 5-8 from 8)
sed -n '1~8{N;N;N;n;N;N;N;p}' $FASTQ | gzip -c - > $FASTQ.2.gz

rm $FASTQ

