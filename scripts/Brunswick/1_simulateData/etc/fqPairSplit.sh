#!/bin/bash
# Split a single, name sorted fasta file
# into seperate 'pairs' /1 and /2 .fa or fq files
#
# sh fqPairSplit.sh <FASTA.fq>
#

FASTQ=$1

# First Read Pair (Print lines 1-4 from 8)
sed -n '1~8{N;N;N;p}' $FASTQ | gzip -c -  > $FASTQ.1.gz

# Second Read Pair (Print lines 5-8 from 8)
sed -n '1~8{N;N;N;n;N;N;N;p}' $FASTQ | gzip -c - > $FASTQ.2.gz

