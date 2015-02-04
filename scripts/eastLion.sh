#!/bin/bash
# eastLion.sh
# 
# USAGE:
#	.eastLion.sh
#	ran from lions.sh
#
# Wrapper for processing .bam file into .fastq
# and running Bowtie2 and Tophat.
# Runs RNAseqPipeline and generates a per-library
# complete lions file and filters to TE-initiated transcripts
# only.
# Output: <library.lions> file
#

# CONTROL PANEL -------------------------------------------


