# sandbox.py
# Create interval trees for Exons/Repeats
# Detect and count chimeric reads
# 
# chimericReadSearch.py <Exons> <Repeats> <input.bam> <temporary.bed> > input.file
#
# ---------------------------------------------------------------------------------

# IMPORT -------------------------------------
import csv
import pysam
import sys
import intervalTree as it
import pickle
import os
import pprint
import threading
from collections import Counter
import multiprocessing
from datetime import datetime
import re


# FUNCTIONS -----------------------------------

# processReads
	# Extracts chimeric reads from bam
	# using exonTree and repeatTree
	# output to bed file (optional)


# buildIntervalTree
def buildIntervalTree( features, min, max ):
	return it.intervalTree(features, 0, 1, min, max)

# buildIntervalTree_tuple
def buildIntervalTree_tuple(tuple):
	return buildIntervalTree(*tuple)



# SCRIPT INITILIZATION ----------------------------------------------------------

	#######################
	# Sort out parameters #
	#######################
	
# Exon File
exon_file_path = '/projects/mbilenky/resources/hg19gc_v14/hg19gc_v14_exons'
resTreeExon = '/home/ababaian/resources/chimeric/trees/test.exon'

# SCRIPT CORE---------------------------------------------------------------------

# Exon File Name
ExonFile = os.path.split(exon_file_path)[1] 

exon_file = open(exon_file_path, 'r')

exonList = []
exonCount = 0

minpoint=int(9999999999999)
maxpoint=int(-9999999999999)
		
transcriptInfo = {}

features = {}
reader = csv.reader(exon_file, delimiter='\t')

# Cycle through each Gene annotation
for row in reader:
	chr = row[2]	# Chromosome
	geneid = row[0]	# GeneID
	transcriptid = row[1]	# TranscriptID
	start = int(row[3])	# Start Location
	end = int(row[4])	# End Location
	strand = int(row[5])	# Strand
	rank_in_transcript = int(row[6])	# Exon #
	gene_biotype = row[7]	# Gene Type
	
	# Add entry for exon to dictionary
	exon = { 'geneid':geneid, 'transcriptid':transcriptid, 'chr':chr, 'start':start, 'end':end, 'strand':strand, 'rank_in_transcript':rank_in_transcript, 'gene_biotype':gene_biotype }
	
	minpoint = min(minpoint, start)
	maxpoint = max(maxpoint, end)
			
	transcript = transcriptInfo.get(transcriptid, [])
	transcript.append( [ rank_in_transcript, start, end ] )
	transcriptInfo[transcriptid] = transcript
	
	feature = [start, end, exonCount]
	exonList.append(exon)
	exonCount = exonCount + 1
	
	flist = features.get(chr, [])
	flist.append(feature)
	features[chr] = flist


# Build a tree for each chromosome
# This is the structure for indexing a line (chromosome coordinates)
print("Building exonTrees ( #" + str(len(features)) + ", min=" + str(minpoint) + ", max=" + str(maxpoint) + ")",file=sys.stderr)
sys.stderr.flush()

exonTrees = {}

# For each chromosome build a separate tree
for k in features.keys():
	print ("     constructing chromosome: " + k, file=sys.stderr)
	exonTrees[k] = buildIntervalTree(features[k], minpoint, maxpoint)
	

# Save exon tree to resource directory
# Open File Object for writing exon tree
outFile = open(resTreeExon, 'wb')
pickle.dump(exonTrees,outFile)
outFile.close()

print("Exon tree building complete",file=sys.stderr)
print("Tree saved to resource directory",file=sys.stderr)
sys.stderr.flush()

# End of Exon Tree Building
