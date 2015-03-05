#!/home/ababaian/bin/python3
# chimericReadSearch.py
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

# FILESTRUCTURES ------------------------------

# ExonList <list>
	# 'transcriptid' : 'ENST00001337.1'
	# 'geneid' : 'FAFF1'
	# 'rank_in_transcript' : 3	(Exon Number)
	# 'chr' : 'X'
	# 'start' : 15000000
	# 'end' : '15000300
	# 'strand' : -1
	# 'gene_biotype' : ''		(not implemented)



# FUNCTIONS -----------------------------------

# processReads
	# Extracts chimeric reads from bam
	# using exonTree and repeatTree
	# output to bed file (optional)

def processReads(samfile_path, exonTrees, repeatTrees, chimericBedFile):


	print("Processing all reads",file=sys.stderr)
	sys.stderr.flush()

	# Input Bam File
	samfile = pysam.Samfile( samfile_path, "rb" )
	
	readIterator = samfile.fetch()
	
	# Horrible code to extract total number of reads in BAM file
	readCount = sum([ eval('+'.join(l.rstrip('\n').split('\t')[2:]) ) for l in pysam.idxstats(samfile_path) ])
	
	print("	" + str(readCount) + " reads in BAM file",file=sys.stderr)
	
	localResults = {}
	count = 0
	nextPerc = 5.0
	
	# Go through every read in the bam file
	for read in readIterator:
		
		# Go through all read pairs that are:
			# Reads are paired
			# Pair on same chromosome
			# map quailty greater then zero (not multi-mapping)
			# [ add check to see read pair map quality as well ?]	
		#if (read.is_proper_pair and read.is_read1 and read.tid==read.mrnm):
		if (read.is_read1 and read.is_paired and read.tid==read.rnext and int(read.mapq)>0):	
			# Get chromosome for reads
			chr = samfile.getrname(read.tid)
			
			# Parse chromosome name to remove chr ('chr3' --> '3')
			# to remove chr (chr3 --> 3)
			chr = chr.replace("chr", "")

			# Skip reads not on chr 1-22,X,Y
			
			valid = list(range(1,23)) # chr1 - chr 22
			valid.append('X') # chrX
			valid.append('Y') # chrY
		
			if (chr not in str(valid)):
			# if read is not canonon chromosome skip
				continue

			# Get start coordinates for both reads
			start1 = read.pos
			start2 = read.mpos
			
			# [ Artem - Try working with spliced reads ]
			
			# Get end coordinates for both reads
			# [ Artem - Check if aligned length is used ]
			end1 = start1 + read.rlen
			end2 = start2 + read.rlen
			
			 
			# At the moment, it will fault here if there is no exon
			# information for the chromosome
			
			# Exon results is a list of exons (rows in the exon file) that
			# intersect with read1/2
			exon_results1 = exonTrees[str(chr)].findRange([start1,end1])
			exon_results2 = exonTrees[str(chr)].findRange([start2,end2])
			
			# Repeat results is a list of repeats (rows in the repeat file) that
			# intersect with read1/2
			repeat_results1 = repeatTrees[str(chr)].findRange([start1,end1])
			repeat_results2 = repeatTrees[str(chr)].findRange([start2,end2])
			
			# Get TRUE/FALSE if reads 1/2 intersect with exons or repeats
			e1 = (len(exon_results1) > 0)
			e2 = (len(exon_results2) > 0)
			r1 = (len(repeat_results1) > 0)
			r2 = (len(repeat_results2) > 0)
			
			# Classify read1 as D/E/R/.
			if (e1 and r1):
				type1 = "D"
			elif (e1):
				type1 = "E"
			elif (r1):
				type1 = "R"
			else:
				type1 = "."
				
			# Classify read2 as D/E/R/.
			if (e2 and r2):
				type2 = "D"
			elif (e2):
				type2 = "E"
			elif (r2):
				type2 = "R"
			else:
				type2 = "."
			
			# Sort (so "RE" becomes "ER")
			type = "".join(sorted(type1 + type2))
			
			# Is Read Chimeric?
			if chimericBedFile != 0 and ((e1 and r2) or (r1 and e2)):
				# Is Chy
				feature_start = min(start1, start2)
				feature_end = max(end1, end2)
				gap = abs(start1-start2)
				line = str(chr) + "\t" + str(feature_start) + "\t" + str(feature_end) + "\tchimericread\t960\t.\t" + str(feature_start) + "\t" + str(feature_end) + "\t0,0,250\t2\t" + str(read.rlen) + "," + str(read.rlen) + "\t0," + str(gap)
				chimericBedFile.write(line + "\n")
			
			# Zips up exon / repeat IDs
			pairs1 = list(zip(exon_results1, repeat_results2))
			pairs2 = list(zip(exon_results2, repeat_results1))
			pairs = pairs1 + pairs2
			
			# Use these (exonID,repeatID) pairs as key to dictionary
			# and store read type in dictionary
			for p in pairs:				
				result = localResults.get(p, [])
				result.append(type)
				localResults[p] = result
				
		count = count + 1
		perc = round((count/float(readCount))*100.0,1)
		
		# Print status to standard output
		if (perc >= nextPerc):
			print("	 " + str(perc) + "% (" + str(datetime.time(datetime.now())) + ")",file=sys.stderr)
			nextPerc = nextPerc + 5.0
		
	
	return localResults
	# end function 

# buildIntervalTree
def buildIntervalTree( features, min, max ):
	return it.intervalTree(features, 0, 1, min, max)

# buildIntervalTree_tuple
def buildIntervalTree_tuple(tuple):
	return buildIntervalTree(*tuple)

# ExistsNonZero
def file_exists(fpath):
	return True if os.path.isfile(fpath) and os.path.getsize(fpath) > 0 else False

#readThread
def readThread( samfile_path, Coords, Strand ):
# <samfile_path>: bam file input for analysis
# <Coords>: A 3-element object with Chromosome, Start, End of element
# <ReadLength>: integer of read length
# <Strand>: True = postive Strand, False = Negative strand
#
#           ___________________
#__________|_____TE____________|_______________ read1 read2	Call
#          |                   |
#          |   1       2       |
#            ;====---;====                      i	i	-
#          |                   | 
#	   |                   |
#        ;====---;====                          -       u	u
#          |                   |
#    ;====---;====                              -	u	u
#          |                   |
#;====---;====                 |                -       -*2     -
#          |                   |
#                            ;====---;====      d	- *1	d
#          |                   |
#                       ;====---;====           d	-	d
#          |                   |
#                   ;====---;====               d       i	d
#          |                   |
#          |                   |
#        ;====--------------;====               -*1	u	u
#          |                   |
#        ;====--------------------;====         -*1	-*2	-
#          |                   |

# Legend:
# ;	Leftmost position of read
# ===	Aligned read length
# ---	Internal sequence
# i	Read Mate is internal to TE
# -	Read Mate is uncounted/ implicit internal read
# d	Read Mate is downstream of TE
# u	Read Mate is upstream of TE

# Error *1, this implicit internal read is a downstream
# Error *2, this implicit internal read is an upstream

	# Input Bam File
	samfile = pysam.Samfile( samfile_path, "rb" )

	# Import reads overlapping the element coordinates
		# [0] = Chromosome , [1] = Start, [2] = End
	ParsedCoord = 'chr{0[0]}:{0[1]}-{0[2]}'.format(Coords)
	readIterator = samfile.fetch(region = ParsedCoord)

	# Initialize output
	upThread = 0
	downThread = 0
	internalThread = 0

	# Iterate through the reads
	# Assume +ve strand, flip results if negative strand
	for read in readIterator:

		# Paired reads on same chromosome only
		if (read.is_paired and read.tid==read.rnext and int(read.mapq)>0):

			# Accessing Mate Information
				# readMate = samfile.mate(read)
				# read.mpost = mate start position

			MateStart = read.mpos
			ReadLength = read.rlen

			if ( MateStart < int(Coords[1])): # pair upstream element boundry
				upThread = upThread + 1

			elif ( (int(MateStart) + int(ReadLength)) > int(Coords[2])): # pair downstream
				downThread = downThread + 1

			else: # internal pair, count but no worries...
				internalThread = internalThread + 1
	# Output
	if (Strand):
		localResults = ( upThread, downThread)
	else:
		localResults = ( downThread, upThread)

	return localResults

# End of readThread functio
# SCRIPT INITILIZATION ----------------------------------------------------------

if __name__ == '__main__':

	print("###########################",file=sys.stderr)
	print("# ChimericReadSearch v1.0 #",file=sys.stderr)
	print("###########################",file=sys.stderr)
	print("",file=sys.stderr)
	print("Data to StdOut, Msgs to StdErr",file=sys.stderr)
	print("",file=sys.stderr)

	#pool = multiprocessing.Pool(processes=6)
	
	if len(sys.argv) != 4 and len(sys.argv) != 5:
		print("USAGE: python " + sys.argv[0] + " <exon_file> <repeat_file> <bam_file> [optional_out_bed_file]",file=sys.stderr)
		print("example: python " + sys.argv[0] + " /projects/mbilenky/resources/hg18_genCodeV3/hg18_genCodeV3_exons /projects/mbilenky/resources/RepeatMasker/ForChimericSearch_hg18 in.bam out.bed",file=sys.stderr)
		sys.exit(1)


	#######################
	# Sort out parameters #
	#######################
	
	# Exon File
	exon_file_path = sys.argv[1]

	# Repeat Annotation File
	annotation_file_path = sys.argv[2]
	
	# Bam File to be analyzed
	samfile_path = sys.argv[3]

	# Optional output bed file
	if len(sys.argv) == 5:
		output_bed = True
		chimeric_bed_file_path = sys.argv[4]
	else:
		output_bed = False
	
	print("Parameters ===================================",file=sys.stderr )
	print(" Exon File: " + exon_file_path,file=sys.stderr)
	print(" Repeat File: " + annotation_file_path,file=sys.stderr)
	print(" Reads File: " + samfile_path,file=sys.stderr)
	print("",file=sys.stderr)	
	print("",file=sys.stderr)
	
	if not (os.path.exists(exon_file_path)):
		print("ERROR: Cannot access exon file",file=sys.stderr)
		print("Exiting...",file=sys.stderr)
		sys.exit(1)
	if not (os.path.exists(annotation_file_path)):
		print("ERROR: Cannot access annotation file",file=sys.stderr)
		print("Exiting...",file=sys.stderr)
		sys.exit(1)
	if not (os.path.exists(samfile_path)):
		print("ERROR: Cannot access reads file",file=sys.stderr)
		print("Exiting...",file=sys.stderr)
		sys.exit(1)
		
# SCRIPT CORE---------------------------------------------------------------------

	
	# RESOURCES -----------------------------------------------

# Pre-computed Exon Data	
	# Directory containing pre-computed IntervalTrees
	resTreeDir = '/home/ababaian/resources/chimeric/precomputed'
	
	# Exon File Name
	ExonFile = os.path.split(exon_file_path)[1] 
	
	# Precomputed IntervalTree (if it exists)
	resTreeExon = (resTreeDir + "/" + ExonFile + ".tree")

	# Precomputed ExonList
	resListExon = (resTreeDir + "/" + ExonFile + ".list")

	# Precomputed Trancsript Info
	resInfoExon = (resTreeDir + "/" + ExonFile + ".info")

# Pre-computed Annotation Data (repeat)
	# Annotation File Name
	AnnFile = os.path.split(annotation_file_path)[1]

	# Precomputed IntervalTree (if it exists)
	resTreeAnn = (resTreeDir + "/" + AnnFile + ".tree")

	# Precomputed AnnotationList
	resListAnn = (resTreeDir + "/" + AnnFile + ".list")

	print("Starting...",file=sys.stderr)
	sys.stderr.flush()

	###########################
	# Build Tree from Exons   #
	###########################

	# Check if ExonList and IntervalTree
	# exists for <exon_file>

	# Exception added for "assembly" runs
	# such that exon trees are computed novo

	if (not "assembly" in ExonFile ) and (file_exists(resTreeExon) and file_exists(resListExon) and file_exists(resInfoExon)):
	
		# Use precomputed interval tree
		
		inFile = open(resTreeExon, 'rb')
		exonTrees = pickle.load(inFile)
		inFile.close()

		# Use precomputed exon list
		inFile = open(resListExon, 'rb')
		exonList = pickle.load(inFile)
		inFile.close()

		# Use precomputed transcript info
		inFile = open(resInfoExon, 'rb')
		transcriptInfo = pickle.load(inFile)
		inFile.close()

		print(" Using precomputed Exon List and Interval Trees.",file=sys.stderr)
		sys.stderr.flush()

	else:
		# Create Interval Tree for Exon File

		# Read exon annotation (.bed)
		exon_file = open(exon_file_path, 'r')
		exonList = []
		exonCount = 0

		minpoint=int(9999999999999)
		maxpoint=int(-9999999999999)
		
		transcriptInfo = {}

		features = {}
		reader = csv.reader(exon_file, delimiter='\t')
		
		# Cycle through each exon annotation
		for row in reader:
			chr = row[2]	# Chromosome
			geneid = row[0]	# GeneID
			transcriptid = row[1]	# TranscriptID
			start = int(row[3])	# Start Location
			end = int(row[4])	# End Location
			strand = int(row[5])	# Strand
			rank_in_transcript = int(row[6])	# Exon number
			exons_in_gene = int(row[8])		# Exons in gene_total
			gene_biotype = row[7]	# Gene Type
			
			# Add entry for exon to dictionary
			exon = { 'geneid':geneid, 'transcriptid':transcriptid, 'chr':chr, 'start':start, 'end':end, 'strand':strand, 'rank_in_transcript':rank_in_transcript, 'gene_biotype':gene_biotype, 'exons_in_gene':exons_in_gene }
			
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
			print ("	 constructing chromosome: " + k, file=sys.stderr)
			exonTrees[k] = buildIntervalTree(features[k], minpoint, maxpoint)
			
		

		print("Exon tree building complete",file=sys.stderr)

		if (ExonFile != "assembly_exons_2"):

			# Save exon tree to resource directory
			# Open File Object for writing exon tree
			outFile = open(resTreeExon, 'wb')
			pickle.dump(exonTrees, outFile)
			outFile.close()		
			print("	Tree saved to resource directory",file=sys.stderr)
		
			# Save exon list to resource directory
			outFile = open(resListExon, 'wb')
			pickle.dump(exonList, outFile)
			outFile.close()		
			print("	List saved to resource directory",file=sys.stderr)
		
			# Save exon info to resource directory
			outFile = open(resInfoExon, 'wb')
			pickle.dump(transcriptInfo, outFile)
			outFile.close()
			print("	Info saved to resource directory",file=sys.stderr)
		
			sys.stderr.flush()
		
		#Else assembly is being run so don't save the exonTree

	# End of Exon Tree Building

	################################
	# Build tree from RepeatMasker #
	################################
	

	# Check if IntervalTree Exists for <ann_file>
	if (not "assembly" in ExonFile ) and (file_exists(resTreeAnn) and file_exists(resListAnn)):
		
		# Use precompiled interval tree
		inFile = open(resTreeAnn, 'rb')
		repeatTrees = pickle.load(inFile)
		inFile.close()

		# Use precompiled repeat
		inFile = open(resListAnn, 'rb')
		repeatList = pickle.load(inFile)
		inFile.close()
		
		print(" Using precomputed Repeat List and Interval Tree",file=sys.stderr)
		sys.stderr.flush()
	else:
		# Compute a new interval tree
		annotations_file = open(annotation_file_path, 'r')

		minpoint=int(9999999999999)
		maxpoint=int(-9999999999999)
		
		repeatList = []
		repeatCount = 0

		features = {}
		reader = csv.reader(annotations_file, delimiter='\t')
		for row in reader:

			chr = row[0].replace('chr','')
			
			if (not chr.isdigit() and chr != 'X' and chr != 'Y'):
				continue
			
			start = int(row[1])
			end = int(row[2])
			strand = row[3]
			r_name = row[4]
			r_class = row[5]
			r_family = row[6]
			
			repeat = {'chr':chr, 'start':start, 'end':end, 'strand':strand, 'name':r_name, 'class':r_class, 'family':r_family }
			
			minpoint = min(minpoint, start)
			maxpoint = max(maxpoint, end)
			
			feature = [start, end, repeatCount]
			repeatList.append(repeat)
			repeatCount = repeatCount + 1
			
			flist = features.get(chr, [])
			flist.append(feature)
			features[chr] = flist
			
		# Build a tree for each chromosome
		# This is the structure for indexing a line (chromosome coordinates)
		# Create one for each chromosome
		print("Building repeatTrees ( #" + str(len(features)) + ", min=" + str(minpoint) + ", max=" + str(maxpoint) + ")",file=sys.stderr)
		sys.stderr.flush()
		
		repeatTrees = {}
		# For each chromosome build a separate tree
		for k in sorted(features.keys()):
			print ("	 constructing chromosome: " + k, file=sys.stderr)
			repeatTrees[k] = buildIntervalTree(features[k], minpoint, maxpoint)
			
		print("Repeat trees complete",file=sys.stderr)
		
		# Open a file object to write tree
		outFile = open(resTreeAnn, 'wb')
		pickle.dump(repeatTrees,outFile)
		outFile.close()

		print("	 Tree saved to resource directory", file=sys.stderr)
		
		# Open a file object to write list
		outFile = open(resListAnn, 'wb')
		pickle.dump(repeatList,outFile)
		outFile.close()

		print("	 List saved to resource directory", file=sys.stderr)
		
		sys.stderr.flush()
	# End of Annotation Tree Interavl Building


	##############################
	# Do read search on BAM File #
	##############################

	results = {}
	
	if (output_bed):
		with open(chimeric_bed_file_path, 'w') as chimericBedFile:
			results = processReads(samfile_path, exonTrees, repeatTrees, chimericBedFile)
	else:
		results = processReads(samfile_path, exonTrees, repeatTrees, 0)
	
	print("Chimeric processing complete",file=sys.stderr)
	sys.stderr.flush()


	
	#############################################################
	# Calculate Exon/Repeat interactions,etc and output results #
	#############################################################
	
		
	print("",file=sys.stderr)
	print("",file=sys.stderr)
	print("Calculate Exon/Repeat Interactions",file=sys.stderr)
	print("",file=sys.stderr)
	
	sys.stderr.flush()
	
	# Iterate over every (exonID,repeatID) interaction in dictionary	
	for r in results.keys():
		
		# Pair of Exon/Repeat
		exon = exonList[r[0]]
		repeat = repeatList[r[1]]

	# General Information
		chr = exon['chr']
		estart = exon['start']
		eend = exon['end']
		rstart = repeat['start']
		rend = repeat['end']
		
		minstart = min(estart, rstart)
		maxend = max(eend, rend)
		coords = "chr" + str(chr) + ":" + str(minstart) + "-" + str(maxend)
		
		rstrand = repeat['strand']
		estrand = exon['strand']
		pos_strand = (estrand == 1)
		

	# Chimeric Read Classification
		# Get the list of read types:
		#	["ER","DR",etc] for this E/R interaction
		rlist = results[r]
		
		# Sum up read types
		type_counts = dict(Counter(rlist).items())
	
		#is_first_exon = (exon['rank_in_transcript'] == 1)
	
	# Exon / Repeat Interaction classification
		if (rstart >= estart and rend <= eend):
			er_interaction = "RInside"
		elif (rstart <= estart and rend >= eend):
			er_interaction = "EInside"
		elif (rstart <= estart and rend >= estart):
			er_interaction = "UpEdge" if (pos_strand) else "DownEdge"
		elif (rstart <= eend and rend >= eend):
			er_interaction = "DownEdge" if (pos_strand) else "UpEdge"
		elif (rend <= estart):
			er_interaction = "Up" if (pos_strand) else "Down"
		elif (rstart >= eend):
			er_interaction = "Down" if (pos_strand) else "Up"
		else:
			er_interaction = "Unknown"
		
	# Does repeat overlap with an exon?
		rexon_results = exonTrees[str(chr)].findRange([rstart,rend])
		is_exonic = "Yes" if (len(rexon_results)>0) else "No"
		
		# Summaries which exons are overlapped with repeat
		transcript_list = []
		for re_id in rexon_results:
			re = exonList[re_id]
			transcript_list.append(re['transcriptid'] + ":" + str(re['rank_in_transcript']))
		
		if len(transcript_list)>0:
			repeat_gene_id = ",".join(transcript_list)
		else:
			repeat_gene_id = '.'
		

	# Transcript Coordinates and Exon Rank Calculations

	# Extract/sort transcript exons coordinates
		transcriptid = exon['transcriptid']
		transcript = transcriptInfo.get(transcriptid, [])
		
		if pos_strand:
			transcript = sorted(transcript)
		else:
			transcript = sorted(transcript, reverse=True)
	
	# Repeat Rank Intron/Upstream
		# determine which intron repeat is in
		# or if it's upstream of exon 1 (-1)
		
		rRank = -1
		for ex in transcript:
			this_ex_rank = ex[0]
			this_ex_start = ex[1]
			this_ex_end = ex[2]
			if pos_strand and this_ex_start < rstart:
				rRank = this_ex_rank
			elif not pos_strand and this_ex_end > rend:
				rRank = this_ex_rank
				
	# Upstream Exon Coordinates
		# rank of current exon
		ExonRank = exon['rank_in_transcript']

		if ExonRank == 1: # First Exon
			# Special case: look for overlapping transcripts
			# to be implemented
			
			# tmp code
			UpExon = transcript[0]
			UpStart = UpExon[1]
			UpEnd = UpExon[2]

		else:
			# Strand Check
			if pos_strand:
				# Plus Upstream Exon Rank
				# (and python 0 index) 
				UpRank = ExonRank - 2

				UpExon = transcript[UpRank]

			else:
				# Minus Upstream Exon Rank
				UpRank = ExonRank - 2
				
				UpExon = sorted(transcript)[UpRank]
			
			# Extract upstream exon coords
			UpStart = UpExon[1]
			UpEnd = UpExon[2]


		# Example code to access previous exon information
		# prevExonInfo = transcript[int(exon['rank_in_transcript'])-1]
		# e.g. prevExonInfo[1] is the start coordinate
		
		# #ER	#DR	#DE	#DD
		ER_count = type_counts.get('ER', 0)
		DR_count = type_counts.get('DR', 0)
		DE_count = type_counts.get('DE', 0)
		DD_count = type_counts.get('DD', 0)
		total_count = ER_count + DR_count + DE_count + DD_count
	
	# Element Read Threading
	# Coutn upstream/downstream direction of read pairs within
	# an element

		# Define Repeat Element Coords
		RepeatCoords = (str(repeat['chr']), str(repeat['start']), str(repeat['end']))
		
		Threads = readThread( samfile_path, RepeatCoords, pos_strand)
		
		upThread = str(Threads[0])
		downThread = str(Threads[1])

	# Count Exons in Gene
	# In assembly; lots of 1 exon transcripts predicted; mostly junk
	# return number of exons in gene for filtration
	# 1: Single Exon Transcript
	# 2: Multi-Exon Transcript
		
		ExonInGene = str(exon['exons_in_gene'])
		
# OUTPUT WRITING ======================================================

		# Output col 1,2: Exon Information
		out = exon['geneid'] + ":" + exon['transcriptid'] + "\t" + str(exon['rank_in_transcript'])
		
		# Output col 3,4: Repeat name / coordinates
		out = out + "\t" + repeat['name'] + ":" + repeat['class'] + ":" + repeat['family'] + "\t" + coords
		
		# Output col 5-7: Interaction information
		out = out + "\t" + er_interaction + "\t" + is_exonic + "\t" + repeat_gene_id
		# Output col 9-12: Read Classifications
		out = out + "\t" + str(ER_count) + "\t" + str(DR_count) + "\t" + str(DE_count) + "\t" + str(DD_count) + "\t" + str(total_count)
		
		# Output col 13-17: Coordinates
		out = out + "\t" + str(chr) + "\t" + str(estart) + "\t" + str(eend) + "\t" + str(rstart) + "\t" + str(rend)
		
		# Output col 18,19: Strands
		out = out + "\t" + str(estrand) + "\t" + str(rstrand)
		
		# Output col 20: Repeat Rank
		out = out + "\t" + str(rRank)

		# Output col 21,22: Upstream Exon Coordinates
		out = out + "\t" + str(UpStart) + "\t" + str(UpEnd)

		# Output col 23,24: UpStream Threads & DownStream Threads
		out = out + "\t" + upThread + "\t" +  downThread
		
		# Output col 25: Exons in Transcript Model
		out = out + "\t" + ExonInGene
		print(out)
			
# Script over	
