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


# Function to extract chimeric reads events from bamFile using exonTree and repeatTree
# and optionally output to a bed file
def processReads(samfile_path, exonTrees, repeatTrees, chimericBedFile):

	print("Processing all reads",file=sys.stderr)
	sys.stderr.flush()

	samfile = pysam.Samfile( samfile_path, "rb" )
	readIterator = samfile.fetch()
	
	# Horrible code to extract total number of reads in BAM file
	readCount = sum([ eval('+'.join(l.rstrip('\n').split('\t')[2:]) ) for l in pysam.idxstats(samfile_path) ])
	print(str(readCount) + " reads in BAM file",file=sys.stderr)
	
	localResults = {}
	count = 0
	nextPerc = 5.0
	for read in readIterator:
		
		# Check if read is in a pair (where both are mapped) to the same chromosome
		# [ Artem - Check if unique ]
		if (read.is_proper_pair and read.is_read1 and read.tid==read.mrnm):
			
			# Get chromosome for reads
			chr = samfile.getrname(read.tid)
			#chr = '1'
			#print chr
			
			# Skip reads on mitochondrial
			if (chr=='M'):
				continue
		
			# Get start coordinates for both reads
			start1 = read.pos
			start2 = read.mpos
			
			# [ Artem - Try working with spliced reads ]
			
			# Get end coordinates for both reads
			# [ Artem - Check if aligned length is used ]
			end1 = start1 + read.rlen
			end2 = start2 + read.rlen
			
			# 
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
		if (perc >= nextPerc):
			print(str(perc) + "% (" + str(datetime.time(datetime.now())) + ")",file=sys.stderr)
			nextPerc = nextPerc + 5.0
		
	
	return localResults
	

def buildIntervalTree( features, min, max ):
	return it.intervalTree(features, 0, 1, min, max)
	
def buildIntervalTree_tuple(tuple):
	return buildIntervalTree(*tuple)
	
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
	
	exon_file_path = sys.argv[1]
	annotation_file_path = sys.argv[2]
	samfile_path = sys.argv[3]
	if len(sys.argv) == 5:
		output_bed = True
		chimeric_bed_file_path = sys.argv[4]
	else:
		output_bed = False
	
	print("Exon File: " + exon_file_path,file=sys.stderr)
	print("Annotation File: " + annotation_file_path,file=sys.stderr)
	print("Reads File: " + samfile_path,file=sys.stderr)
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
		
	
	print("Starting...",file=sys.stderr)
	sys.stderr.flush()

	# /projects/mbilenky/mlc/Jake/chr1_only (SAM format)
	#exon_file_path='/projects/mbilenky/resources/hg18_genCodeV3/hg18_genCodeV3_exons'
	#annotation_file_path='/projects/mbilenky/resources/RepeatMasker/ForChimericSearch_hg18'
	#samfile_path = '/projects/mbilenky/mlc/Jake/colon_data/HS0988.bam'

	# less /projects/mbilenky/mlc/RepeatMasker/RetroTransposons	 | sed -e 's/\t/ /g' | grep "chr1 " | sed -e 's/chr//' > repeat_chr1
	# less /projects/03/genereg/projects/SOLEXA/resources/Ensembl_v65/mm9v65/mm9v65_exons | awk ' { if ( $3 == 1 ) print $ 0 }' > chr1_exons

	preloadRepeats = False
	preloadExons = False
	
	###########################
	# Build tree from Ensembl #
	###########################

	exon_file = open(exon_file_path, 'r')
	
	exonList = []
	exonCount = 0

	minpoint=int(9999999999999)
	maxpoint=int(-9999999999999)
	
	transcriptInfo = {}

	features = {}
	reader = csv.reader(exon_file, delimiter='\t')
	for row in reader:
		chr = row[2]
		
		geneid = row[0]
		transcriptid = row[1]
		start = int(row[3])
		end = int(row[4])
		strand = int(row[5])
		rank_in_transcript = int(row[6])
		gene_biotype = row[7]
		
		# Define dictionary for exon type
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
	# Create one for each chromosome
	print("Building exonTrees ( #" + str(len(features)) + ", min=" + str(minpoint) + ", max=" + str(maxpoint) + ")",file=sys.stderr)
	sys.stderr.flush()
	
	exonTrees = {}
	# For each chromosome build a separate tree
	for k in features.keys():
		exonTrees[k] = buildIntervalTree(features[k], minpoint, maxpoint)
		
	print("Exon trees complete",file=sys.stderr)
	sys.stderr.flush()

	################################
	# Build tree from RepeatMasker #
	################################
	

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
		repeatTrees[k] = buildIntervalTree(features[k], minpoint, maxpoint)
		
	print("Repeat trees complete",file=sys.stderr)
	sys.stderr.flush()

	##############################
	# Do read search on BAM File #
	##############################

	results = {}
	
	if (output_bed):
		with open(chimeric_bed_file_path, 'w') as chimericBedFile:
			results = processReads(samfile_path, exonTrees, repeatTrees, chimericBedFile)
	else:
		results = processReads(samfile_path, exonTrees, repeatTrees, 0)
	
	print("Processing complete",file=sys.stderr)
	sys.stderr.flush()


	
	#############################################################
	# Calculate Exon/Repeat interactions,etc and output results #
	#############################################################
	
	# Iterate over every (exonID,repeatID) interaction in dictionary
	for r in results.keys():
		exon = exonList[r[0]]
		repeat = repeatList[r[1]]
		
		# Get the list of chimeric read types: ["ER","DR",etc] for this E/R interaction
		rlist = results[r]
		
		# Sum up types
		type_counts = dict(Counter(rlist).items())
		
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
		
		#is_first_exon = (exon['rank_in_transcript'] == 1)
		
		# Pull out exon/repeat interaction
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
		
		# Transcribed repeat checker?
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
			
		# Calculate intron number of repeat
		transcriptid = exon['transcriptid']
		transcript = transcriptInfo.get(transcriptid, [])
		if pos_strand:
			transcript = sorted(transcript)
		else:
			transcript = sorted(transcript, reverse=True)
			
		prevExon = -1
		for e in transcript:
			this_e_id = e[0]
			this_e_start = e[1]
			this_e_end = e[2]
			if pos_strand and this_e_start < rstart:
				prevExon = this_e_id
			elif not pos_strand and this_e_end > rend:
				prevExon = this_e_id
				
		# Example code to access previous exon information
		# prevExonInfo = transcript[int(exon['rank_in_transcript'])-1]
		# e.g. prevExonInfo[1] is the start coordinate
		
		# #ER	#DR	#DE	#DD
		ER_count = type_counts.get('ER', 0)
		DR_count = type_counts.get('DR', 0)
		DE_count = type_counts.get('DE', 0)
		DD_count = type_counts.get('DD', 0)
		total_count = ER_count + DR_count + DE_count + DD_count
		
		out = exon['geneid'] + ":" + exon['transcriptid'] + "\t" + str(exon['rank_in_transcript'])
		out = out + "\t" + repeat['name'] + ":" + repeat['class'] + ":" + repeat['family'] + "\t" + coords
		out = out + "\t" + er_interaction + "\t" + is_exonic + "\t" + repeat_gene_id
		out = out + "\t" + str(ER_count) + "\t" + str(DR_count) + "\t" + str(DE_count) + "\t" + str(DD_count) + "\t" + str(total_count)
		out = out + "\t" + str(chr) + "\t" + str(estart) + "\t" + str(eend) + "\t" + str(rstart) + "\t" + str(rend)
		out = out + "\t" + str(estrand) + "\t" + str(rstrand) + "\t" + str(prevExon)
		print(out)
			
	
