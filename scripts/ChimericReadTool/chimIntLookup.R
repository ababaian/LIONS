# chimIntLookup.R
#
# Usage chimIntLookup.R 
# integrated into chimIntersect.sh
# INPUT: <chimID.tmp> <assBed.tmp>
# OUTPUT: <transcripts.bed>


# INPUT ===================================================

# Chimeric input list in transcriptID
# ordered from ChimericReadTool output
ChimID_path='chimID.tmp'

# Assembly Transcripts
# parsed into bed file
Assembly_path='assBed.tmp'


# IMPORT ==================================================
# Chimeric Output Ordered list of transcriptIDs
ChimOrder = read.csv(file = ChimID_path,
		header = F,
		sep = '\t',
		quote = '')

# Chimeric Output; unique transcript IDs
ChimID = unique(as.character(ChimOrder[,]))
ChimOrder = as.character(unlist(ChimOrder))

Assembly = read.csv(file = Assembly_path,
		header = F,
		sep = '\t',
		quote = '')

AssemblyID = as.character(Assembly[,4])

# Bed Format
#	1: chromosome
#	2: start
#	3: end
#	4: transcriptID
#	5: FPKM (score)
#	6: strand

# Match Chimeric TranscriptID to Assembly
# and return a bed file of ChimID ordered
# entries for assembly transcripts

# Chimeric Transcripts Only
	#BedOrder = match(ChimID, AssemblyID)
	#AssemblyBed = Assembly[BedOrder,]

# Entire Assembly
	AssemblyBed = unique(Assembly)

write.table(x = AssemblyBed,
	file = 'transcripts.bed',
	col.names = F,
	row.names = F,
	sep = '\t',
	quote = F)


# Run Bed Tools
# Intersects usign left outer join onto the transcripts.bed file
# output is itersect.bed
	system("bedtools intersect -loj -a transcripts.bed -b refGene.bed > intersectBed.tmp")

# Import Data Intersection
# Data is two left joined Bed files of 6x2 length
Intersection = read.csv(file = 'intersectBed.tmp',
		header = F,
		sep = '\t',
		quote = '')

# Parse and remove duplicate entries
#V1	1: chromosome 	[Assembly]
#V2	2: start	[Assembly]
#V3	3: end		[Assembly]
#V4	4: transcriptID	[Assembly]
#V5	5: FPKM (score)	[Assembly]
#V6	6: strand	[Assembly]
#V10	7: GeneID	[Reference]
#V12	8: Strand	[Reference]

Intersection = unique( Intersection[,c(1:6,10,12)] )

# Parse GeneID and GeneStrand to characters for aggregation
Intersection$V10 = as.character(Intersection$V10)
Intersection$V12 = as.character(Intersection$V12)


#Make it a toy
#Intersection = Intersection[500:600,]

# Aggregate GeneID and GeneStrand to unique 
Intersection = aggregate(cbind(V10,V12) ~ ., Intersection, paste, collapse = "; ")


compare = function ( A, B){
		output = (A == B)
		return(output)}

StrandCompare = function( InDF ){
	# Compare a single strand orientation
	# to multiple other strand orientations as list
	# to determine if they are Sense, Antisense,
	# or complex relationships

	# Origin	Comparison	Relationship
	# +		+		Sense
	# +		-		AntiSense
	# +		+ -		Complex
	# + 		+ +		S
	# +		- -		AS
	# +		+ + +		S
	# + 		+ + -		C
	# .		any		C

	# Coorce Origin to a vector
	Origin = as.character(InDF[1])
	
	# Split Comparison to list
	Comparison = strsplit( x = as.character(InDF[2]), split = "; ")

	if (Origin == "." ){
		# Strand of origin unknown
		return('u')

	}else if (Comparison == ".") {
		# Intergenic Element
		return('i')
	} else {
		# Create logical vector of Sense orientations
		Sense = lapply(Comparison, compare, A = Origin)
		# If it's not sense then it's antisense
		Antisense = lapply(Sense, "!")

		# Check to see if all orientations are the same
		Sense = unlist(lapply(Sense,all))
		Antisense = unlist(lapply(Antisense,all))
		Complex = !Sense & !Antisense
		
		if (Sense){
			return('s')
		}else if (Antisense) {
			return('as')
		}else if (Complex) {
			return('c')
		}else {
			return('?')
		}
}}

# Annotate Assembly / Reference Comparison using StrandCompare 
V13 = apply(Intersection[,c(6,8)], 1, StrandCompare)

Intersection = cbind(Intersection, V13)


# Output the intersection of the assembly and pcGene set
write.table(x = Intersection,
	file = 'Assembly.int.refGene',
	col.names = F,
	row.names = F,
	sep = '\t',
	quote = F)


# Order output in the original Chimeric.Output (ChimID) order
#
# if Xo is an ordered vector
#  and Yu is an unordered vector of corresponding elements
#  Yu[match(Xo,Yu)] restores order to Yu --> Yo

# Chimeric Order (ChimOrder)
# Gene/Info Order (IntIndx)

IntIndex = as.character(Intersection$V4)
sortOrder = match(ChimOrder, IntIndex)


# Use the sortOrder for rows

#Output = Intersection[sortOrder,c('V4','V10','V12','V13')]
#OutNames = c('AssemblyID','RefID', 'RefStrand', 'assXref')

Output = Intersection[sortOrder,c('V10','V12','V13')]
OutNames = c('RefID', 'RefStrand', 'assXref')

write.table(x = Output,
	file = 'appendCol.tmp',
	col.names = OutNames,
	row.names = F,
	sep = '\t',
	quote = F)

# End of script :D
