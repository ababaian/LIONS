# chimGroup.R
# 
# Usage:
# Rscript chimGroup.R <project.lions> <input.list> <rmStats.Rdata Path> <Recurruance Cutoff> <Specificty Cutoff>
#
# Group chimera detected from each library into biological groups
# and report overlapping/depleted events between the two:
# 

# CONTROL PANEL ======================================================

  STDIN = commandArgs(trailingOnly = TRUE)

  pLIONS = STDIN[1] # Input LIONS project file
   pNAME = as.character(strsplit(pLIONS, '.lions')) # project+run name

  stdOUTPUT = paste(pNAME, '.rslions', sep='') # Grouped Project Lions File
  
  invOUTPUT = paste(pNAME, '.inv.rslions', sep='') # Inversed Grouped File

  GROUP_LIST = STDIN[2]
    # by convention this file is the 3 column CSV file with:
    # 1: Library Name, 2: Bam Path (on GSC), 3: Group Number
    # Group Numbers: [Although assigned in parameters explicitely]
    #   1: Normal/Baseline
    #   2: Cancer/Group I
    #   3: Other/ non-grouped samples
    #   ...

  RMDATA_PATH= STDIN[3]

# Parameters

  # Multiple Cancer Library - retain
    # How many cancer libraries (2) is the read in to be retained
    # Equal to or greater then this value
    # parameter.ctrl variable: cgGroupRecurrence
  
    MultiCan = as.numeric(STDIN[4])

  # Non Normal - deplete
    # Present in how many normal libraries to be depleted
    # Greater then or equal to this value
    # parameter.ctrl variable: cgSpecificity

    NonNormal = as.numeric(STDIN[5])

# BEGIN SCRIPT LOOPING ===============================================
# two iteration for loop to run
# analysis of Group 1 v 2
# and         Group 2 v 1 (inverse)

for  (RUN in c(1:2)) {

	if (RUN == 1){

	print(" Running standard analysis" )
	# Standard:
	OUTPUT = stdOUTPUT

	# Cancer specific vs normal control
	 # Normal Group (control)
	 normGroup = 1
	 
	 # Cancer Group
	 canGroup = 2

	 } else {
	 print(" Running inverse analysis" )
	# Inverse:
	OUTPUT = invOUTPUT

	# Normal specific vs. cancer control
	 # Normal Group
	 normGroup = 2

	 # Cancer Group (control)
	 canGroup = 1
	}

# INITIALIZATION ====================================================

# Import

  Chimera = read.csv(file=pLIONS, header=T, sep='\t')
  
  # Sub-sample Chimera
    #Sample = rep(0,times = nrow(Chimera))
    #Sample = unlist(lapply(Sample,FUN=sample,x=1:10,size=1))
    #Sample = which(Sample == 1)
    #Chimera = Chimera[Sample,]

  Groups = read.csv(file=GROUP_LIST, header=F, sep='\t')

# FUNCTIONS =========================================================

GroupAssign = function(Library){
  # INPUT:
  # Library is the library name from Chimera['LIBRARY']
  # GroupMat is the Groups matrix; col 1: libary, col 3: group number
  # Use sapply to work it
  
  Match = which(as.character(Library) == as.character(Groups[,1]))
  Assign = as.numeric(Groups[Match,3])
  
  return(Assign)
}

MatchEntry = function (ChimeraRow, GroupN){
  # INPUT:
  # For a given ChimeraRow, find a TE/Exon interaction match
  # within the Chimera table and count the occurances
  # Assembly: use RepeatID to identify commonalities
  
  # Match Genes
  #ExonMatch = which(ChimeraRow['transcriptID'] == Chimera[,'transcriptID'])
  
  # Match Coords (TE)
  #TEMatch = which(ChimeraRow['coordinates'] == Chimera[,'coordinates'])
  
  # Match RepeatID
  TEMatch = which(ChimeraRow['RepeatID'] == Chimera[,'RepeatID'])
  
  # Match Group
  GMatch = which(GroupN == as.numeric(Chimera[,'Group']) )
  
  # Intersect sets to find unique TE/Exon interactions per group
  #  Interactions = intersect(ExonMatch,TEMatch)
  #  GroupInteractions = intersect(Interactions,GMatch)
  
  GroupInteractions = intersect(TEMatch,GMatch)
  
  # Return count matching all 3 criteria
  HitCount = length(GroupInteractions)
  
  return(HitCount)
  
}

OrientAssembly = function (AXR){
  # Input a character vector of assXref
  # orientations for a group of chimera
  # Implement the following dictionary
    # [i/u], [s,as,c] = [s,as,c]
    # [s] , [as] = c
    # [c] , [s,as] = c
    # [i] , [u] = i
  
    if ( "c" %in% AXR ) {
      # If there is a complex interaction
      # it is dominant to all others
      interaction = "c"
      
    } else if ( "s" %in% AXR ) {
      # Contains a sense annotaiton
      if ( "as" %in% AXR ) {
        # Contains both s & as: complex
        interaction = "c"
      } else {
        # Contains s but not as
        interaction = "s"
      }
    } else if ("as" %in% AXR ) {
      # Contains as but not s
      interaction = "as"
    } else if ("i" %in% AXR ) {
      # Contains i but not c,s,as
      interaction = "i"
    } else {
      # contains only u
      interaction = "u"
    }
    return(interaction)
}

Recurrence = function(Chimera) {
  
  # Recurrence Matrix
  # For each library in the Chimera Table,
  # find the recurrence of it's chimera
  # with each other library.
  
  # Library Names
  LIB_VECTOR = c(as.vector(unique(Chimera[,'LIBRARY'])))
  nLIB = length(LIB_VECTOR)
  
  # Chimera for each library
  # For every Library (col) how many of the chimeric
  # calls are inLibrary (row)?
  
  RecMatrix = matrix(0, nrow = nLIB, ncol = nLIB)
  colnames(RecMatrix) = LIB_VECTOR
  rownames(RecMatrix) = paste('in',LIB_VECTOR,sep='')
  
  # Establish a list of which chimera are in which library
  # for iterative searching
  ChimeraLookup = list()
  for (Count in c(1:nLIB)){
    LibLookup = which(Chimera$LIBRARY == LIB_VECTOR[Count])
    ChimeraLookup = c(ChimeraLookup, list(LibLookup))
  }
  
  # Populate Recurrence Matrix
  for (Key in c(1:nLIB)) { # for every column (Key Library)
    
    KeyChimera = as.vector(Chimera$RepeatID[unlist(ChimeraLookup[Key])])
    
    for (Query in c(1:nLIB)){ # go through each row (Query library)
      QueryChimera = as.vector(Chimera$RepeatID[unlist(ChimeraLookup[Query])])
      
      Hits = length(which(KeyChimera %in% QueryChimera))
      
      RecMatrix[Key,Query] = Hits
      
    }
  }
  # Done all that silly looping!
  return(RecMatrix)
}

# SCRIPT CORE =======================================================

# Find recurrent elements -------------------------------------------

# Remove Duplicate Entries within a library from Chimera Table

  # Unique Entries
  ID = paste(Chimera[,'RepeatID'],Chimera[,'LIBRARY'],sep='~')
  # Check for Duplicates, True = Unique or First Entry, False = Duplicate
  IDdup = !duplicated(ID)
  # Remove Duplicates
  Chimera = Chimera[IDdup,]

# Group Chimera
  # Assign biological catagories to Chimera from Groups

  Chimera['Group'] = sapply(Chimera[,'LIBRARY'], FUN = GroupAssign)
  
# Count Groups in Chimera
  Chimera['Normal_Occ'] = apply(X = Chimera, MARGIN = 1, FUN = MatchEntry, GroupN = normGroup)
  Chimera['Cancer_Occ'] = apply(X = Chimera, MARGIN = 1, FUN = MatchEntry, GroupN = canGroup)

# Filter Chimeric List for recurrent elements

  Filter_Normal = which( Chimera[,'Normal_Occ'] <= NonNormal)
  Filter_Cancer = which( Chimera[,'Cancer_Occ'] >= MultiCan)

  Filter = intersect(Filter_Normal,Filter_Cancer)

  ChimFiltered = Chimera[Filter,]

# Parse for output
# Parse ChimeraOutput such that every TE/Exon is 1 row, INPUT changed to
  # comma seperated list of libraries containing this TE/Exon

  ChimeraOutput = data.frame()

# Library Parse
# Create a list of libraries with
# RepeatID Present

Count = 0
Go = T

while (Go == T){
  # Add to Counter
  Count = Count + 1

  # Initialize Output Column
  #Chimera[1,c('transcriptID','exonRankInTranscript','repeatName','RepeatRank','coordinates','ER_Interaction','Normal_Occ','Cancer_Occ')]
  OutCol = c('transcriptID', 'exonRankInTranscript','EStrand',
             'RefID','RefStrand','assXref', 
             'repeatName', 'RepeatRank',
             'coordinates', 'ER_Interaction','RepeatID',
             'Normal_Occ','Cancer_Occ')
  
  OutputRow = ChimFiltered[1,OutCol]
  
  # Match RepeatID
  Matches = which(ChimFiltered[1,'RepeatID'] == ChimFiltered[,'RepeatID'])
  
#  # Match Genes
#  ExonMatch = which(ChimFiltered[1,'transcriptID'] == ChimFiltered[,'transcriptID'])
#  # Match Coords (TE)
#  TEMatch = which(ChimFiltered[1,'coordinates'] == ChimFiltered[,'coordinates'])
#  Matchs = intersect(ExonMatch, TEMatch)
  
  # Parse Libraries
  Libraries = ChimFiltered[Matches, 'LIBRARY']
  OutputRow['Library'] = paste(Libraries,collapse=";")
  
  # Parse RefID to unique RefID entries
  RefID = as.character(unlist(ChimFiltered[Matches, 'RefID']))
  RefID = unique(as.character(unlist(strsplit(RefID, split="; "))))
  RefID = paste(RefID, collapse="; ")
  OutputRow['RefID'] = RefID

  # Parse assXref interaction
  AXR = as.character(unlist(ChimFiltered[Matches, 'assXref']))
  AXR = unique(as.character(unlist(strsplit(AXR, split="; "))))
  AXR = OrientAssembly(AXR)
  OutputRow['assXref'] = AXR
  
  
  # Write to output
  ChimeraOutput = rbind(ChimeraOutput,OutputRow)
  
  # Remove matches from list
  ChimFiltered = ChimFiltered[-Matches,]
  
  if ( nrow(ChimFiltered) == 0){ Go = F } # When all entries are ran through stop
}

# Sort Output by coordinates
  ChimeraOutput = ChimeraOutput[order(ChimeraOutput[,'coordinates']),]
  # Remove Duplicate Entries

  # Unique Entries
  ID = ChimeraOutput[,'RepeatID']
  # Check for Duplicates, True = Unique or First Entry, False = Duplicate
  IDdup = !duplicated(ID)
  # Remove Duplicates
  ChimeraOutput = ChimeraOutput[IDdup,]

# Clean Up
#  \    @ /    
#   \ @  /
#    ----
# rm( ChimFiltered, Count, AXR, Filter, Filter_Cancer, Filter_Normal,
#      Go, ID, IDdup, Matches, MultiCan, RefID, NonNormal)

RecMatrix = Recurrence(Chimera) # custom function
# Normalize Recurrence Matrix based on Total number of hits
  RecMax = diag(RecMatrix)
  RecMatrixNorm = 100*(RecMatrix / RecMax)

# Subset LTR elements
  # Which elements are LTR
  LTR_MATCH = which("LTR" == matrix(unlist(strsplit(as.vector(
    Chimera$repeatName),split=":")),ncol=3,byrow=T)[,2])
  
#  LtrMatrix = Recurrence(Chimera[LTR_MATCH,])
#  LtrMax = diag(LtrMatrix)
#  LtrMatrixNorm = 100*(RecMatrix / RecMax)

# Library Statistics ------------------------------------------------

# Load Repeat Masker Database
  load(RMDATA_PATH)

# Parse RMDB to RepeatID Vector
# chrX:<start>
  RMDB_id = paste(RMDB$genoName, RMDB$genoStart, sep=":")

# Link Chimera Table Elements with their RMDB element
  Chimera$RMDB = unlist(lapply(Chimera$RepeatID,
                               FUN=function(Key,Query=RMDB_id){
                                 return(which(Key == Query))}
                ))

# Chimera in Library List
  # Library Names
  LIB_VECTOR = as.vector(unique(Chimera[,'LIBRARY']))

  # Repeat Masker Names
  repName = sort(unique(RMDB$repName))
  repClass = sort(unique(RMDB$repClass))

  
  FIRST = T # Initialization

# # Iterate through each library and count 
# # how much each element appears
#   for (LIB in LIB_VECTOR) {
#     
#     if (LIB == 'Recurrent') {
#       # Use Recurrent Library
#       ChimLib = ChimeraOutput
#     } else {
#       
#       # For each library create submatrix of only that libraries entries
#       ChimLib = Chimera[ which(LIB == Chimera['LIBRARY']),]
#     }
#     
# #     Repeats = matrix(unlist(
# #       strsplit(as.character(
# #         ChimLib[,'repeatName']),':'))
# #                      ,ncol=3,byrow=T)
#     
#   # Initialize Outpute
#     Hits = nrow(Repeats)
#     hitName = length(repName)
#     hitClass = length(repClass)
#     OutLib = data.frame(Hits) # Hits in Library
#     
#   # TE Class
#     OutLib$'LINE' = length(which(Repeats[,2] == 'LINE'))
#     OutLib$'SINE' = length(which(Repeats[,2] == 'SINE'))
#     OutLib$'LTR' = length(which(Repeats[,2] == 'LTR'))
#     OutLib$'DNA' = length(which(Repeats[,2] == 'DNA'))
#     OutLib$'Other' = length(which(Repeats[,2] == 'Other'))
#     
#   # LTR Family
#     OutLib$'ERV1' = length(which(Repeats[,3] == 'ERV1'))
#     OutLib$'ERVL' = length(which(Repeats[,3] == 'ERVL'))
#     OutLib$'ERVL-MaLR' = length(which(Repeats[,3] == 'ERVL-MaLR'))
#     
#   # THE Family
#     OutLib$'THE' = length(grep("THE*",Repeats[,1]))
#     OutLib$'nonTHE' = OutLib$'ERVL-MaLR' - OutLib$'THE'
#     
#     OutLib$'THE1A' = length(which("THE1A" == Repeats[,1]))
#     OutLib$'THE1B' = length(which("THE1B" == Repeats[,1]))
#     OutLib$'THE1C' = length(which("THE1C" == Repeats[,1]))
#     OutLib$'THE1D' = length(which("THE1D" == Repeats[,1]))
#     OutLib$'THEint' = length(grep("THE.*int",Repeats[,1],perl=T))
#     
#     if (FIRST) {
#       LIB_HITS = OutLib
#       FIRST = F
#     } else {
#       LIB_HITS = rbind(LIB_HITS,OutLib)
#     }}

# Chimera Sense/Antisense Binary Model (presense absence) ----------

ChimIO = data.frame(Type = c('Sense','AntiS','Inter','Cmplx'))
  
for (LIB in as.character(Groups[,1])) {
  COL = colnames(ChimIO)
  Chimera.libary = Chimera[which(Chimera$LIBRARY == LIB),]
  
  Sense = Chimera.libary$RefID[which(Chimera.libary$assXref == 's')]
    Sense = unique(unlist(strsplit(as.character(Sense),split='; ')))
  
  AntiS = Chimera.libary$RefID[which(Chimera.libary$assXref == 'as')]
    AntiS = unique(unlist(strsplit(as.character(AntiS),split='; ')))
  
  Inter = Chimera.libary$RefID[which(Chimera.libary$assXref == 'i')]
    Inter = unique(unlist(strsplit(as.character(Inter),split='; ')))
  
  Cmplx = Chimera.libary$RefID[which(Chimera.libary$assXref == 'c')]
    Cmplx = unique(unlist(strsplit(as.character(Cmplx),split='; ')))

  ChimIO[,LIB] = list(rbind(list(Sense),list(AntiS),list(Inter),list(Cmplx)))
}

# Write Output ------------------------------------------------------

# Write output
write.table(ChimeraOutput,
	file = OUTPUT,
      	quote = F,
	sep = '\t',
	row.names = F,
	col.names = T)

# Write binary output
	save(ChimIO,file=paste(OUTPUT,'IO.Rdata',sep=''))
	
} # End for loop

# 
# write.table(LIB_HITS,
#             file = paste(OUTPUT,".stats.csv", sep=''),
#             quote = F,
#             sep = '\t',
#             row.names = F,
#             col.names = T)

# End of Script
