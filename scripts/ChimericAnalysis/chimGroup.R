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

  if ( length(Match) == 0){
    Assign = '0'
  } else {
    Assign = as.numeric(Groups[Match,3])
  }
  
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


Fishers.Matrix = function(X,Y){
  
  # Fisher's Exact Test (P-value return)
  # X = number of Group A libraries
  # Y = number of Group B libraries
  # returns a (X+1)x(Y+1) matrix to include 0
  # corresponding to the p-value of how many
  # librares in Group A or Group B are chimera
  # positive
  
  VALUE = c() #initialize p-value output vector
  
  for (ITERATION in (sequence(X + 1) - 1) ) {
    
    X_positive = ITERATION # start at zero positive
    X_negative = X - ITERATION # and all are negative
    
    for (ITERATION2 in (sequence(Y+1) - 1) ) {
      
      Y_positive = ITERATION2
      Y_negative = Y - ITERATION2
      
      Xvector = c(X_positive,X_negative)
      Yvector = c(Y_positive,Y_negative)
      
      VALUE = cbind( VALUE, 
                     fisher.test( rbind(c(X_positive,X_negative),
                                        c(Y_positive,Y_negative)  ),
                                  alternative = 'greater')$p.value)
    } }
  
  MAT = matrix(VALUE, nrow = X + 1, byrow = T, ncol = Y + 1)
  
  return(MAT)
  
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
Go = T # Keep going through libraries

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

# Clean Upmas.sou
#  \    @ /    
#   \ @  /
#    ----
 rm( ChimFiltered, Count, AXR, Filter, Filter_Cancer, Filter_Normal,
      Go, ID, IDdup, Matches, RefID)
  gc()

# Library Statistics ------------------------------------------------

# Load Repeat Masker Database
  load(RMDATA_PATH)

# Parse RMDB to RepeatID Vector
# chrX:<start>
  RMDB_id <- unique( paste( RMDB$genoName, RMDB$genoStart, sep=":") )


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


# Chimera Sense/Antisense Binary Model (presense absence) ----------
# Protein Coding Genes and their interacting transcript orientation
ChimProtIO = data.frame(Type = c('Sense','AntiS','Inter','Cmplx'))

for (LIB in as.character(Groups[,1])) {
  COL = colnames(ChimProtIO)
  Chimera.libary = Chimera[which(Chimera$LIBRARY == LIB),]
  
  Sense = Chimera.libary$RefID[which(Chimera.libary$assXref == 's')]
    Sense = unique(unlist(strsplit(as.character(Sense),split='; ')))
  
  AntiS = Chimera.libary$RefID[which(Chimera.libary$assXref == 'as')]
    AntiS = unique(unlist(strsplit(as.character(AntiS),split='; ')))
  
  Inter = Chimera.libary$RefID[which(Chimera.libary$assXref == 'i')]
    Inter = unique(unlist(strsplit(as.character(Inter),split='; ')))
  
  Cmplx = Chimera.libary$RefID[which(Chimera.libary$assXref == 'c')]
    Cmplx = unique(unlist(strsplit(as.character(Cmplx),split='; ')))

  ChimProtIO[,LIB] = list(rbind(list(Sense),list(AntiS),list(Inter),list(Cmplx)))
}

# Chimera Presense/Absense Vector
# Protein Coding Genes and their interacting transcript orientation
if (RUN == 1){
  #print(" Running standard analysis" )

  # Chimeric Identifiers
  ChimID = paste(as.character(RMDB$genoName),':',as.character(RMDB$genoStart),sep='')
  
  # Chimeric Bionary (IO) dataframe based on RMDB order
  ChimIO = data.frame(ChimID)
  
  for (LIB in as.character(Groups[,1])) {
    
    # Initialize a FALSE vector for all RMDB
    entryIO = vector('logical', length = nrow(RMDB))
    
    # Which Chimera in Library are present
    # Assign those REFID true in EntryIO
    Chimera.library = Chimera[which(Chimera$LIBRARY == LIB),]$RepeatID
    
    entryIO[which(ChimID %in% Chimera.library)] = TRUE
    
    # Pass this onto the ChimIO logical matrix
    ChimIO[[LIB]] = entryIO
  }
  
  # Recurrent and Specific Chimera as defiend by LIONS
  # Initialize a FALSE vector for all RMDB
  entryIO = vector('logical', length = nrow(RMDB))
  
  # Which Chimera in Library are present
  # Assign those REFID true in EntryIO
  Chimera.library = as.character(ChimeraOutput$RepeatID)
  
  entryIO[which(ChimID %in% Chimera.library)] = TRUE
  
  # Append recurrent and specific entries
  ChimIO[['Recurrent']] = entryIO
  
} else {
  #print(" Running inverse analysis" )
  
  # Inverse of Recurrent and Specific Chimera as defiend by LIONS
  # Initialize a FALSE vector for all RMDB
  entryIO = vector('logical', length = nrow(RMDB))
  
  # Which Chimera in Library are present
  # Assign those REFID true in EntryIO
  Chimera.library = as.character(ChimeraOutput$RepeatID)
  
  entryIO[which(ChimID %in% Chimera.library)] = TRUE
  
  # Append recurrent and specific entries
  ChimIO[['invRecurrent']] = entryIO
}

# Write Output ------------------------------------------------------

# Write output
write.table(ChimeraOutput,
	file = OUTPUT,
      	quote = F,
	sep = '\t',
	row.names = F,
	col.names = T)

} # End for loop

# Write binary output
save(ChimIO,file=paste(pNAME,'IO.Rdata',sep=''))

# End of Script
