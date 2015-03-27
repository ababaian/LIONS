# rmResource.R <RepeatMasker.ucsc>
#
# Input standard RepeatMasker table and
# process it into a table of element-wise statistics
# 

# INPUT =============================================================

  # Input UCSC RepeatMasker Table (also rmdb files should work)
  # Downloaded from
  # http://genome.ucsc.edu/cgi-bin/hgTables
  # RepeatMasker hg19 - all genome
  # Check README for more details

  STDIN = commandArgs(trailingOnly = TRUE)

  RMDB_PATH = STDIN[1] # Input LIONS project file
  # Columns of input <rm>.ucsc file
  # See LIONS README for information on file columns

  # Output Rdata files (binary file)
  OUTPUT = 'rmStats.Rdata'

  # Repeat Masker Database
  RMDB = read.csv(RMDB_PATH,header=T, sep = '\t')[,c(-1,-17)]

# FUNCTIONS =========================================================

  # SearchWhich
  SearchWhich = function(Key, Query) {
    # Use which to lookup a Key string
    # in Query vector for lapply
    return = which(Key == Query) }

  # ListSum
  ListSum = function(Key,Query) {
    # For a list <key> look-up those elements
    # in Query numeric vector and take their sum
    return = sum(Query[unlist(Key)]) }

  ListMean = function(Key,Query) {
    # For a list look-up Query and return
    # a list of c(Mean,Standard Deviation)
    return = paste(   mean( Query[unlist(Key)] ) ,
                      sd(   Query[unlist(Key)] ) ,
                      sep = ';')
  }

# Chimeric Read Search CSV File =====================================

# RepeatMasker.Rdata ================================================
# table Name/Class/Family statistics for analysis in West Lion

# RepeatNames 
  Name = sort(unique(RMDB$repName))
# RepeatFamily
  Family = sort(unique(RMDB$repFamily))
# RepeatClass
  Class = sort(unique(RMDB$repClass))

# Element Size calculation (bp)
  RMDB$size = RMDB$genoEnd - RMDB$genoStart

# RepeatStats
  # statistical analysis of each repeat across the genome
  rmStats = data.frame(Name)
    famStats = data.frame(Family)
    claStats = data.frame(Class)
  
  # Repeat Family
  rmStats$'Family' = RMDB$repFamily[unlist(lapply(Name,match,table=RMDB$repName))]

  # Repeat Class
  rmStats$'Class' = RMDB$repClass[unlist(lapply(Name,match,table=RMDB$repName))]
    famStats$'Class' = RMDB$repClass[unlist(lapply(Family,match,table=RMDB$repFamily))]

# For each repName find the matching entries in RMDB
  nameMatch = lapply(rmStats$Name, FUN=SearchWhich, Query = RMDB$repName)
    famMatch = lapply(famStats$Family, FUN=SearchWhich, Query = RMDB$repFamily)
    claMatch = lapply(claStats$Class, FUN=SearchWhich, Query = RMDB$repClass)

# CoverageStats
  # Number of distinct elements in the Genome
  rmStats$'Number' = unlist(lapply(nameMatch,length))
    famStats$'Number' = unlist(lapply(famMatch,length))
    claStats$'Number' = unlist(lapply(claMatch,length))
  
  # The total bases covered by each distinct element in the genome
  rmStats$'Bases' = unlist(lapply(nameMatch,ListSum,Query=RMDB$size))
    famStats$'Bases' = unlist(lapply(famMatch,ListSum,Query=RMDB$size))
    claStats$'Bases' = unlist(lapply(claMatch,ListSum,Query=RMDB$size))

  # The average size of an element and it's standard deviation
  rmStats$'Size_mean' = unlist(lapply(nameMatch,ListMean,Query=RMDB$size))

  # The average divergence of an element and it's standard deviation
  rmStats$'Div_mean' = unlist(lapply(nameMatch,ListMean,Query=RMDB$milliDiv))

  # The average deletion of an element and it's standard deviation
  rmStats$'Del_mean' = unlist(lapply(nameMatch,ListMean,Query=RMDB$milliDel))

# ELEMENT/ CLASS REPRESENTATION =====================================

# Create a data.frame of only Transposable Element data
# and create Class/Family/Repeat -wise representation values

# Transposable Element List
  TE = c("LINE", "SINE", "LTR", "DNA", "Other")

# TE Data Frame
  teStats = rmStats[ which( rmStats$Class %in% TE ),]
    teSum = sum(teStats$Bases) # Total TE bases

# TE Representation (% of Total)
  # There are 3 levels of organization of Transposable Elements
  # in Repeat Masker
  # Transposable Element Class > Family > Element Name
  #
  # Representation is how much of one catagory is composed of a lower
  # catagory of organization
  # i.e. totalTE = Percent of all TE made up by TE Class ALPHA
  # i.e. totalFAM = Percent of all TE made up by TE Family Alphie
  # i.e. classFAM = Percent of all class ALPHA made up by family Alphie
  # i.e. famNAME = Percent of ERV-MaLR Family made up by THE1B elements

# Count Bases at each level
  teStats$totalBP = sum(teStats$Bases) # Total Bases covered by all TE

  teStats$classBP = unlist(lapply(teStats$Class, FUN = function(CLASS){
    claStats$Bases[which( CLASS == claStats$Class )] }))

  teStats$famBP = unlist(lapply(teStats$Family, FUN = function(FAMILY){
    famStats$Bases[which( FAMILY == famStats$Family )] }))
  
  #teStats$nameBP = teStats$Bases

# Representation Calculation
  teStats$classTOTAL = 100 * (teStats$classBP / teStats$totalBP)
  teStats$familyTOTAL = 100 * (teStats$famBP / teStats$totalBP)
  teStats$teTOTAL = 100 * (teStats$Bases / teStats$totalBP)

  teStats$familyCLASS = 100 * (teStats$famBP / teStats$classBP)
  #teStats$teCLASS = 100 * (teStats$Bases / teStats$classBP  ## Not used

  teStats$teFAMILY = 100 * (teStats$Bases / teStats$famBP)

  teStats = subset(teStats, select = -c(totalBP,classBP,famBP) )


# OUTPUT ============================================================
  # Save output variables as binary
  save(RMDB,rmStats,famStats,claStats,teStats,file=OUTPUT)
  
  
