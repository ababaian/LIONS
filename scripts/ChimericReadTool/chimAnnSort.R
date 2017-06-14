# chimAnnSort.R
#
# chimAnnSort.R <ChimericResults.pc.lcsv> <FilteredOutputName> <ANN Model Path>
#
#
# ---------------------------------------------------------
# Reads ChimericReadTool.sh output
#     Parse Rdata
#
#

# Libraries
  library(neuralnet)

# ==============================================================================
# Control Panel ================================================================
# ==============================================================================

 STDIN = commandArgs(trailingOnly = TRUE)
 
 INPUT = STDIN[1] # Chimera_Results
   INPUTNAME = unlist(strsplit(as.character(INPUT), split = "\\."))[1]
 
 OUTPUT = STDIN[2] # Chimeric Output
 
 # ANN Model Paths
 MODEL_PATH = STDIN[3]

# Manual Entry Parameters ----------------------

#INPUT = 'A05254r.lcsv'
#  INPUTNAME = 'A05254r'

#OUTPUT = 'A05254r.ann.lion'

#MODEL_PATH = 'ANN_Model_160919.Rdata'
  # list of the types of models used.
  # this is 'frail' in the sense it should be in the .Rdata
  # file in the future. If another ANN model is being used
  # for another interaction it could break the script.

  
modelTypes = c('up', 'upedge', 'einside')
interactionNames = c('Up', 'UpEdge', 'EInside')
scTHREAD = 10 # No Thread Ratio cut-off is set, simply write '10' when INF

# IMPORT ======================================================================

# Read output from ChimericReadTool.sh
ChimeraTable = read.csv(file=INPUT,
                        header=T,
                        sep='\t')

# Set Initial Chimera Table to DATA for ANN-classification
  DATA = ChimeraTable 

# Load Model Data
  load(MODEL_PATH)
  
# Functions ====================================================================

normalize = function(x)
{
  # Normalize a vector to the range [0,1]
  # Remove edge cases
    x[is.nan(x)] = 0
    x[is.na(x)] = 0
    x[is.infinite(x)] = 0
  
    x = x - min(x)
    x = x/max(x)
    
    return(x)
}

makeRandomString <- function(n=1, lenght=12)
{
  randomString <- c(1:n)                  # initialize vector
  for (i in 1:n)
  {
    randomString[i] <- paste(sample(c(0:9, letters, LETTERS),
                                    lenght, replace=TRUE),
                             collapse="")
  }
  return(randomString)
}

# CALCULATED VALUES ===========================================================
# For ChimeraTable and it's output

# CONTRIBUTION (via Max)
  # = Repeat(max reads) / Exon(max reads)
  # Relative expression of Repeat to interacting Exon

  ChimeraTable$Contribution = (ChimeraTable$RepeatMaxCoverage /
                               ChimeraTable$ExonMax )

# UpstreamCoverageRatio (UpCov)
  # = Repeat (max) / UpstreamRepeat (max)
  # Look upstream of the repeat to see if there is expression
  
  ChimeraTable$UpCov = (ChimeraTable$RepeatMaxCoverage /
                          ChimeraTable$UpstreamRepeatMaxCoverage)

# Upstream Exon Experssion
  # If Exon 1, set value to Inf
  # = Exon (RPKM) / UpstreamExon (RPKM)
  
  ChimeraTable$UpExonRatio = (ChimeraTable$ExonRPKM /
                                ChimeraTable$UpExonRPKM)
  
  EX1 = which( ChimeraTable[,'exonRankInTranscript'] == 1 )
  ChimeraTable[EX1,'UpExonRatio'] = Inf

# Thread Ratio
  # DownThread / UpThread
  # If UpThread = 0, set value to cutoff
  ChimeraTable$ThreadRatio = (ChimeraTable$DownThread /
                              ChimeraTable$UpThread)
  
  UpThreadZero = which(ChimeraTable$UpThread == 0)
  ChimeraTable[UpThreadZero,'ThreadRatio'] = scTHREAD

# Multi-Exonic Transcripts
  MultiEx = which( ChimeraTable[,'ExonInGene'] == 2)  

# ANN - Normalize and Classify ========================================

# For each type of classification, there is a seperate ANN Model
# for classification. All output pipes into one vector which is saved.

annClassification = rep(-1, nrow(DATA))
  # {1,0} = classified, where 1 = Chimeric, 0 = Non-Chimeric
  # -1 = non-classified
annResult = rep(-1, nrow(DATA))

# Iterate through each interaction Type and load the respective model
for (N in 1:length(modelTypes)){
  
  # Iteration Parameters
  INTERACTION = modelTypes[N]       # Interaction Type Vector
  intNAME     = interactionNames[N] # Interaction Names
  
  MODEL = eval(parse(text=paste0('nn_',INTERACTION)))  # Respective ANN Model
 
  CUTOFF = eval(parse(text=paste0( INTERACTION, 'CUT'))) # Cutoff Threshold for this model
    
  # DATA TRANSFORMATION ---------------------------------------
  # Log and Scale Data to [0,1]
  DATA.ue = DATA[which(DATA$ER_Interaction == intNAME),]
  
  # Exon Rank in Transcript
    DATA.ue$exonRankInTranscript[ DATA.ue$exonRankInTranscript == 1] = 0
    DATA.ue$exonRankInTranscript[ DATA.ue$exonRankInTranscript > 1] = 1
  
  # Repeat Rank in Transcript
    DATA.ue$RepeatRank[ DATA.ue$RepeatRank == -1] = 0
    DATA.ue$RepeatRank[ DATA.ue$RepeatRank > 0] = 1
  
  # Thread Ratio
    DATA.ue$ThreadRatio = DATA.ue$DownThread/DATA.ue$UpThread
      zeroHold = is.nan(DATA.ue$ThreadRatio) # 0/0 set to min
      oneHold = is.infinite(DATA.ue$ThreadRatio) # x/0 set to max
  
    DATA.ue$ThreadRatio = log10(DATA.ue$ThreadRatio) 
    DATA.ue$ThreadRatio = normalize(DATA.ue$ThreadRatio)
    DATA.ue$ThreadRatio[zeroHold] = 0
    DATA.ue$ThreadRatio[oneHold] = 1
  
  # Total
    DATA.ue$Total = log10(DATA.ue$Total)
    DATA.ue$Total = normalize(DATA.ue$Total)
  
  # Read Types
    DATA.ue$ER = normalize(log10(DATA.ue$ER))
    DATA.ue$DR = normalize(log10(DATA.ue$DR))
    DATA.ue$DE = normalize(log10(DATA.ue$DE))
    DATA.ue$DD = normalize(log10(DATA.ue$DD))
  
  # Exon Max Coverage
    DATA.ue$ExonMax = log10(DATA.ue$ExonMax)
    DATA.ue$ExonMax[is.infinite(DATA.ue$ExonMax)] = 0
    DATA.ue$ExonMax = normalize(DATA.ue$ExonMax)
  
  # UpExon Max Coverage
    DATA.ue$UpExonMax = log10(DATA.ue$UpExonMax)
    DATA.ue$UpExonMax[is.infinite(DATA.ue$UpExonMax)] = 0
    DATA.ue$UpExonMax = normalize(DATA.ue$UpExonMax)
  
  # Repeat Max Coverage
    DATA.ue$RepeatMaxCoverage = log10(DATA.ue$RepeatMaxCoverage )
    DATA.ue$RepeatMaxCoverage [is.infinite(DATA.ue$RepeatMaxCoverage )] = 0
    DATA.ue$RepeatMaxCoverage  = normalize(DATA.ue$RepeatMaxCoverage )
  
  # Exon Max Coverage
    DATA.ue$UpstreamRepeatMaxCoverage = log10(DATA.ue$UpstreamRepeatMaxCoverage)
    DATA.ue$UpstreamRepeatMaxCoverage[is.infinite(DATA.ue$UpstreamRepeatMaxCoverage)] = 0
    DATA.ue$UpstreamRepeatMaxCoverage = normalize(DATA.ue$UpstreamRepeatMaxCoverage)

  # Classification ---------------------------------------
  
  # Load ANN Model
    WEIGHTS = MODEL$weights
  
  DATA.model = data.frame(DATA.ue[,c('Total','exonRankInTranscript','RepeatRank','ThreadRatio',
                                     'ExonMax','UpExonMax','RepeatMaxCoverage')])
    
    nn.df = compute(MODEL, covariate=DATA.model, 1)
      # It would be a good idea to locally define 'compute' function
      # so that the dependency on library(neuralnet) can be removed
      # from the main LIONS pipeline. (still would be required for training)
  
  Classification = rep(-1,nrow(nn.df$net.result))
    Classification[which(nn.df$net.result >= CUTOFF)] = 1
    Classification[which(nn.df$net.result < CUTOFF)] = 0
  
  
  # For each intNAME in annClassification;
    annClassification[which(DATA$ER_Interaction == intNAME)] = Classification
    annResult[which(DATA$ER_Interaction == intNAME)] = nn.df$net.result
  
}

# CHIMERA FILTRATION ==========================================================
# (via pre-trained ANN Models)

ChimeraOut = ChimeraTable[which(annClassification == 1),]

  print(paste('     Output Entries   : ', nrow(ChimeraOut)))

# Start of Repeat (For comparisons)
  ChimeraOut$RepeatID = paste('chr',ChimeraOut$Chromosome,':',ChimeraOut$RStart,sep='')

# Sort Rows
  SORT=4 # Coordinates
  ChimeraOut = ChimeraOut[order(ChimeraOut[,SORT]),]


# Add input to last column
#  (for comparing multiple lists)

# LIBRARY = unlist(strsplit(INPUT,split = '/'))[1]
  LIBRARY = INPUTNAME
  ChimeraOut = cbind(ChimeraOut, LIBRARY)

# To each entry; add ann Classification result
  ChimeraOut$annResult = annResult[which(annClassification == 1)]


# PARSE and OUTPUT ============================================================

# Sort by ANN classification
  ChimeraOut = ChimeraOut[order(ChimeraOut$annResult, decreasing = T),]
  
# Write output CSV
write.table(ChimeraOut,
            file = OUTPUT,
            quote = F,
            sep = '\t',
            row.names = F,
            col.names = T)

# End of Script :D
