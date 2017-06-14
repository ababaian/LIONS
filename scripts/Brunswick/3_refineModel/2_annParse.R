# annParse.R
#
# From Standard Input of the matches.sh script: read
# Exon / Repeat Interaction Type:
#
# This script reduces the size of the ANN_*.Rdata models
# by 'retraining' it against a small data set.
# The model doesn't change but some unneccesary data is removed
# 
# Saves a master output file with all values needed for LIONS useage
# --> ANN_MODEL_######.Rdata

# Libraries
library(neuralnet)
library(GenomicRanges)

# ==============================================================================
# Control Panel ================================================================
# ==============================================================================
# Unique Model Identifier
DATE = '_160919'

# Best models for input
MODELS = c('ANN_Up_OuEDtG.Rdata',
           'ANN_UpEdge_mqm5M0.Rdata',
           'ANN_EInside_qblwie.Rdata')

CUTOFFS = c('Cutoff_ROC01_Up_OuEDtG_1_.Rdata',
            'Cutoff_ROC01_UpEdge_mqm5M0_2_.Rdata',
            'Cutoff_ROC01_EInside_qblwie_1_.Rdata')

REP_N = c(1,2,1) # Reps within the above models to be used in final model

count = 1

for (M in MODELS){
  
# Scripting Inputs
  ANN_TRAINED = M
  INTERACTION = unlist( strsplit(ANN_TRAINED,split = '_'))[2]
  TRAINED_REP = REP_N[count] # Model number in ANN_Trained to use (1 is default)
    count = count + 1

# TRAINING PARAMETERS for Reducing ANN file sizes
  DATA_SET = 'DATA_train.Rdata'
  TRAIN = T # run training [T] or run visualization [F]
  ITER = 2 # Number of iterations per model [1e5]
  THRESH = 100 # Threshold for model cutoff [0.001]
  REPS = 1 # number of models to train per file
  SINGLE = T # Single Iteration or Infinite loop of training?
  RAND_START = F # Random initial weights (T) or known start model (F)

# Fixed initial weights
    load(DATA_SET)
    load(ANN_TRAINED) # load a trained neural network
    WEIGHTS = nn$weights
    rm(nn)

# Reduce size of DATA set
# to make the model files smaller.
  miniDATA = DATA[1:1000,]
  DATA = miniDATA


# Functions --------------------------------------------------------------------

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

# Load and Initialize ----------------------------------------------------------

# Transform the Data for the ANN
# Log and Scale Data to [0,1]
DATA.ue = DATA[which(DATA$ER_Interaction == INTERACTION),]

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

# ==============================================================================
# Reduce ANN Size ==============================================================
# ==============================================================================

  # Known Start Weights
  nn = neuralnet(
    data = DATA.ue,
    formula = TruePos ~ Total + exonRankInTranscript + RepeatRank + ThreadRatio + ExonMax + UpExonMax + RepeatMaxCoverage,
    hidden = 7,
    err.fct = "ce",
    linear.output = F,
    algorithm = 'rprop+',
    stepmax = ITER,
    threshold = THRESH,
    rep = REPS,
    startweights = WEIGHTS[TRAINED_REP])
    

# Save Progress -----------

if (INTERACTION == 'Up'){
  nn_up = nn
} else if (INTERACTION == 'UpEdge'){
  nn_upedge = nn
} else if (INTERACTION == 'EInside'){
  nn_einside = nn
}

# # Write output
# 
# # Unique Run ID
#   runID = makeRandomString(1,6)   # Unique Run Identifer
# 
# # Summary Statistics
#   NUM_ANN = length(nn$result.matrix[1,])
#   MIN_ERROR = min(nn$result.matrix[1,])
#   MIN_I = which(as.vector( MIN_ERROR == nn$result.matrix[1,]))
#   MIN_STEPS = nn$result.matrix[3,MIN_I]
#   
#   nn.summary = paste('runID:',runID,' Interaction:',INTERACTION,'  Models_trained:',NUM_ANN,
#                      '  Minimum_Error:', MIN_ERROR, '  in:', MIN_STEPS, 'steps')
# 
# # Write output
# #  write.table(nn.summary,paste0('Model.txt'),append=T,quote=F,row.names=F,col.names=F)
#   
#   OUTPUT = paste('ANN_', INTERACTION, DATE, '.Rdata', sep = '') # Output Name
#   OUTPUT.data = c('nn','runID')
#   
#   save(nn, file=OUTPUT)

}

# ==============================================================================
# Single Model Output ==========================================================
# ==============================================================================


# # Clean-up workspace
# # rm(list = ls())
# rm(DATA, DATA.ue, miniDATA, ANN_TRAINED, count, DATA_SET)
# rm(ITER, M, MIN_ERROR, MIN_I, MIN_STEPS, nn, nn.summary, NUM_ANN, oneHold)
# rm(OUTPUT, OUTPUT.data, RAND_START, REPS, runID, SINGLE, THRESH, TRAIN, TRAINED_REP)
# rm(zeroHold, WEIGHTS, INTERACTION)
# 
# # ANN Models
#   M_up = paste('ANN_Up', DATE, '.Rdata', sep = '')
#   M_upedge = paste('ANN_UpEdge', DATE, '.Rdata', sep = '')
#   M_einside = paste('ANN_EInside', DATE, '.Rdata', sep = '')

# Load ROC Cutoff files
  C_up = CUTOFFS[1]
  C_upedge = CUTOFFS[2]
  C_einside = CUTOFFS[3]

# Load ANN Models into independent objects
  load(C_up)
    cut_up = cut01
    upCUT = cut_up$ROC01$Global$optimal.cutoff$cutoff
    

  load(C_upedge)
    cut_upedge = cut01
    upedgeCUT = cut_upedge$ROC01$Global$optimal.cutoff$cutoff
    
  load(C_einside)
    cut_einside = cut01
    einsideCUT = cut_einside$ROC01$Global$optimal.cutoff$cutoff
   
    
save(MODELS, nn_up, nn_upedge, nn_einside, cut_up, cut_upedge, cut_einside, upCUT, upedgeCUT, einsideCUT, file = paste(sep='','ANN_Model_',DATE,'.Rdata')) 
    
    
