# chimNSort.R - <ER_Interaction> <Brunswick Directory>
#
# From Standard Input of the matches.sh script: read
# Exon / Repeat Interaction Type:

# Libraries
  library(neuralnet)
  library(GenomicRanges)

# ==============================================================================
# Control Panel ================================================================
# ==============================================================================

# Scripting Inputs
	#STDIN = commandArgs(trailingOnly = TRUE)
	#INTERACTION = STDIN[1]
	#BRUNSWICK = STDIN[2]
		#cd $BRUNSWICK

# Manual Inputs
  # Set type of Interaction manually
    #INTERACTION = 'UpEdge'

# TRAINING or DATA EXPLORATION PROTOCOL
  DATA_SET = 'DATA_train.Rdata'
  TRAIN = T # run training [T] or run visualization [F]
    ITER = 5e6 # Number of iterations per model [1e5]
    THRESH = 0.00001 # Threshold for model cutoff [0.001]
    REPS = 1 # number of models to train per file
    SINGLE = T # Single Iteration or Infinite loop of training?
  
  RAND_START = F # Random initial weights (T) or known start model (F)
    ANN_TRAINED = 'ANN_Up_rOOjTx.Rdata'

  if (TRAIN) {
# Training Initialization Parameters      
    
    if (RAND_START) {
      # Random Start
        load(DATA_SET)
        
    } else {
      # Fixed initial weights
        load(DATA_SET)
        load(ANN_TRAINED) # load a trained neural network
        WEIGHTS = nn$weights
        rm(nn) }
    
  } else {    
# Data Exploration Initialization
    
    library(ggplot2)
    # Load a trained neural network
      load(ANN_TRAINED)
      
    # Load a Master Value Table from Chimeric Pipeline
      # can be ALL data, TRAIN data or TEST data
      load(DATA_SET)
}


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
# TRAINING ANN =================================================================
# ==============================================================================
    
COUNT = 1

RUNTRAINING = TRAIN # Set while command to control panel setting

while (RUNTRAINING) {
  
  if (RAND_START) {
  
  # Random Start Weights
    nn = neuralnet(
      data = DATA.ue,
      formula = TruePos ~ Total + exonRankInTranscript + RepeatRank + ThreadRatio + ExonMax + UpExonMax + RepeatMaxCoverage,
      hidden = 7,
      err.fct = "ce",
      linear.output = F,
      algorithm = 'rprop+',
      stepmax = ITER,
      threshold = THRESH,
      rep = REPS)
  
  } else {
    
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
      startweights = WEIGHTS)
    
  }
  
# Save Progress -----------            

  # Unique Run ID
    runID = makeRandomString(1,6)   # Unique Run Identifer
  
  # Summary Statistics
    NUM_ANN = length(nn$result.matrix[1,])
    MIN_ERROR = min(nn$result.matrix[1,])
    MIN_I = which(as.vector( MIN_ERROR == nn$result.matrix[1,]))
    MIN_STEPS = nn$result.matrix[3,MIN_I]
    
  nn.summary = paste('runID:',runID,'  Models_trained:',NUM_ANN,
                     '  Minimum_Error:', MIN_ERROR, '  in:', MIN_STEPS, 'steps')
  
  # Write output
  write.table(nn.summary,paste0(INTERACTION,'/Summary.txt'),append=T,quote=F,row.names=F,col.names=F)
                   
  OUTPUT = paste('ANN_',INTERACTION,'_',runID,'.Rdata', sep = '') # Output Name
    OUTPUT.data = c('nn','runID')
  
  save(nn, file=OUTPUT)
  
print(paste('Finished Run:',COUNT))
  COUNT = COUNT + 1
  
  # Exit loop or continue with a new training set?
  if (SINGLE) {
    RUNTRAINING = F
  }
}


# ==============================================================================
# Data Exploration =============================================================
# ==============================================================================
    
if (TRAIN == F) {  

  library(ggplot2)
#   #Test Neuronal Network
#     NN = data.frame(nn$response)
#     NN = cbind(NN,as.numeric(unlist(nn$net.result[[1]])))
#     colnames(NN) = c('Response','ANN')
    
  # # Plot an input variable between the two groups
  #   PLOT = ggplot(data = DATA.ue, aes(x = TruePos, y = RepeatRank))
  #     PLOT = PLOT + geom_violin(aes(fill=factor(TruePos)), alpha = 0.5)
  #     PLOT
  
  # ANN
    DATA.model = data.frame(DATA.ue[,c('Total','exonRankInTranscript','RepeatRank','ThreadRatio',
                                  'ExonMax','UpExonMax','RepeatMaxCoverage')])
    nn.df = compute(nn,covariate=DATA.model,1)
  
    nn.df = data.frame(cbind(unlist(DATA.ue$TruePos), unlist(nn.df$net.result)))
      colnames(nn.df) = c('response','net.result')
  
  #AUC plot
  TruePos = which(nn.df$response == 1)
    Pos = length(TruePos)
  
  TrueNeg = which(nn.df$response !=1)
    Neg = length(TrueNeg)
  
  # Initialize Output
    bins = 10000
    OUTPUT = t(matrix(c(0,0,0,0)))
  
  for (N in 0:bins){
    # Cycle through net.result = 0.0 to 1.0
      M <<- N / bins
      # Number of called positive/Negative at this M
        tPos <<- length(which(nn.df$net.result[TruePos] >= M))
        fNeg <<- length(which(nn.df$net.result[TruePos] < M))
      
        tNeg <<- length(which(nn.df$net.result[TrueNeg] < M))
        fPos <<- length(which(nn.df$net.result[TrueNeg] >= M))
      
          sumPos <<- tPos + fPos
          sumNeg <<- tNeg + fNeg
      
      Cutoff <<- M
      TPR <<- tPos / Pos
      FPR <<- fPos / Neg
      FDR <<- fPos / (tPos + fPos)
    
      OUTPUT <<- rbind(OUTPUT, c(Cutoff,TPR,FPR,FDR))
  }
  
    OUTPUT = data.frame(OUTPUT)
      colnames(OUTPUT) = c('Cutoff','TPR','FPR','FDR')
    
  PLOT = ggplot(data = OUTPUT, aes(x = FPR, y = TPR))
    PLOT = PLOT + ggtitle('Chimera ROC - "-" Classification')
    PLOT = PLOT + geom_line(aes(colour='red'))
    PLOT = PLOT + geom_abline(intercept = 0, slope = 1)
    PLOT = PLOT + theme(legend.position = "none")
    PLOT
  
    PLOT = ggplot(data = nn.df, aes(x = net.result))
      PLOT = PLOT + geom_density(aes(fill=factor(response)), alpha = 0.5)
      PLOT
  
  ## Parameter Plots
  # par(mfrow=c(2,2))
  #   gwplot(nn,selected.covariate="Total", min=-10, max=10)
  #   gwplot(nn,selected.covariate="exonRankInTranscript", min=-10, max=10)
  #   gwplot(nn,selected.covariate="RepeatRank", min=-10, max=10)
  #   gwplot(nn,selected.covariate="ThreadRatio",min=-10, max=10)
  
} # end visualization
