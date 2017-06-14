# annROC.R - <DATA_set.Rdata> <ANN_Model_ID.Rdata>
#
# After training Neural Networks with chimNeuro / matches
# 
# This script plots the ANN model ROC curves (vs. DATA_SET)
# and calcualtes optimal cut-off thresholds
# 
# Saves a plot and cut-off file

# Libraries
  library(neuralnet)
  library(GenomicRanges)
  library(ggplot2)
  library(OptimalCutpoints)
  

# ==============================================================================
# Control Panel ================================================================
# ==============================================================================

# Scripting Inputs
	STDIN = commandArgs(trailingOnly = TRUE)
  DATA_SET = STDIN[1]
  ANN_TRAINED = STDIN[2]

# Manual Inputs  
  #DATA_SET = 'DATA_test.Rdata'
  #ANN_TRAINED = 'ANN_Up_OuEDtG.Rdata'
  #ANN_TRAINED = 'ANN_UpEdge_7DZpcf.Rdata'
  #ANN_TRAINED = 'ANN_UpEdge_mqm5M0.Rdata'
  #ANN_TRAINED = 'ANN_EInside_qblwie.Rdata'
  #ANN_TRAINED = 'ANN_RInside_

  CUTOFF = 'ROC01' # Cutoff type (From optimal.cutoffs)
  
# Parse Inputs
  INTERACTION = unlist(strsplit(ANN_TRAINED,split = "_"))[2]
  RUNID = unlist(strsplit(ANN_TRAINED,split = "_"))[3]
    RUNID =  unlist(strsplit(RUNID, split = '.Rdata'))

# Load Data sets and ANN Trained
  load(DATA_SET)
  load(ANN_TRAINED) # load a trained neural network
  WEIGHTS = nn$weights
  REPS = length(WEIGHTS) # Number of models trained in nn file
  
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
# Data Exploration =============================================================
# ==============================================================================
    
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
      
for (nREP in 1:REPS){
  
  # ANN
    nn.df = compute(nn, covariate=DATA.model, rep = nREP)
  
    nn.df = data.frame(cbind(unlist(DATA.ue$TruePos), unlist(nn.df$net.result)))
      colnames(nn.df) = c('response','net.result')
  
  #AUC plot
  TruePos = which(nn.df$response == 1)
    Pos = length(TruePos)
  
  TrueNeg = which(nn.df$response !=1)
    Neg = length(TrueNeg)
  
  # Initialize Output
    bins = 10000
    OUTPUT = t(matrix(c(0,0,0,0,0,0,0,0)))
  
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
    
      OUTPUT <<- rbind(OUTPUT, c(Cutoff,TPR,FPR,FDR,tPos,fNeg,tNeg,fPos))
  }
  
    OUTPUT = data.frame(OUTPUT)
      colnames(OUTPUT) = c('Cutoff','TPR','FPR','FDR','tPos','fNeg','tNeg','fPos')

    cut01 = optimal.cutpoints(data = nn.df, X = "net.result", status = "response", tag.healthy = 0,
                              methods = CUTOFF) # Minimum Eucledian distance to (0,1)
    
      
png(filename = paste(sep = '_',"ROC", INTERACTION, RUNID, nREP,".png"), width = 800, height = 800)
  PLOT = ggplot(data = OUTPUT, aes(x = FPR, y = TPR))
    PLOT = PLOT + theme_bw(base_size = 18) 
    PLOT = PLOT + ggtitle(paste('Chimera ROC -',INTERACTION,'Classification'))
    PLOT = PLOT + xlab("1 - Specificity")
    PLOT = PLOT + ylab("Sensitivity")
    PLOT = PLOT + geom_line(aes(colour='red'))
    # X = Y Line
    PLOT = PLOT + geom_abline(intercept = 0, slope = 1)
    # Optimal Cutoff Based
    PLOT = PLOT + geom_vline(xintercept = as.numeric(1-round( cut01[[CUTOFF]]$Global$optimal.cutoff$Sp,4)) , colour = 'green', linetype = "longdash")
    # Number of TP in set and AUC
    PLOT = PLOT + annotate("text", x = 0.90, y = 0.05, size = 8, label = paste("N TP = ",length(which( DATA.ue$TruePos == 1) )))
    PLOT = PLOT + annotate("text", x = 0.90, y = 0.1, size = 8, label = paste("AUC = ", round(cut01[[CUTOFF]]$Global$measures.acc$AUC[1],4) ))
    # Sensitivty and Specificity at Cutoff
    PLOT = PLOT + annotate("text", x = 0.90, y = 0.25, size = 8, label = paste("SE = ",round( cut01[[CUTOFF]]$Global$optimal.cutoff$Se, 4 )))
    PLOT = PLOT + annotate("text", x = 0.90, y = 0.3, size = 8, label = paste("SP = ", round( cut01[[CUTOFF]]$Global$optimal.cutoff$Sp, 4 )))
    PLOT = PLOT + annotate("text", x = 0.90, y = 0.20, size = 8, label = paste("Cutoff = ", round(cut01[[CUTOFF]]$Global$optimal.cutoff$cutoff,4)))
    
    PLOT = PLOT + theme(legend.position = "none")
    print(PLOT)
    
dev.off()

# Save cutoff value file
# Cutoff_METHOD_INTERACTION_nREP_RUNID.Rdata
save(cut01, file = paste(sep = '_', "Cutoff",CUTOFF,INTERACTION,RUNID,nREP,'.Rdata'))

}

      
# Other / Older Code -----------------------------------------------           
# cut01
# cut01[[CUTOFF]]$Global$optimal.cutoff$cutoff

# # False Discovery Rate over Response Curve
#   PLOT = ggplot(data = OUTPUT, aes(x = Cutoff, y = FDR))
#   PLOT = PLOT + theme_grey(base_size = 18) 
#   PLOT = PLOT + ggtitle(paste('Chimera FDR vs Cutoff -',INTERACTION,'Classification'))
#   PLOT = PLOT + geom_line(aes(colour='red'))
#   PLOT = PLOT + scale_x_continuous(limits = c(0,1)) + scale_y_continuous(limits = c(0,1))
#   #PLOT = PLOT + geom_vline(xintercept = 0.05, colour="red", linetype = "longdash")
#   PLOT = PLOT + theme(legend.position = "none")
#   PLOT
#   
#   # Classification over Cutoff
#   PLOT = ggplot(data = OUTPUT, aes(x = Cutoff))
#   PLOT = PLOT + geom_line(aes(y = tPos, colour = "tPos"))
#   #PLOT = PLOT + geom_line(aes(y = fNeg, colour = "fNeg"))
#   #PLOT = PLOT + geom_line(aes(y = tNeg, colour = "tNeg"))
#   PLOT = PLOT + geom_line(aes(y = fPos, colour = "fPos"))
#   #PLOT = PLOT + geom_line(aes(y = FDR, colour = "FDR"))
#  # PLOT = PLOT + scale_y_continuous(limits = c(0,250))
#   PLOT
#   
#   PLOT = PLOT + theme_grey(base_size = 18) 
#   PLOT = PLOT + ggtitle(paste('Chimera FDR vs Cutoff -',INTERACTION,'Classification'))
#   PLOT = PLOT + geom_line(aes(colour='red'))
#   PLOT = PLOT + scale_x_continuous(limits = c(0,1)) + scale_y_continuous(limits = c(0,1))
#   #PLOT = PLOT + geom_vline(xintercept = 0.05, colour="red", linetype = "longdash")
#   PLOT = PLOT + theme(legend.position = "none")
#   PLOT
#   
#   
# #     PLOT = ggplot(data = nn.df, aes(x = net.result))
# #       PLOT = PLOT + geom_density(aes(fill=factor(response)), alpha = 0.5)
# #       PLOT
#   
# ## Parameter Plots
#   # par(mfrow=c(2,2))
#   #   gwplot(nn,selected.covariate="Total", min=-10, max=10)
#   #   gwplot(nn,selected.covariate="exonRankInTranscript", min=-10, max=10)
#   #   gwplot(nn,selected.covariate="RepeatRank", min=-10, max=10)
#   #   gwplot(nn,selected.covariate="ThreadRatio",min=-10, max=10)
# 
# # end visualization

  
# CUTOFFS =====================================================================
#library(OptimalCutpoints)  

  # cut01 = optimal.cutpoints(data = nn.df, X = "net.result", status = "response", tag.healthy = 0,
  #                   methods = "ROC01") # Minimum Eucledian distance to (0,1) 
  # summary(cut01)
  #plot(cut01)
  
  # 
  # cutPPV = optimal.cutpoints(data = nn.df, X = "net.result", status = "response",  tag.healthy = 0,
  #                   methods = "MinValuePPV", valuePPV = 0.50) # Minimize PPV (1 - FDR)
  # summary(cutPPV)
  # plot(cutPPV)
  # 
  # cutMAX = optimal.cutpoints(data = nn.df, X = "net.result", status = "response",  tag.healthy = 0,
  #                            methods = "ValuePPV", valuePPV = 0.75)
  # summary(cutMAX)
  # plot(cutMAX)  
  # 
  # cutMAX = optimal.cutpoints(data = nn.df, X = "net.result", status = "response",  tag.healthy = 0,
  #                            methods = "MaxProdNPVPPV")
  # summary(cutMAX)
  # plot(cutMAX)   
  
  
  # Reducing Model Reps (Reduce NN model to single, best rep)
  
  # nn.reduced = neuralnet(
  #   data = DATA.ue,
  #   formula = TruePos ~ Total + exonRankInTranscript + RepeatRank + ThreadRatio + ExonMax + UpExonMax + RepeatMaxCoverage,
  #   hidden = 7,
  #   err.fct = "ce",
  #   linear.output = F,
  #   algorithm = 'rprop+',
  #   stepmax = 100,
  #   threshold = 10000,
  #   rep = 1,
  #   startweights = nn$weights[10])
  # 
