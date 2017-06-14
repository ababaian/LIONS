# ANN_graph.R
#
# Graph the 'LIONS' parameters to distinguish
# 
#

  library('ggplot2')

# # Load DATA (LCSV) for visualization
#   load('DATA_train.Rdata')
#     DATA1 = DATA
#   load('DATA_test.Rdata')
#     DATA = rbind(DATA1,DATA)
#     rm(DATA1)

# Load DATA (LCSV) for visualization (all)
  load("DATA_all.Rdata")

# Breakdown DATA by classification
  
  upDATA = DATA[which(DATA$ER_Interaction == "Up"),]
  ueDATA = DATA[which(DATA$ER_Interaction == "UpEdge"),]
  dwDATA = DATA[which(DATA$ER_Interaction == "Down"),]
  deDATA = DATA[which(DATA$ER_Interaction == "DownEdge"),]
  eiDATA = DATA[which(DATA$ER_Interaction == "EInside"),]
  riDATA = DATA[which(DATA$ER_Interaction == "RInside"),]

print('ANN DATA Statistics -------------')
print('')
print(' Number of True Positives in each set')
  print( paste(" Up: ", sum(upDATA$TruePos)) )
  print( paste(" UpEdge: ", sum(ueDATA$TruePos)) )
  print( paste(" Down: ", sum(dwDATA$TruePos)) )
  print( paste(" DownEdge: ", sum(deDATA$TruePos)) )
  print( paste(" EInside: ", sum(eiDATA$TruePos)) )
  print( paste(" RInside: ", sum(riDATA$TruePos)) )
 
# Review how TruePos are being selected; this may be
# non-optimal and thus the entire thing is a bit screwy

# DATA NORMALIZATION =========================================
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

# ================================================================
  
DATA.ue = upDATA
  
# Thread Ratio
  DATA.ue$ThreadRatio = DATA.ue$DownThread/DATA.ue$UpThread
    zeroHold = is.nan(DATA.ue$ThreadRatio) # 0/0 set to min
    oneHold = is.infinite(DATA.ue$ThreadRatio) # x/0 set to max
  
  DATA.ue$ThreadRatio = log10(DATA.ue$ThreadRatio) 
  DATA.ue$ThreadRatio = normalize(DATA.ue$ThreadRatio)
    DATA.ue$ThreadRatio[zeroHold] = 0
    DATA.ue$ThreadRatio[oneHold] = 1
  
# Exon Rank in Transcript
  DATA.ue$exonRankInTranscript[ DATA.ue$exonRankInTranscript == 1] = 0
  DATA.ue$exonRankInTranscript[ DATA.ue$exonRankInTranscript > 1] = 1
  
# Total
  DATA.ue$Total = log10(DATA.ue$Total)
  DATA.ue$Total = normalize(DATA.ue$Total)

# ================================================================
  
# Thread Ratio Density Plot 
  PLOT = ggplot(data = DATA.ue, aes(x = Total, fill = factor(DATA.ue$TruePos), alpha = 0.1))
  PLOT = PLOT + geom_density()
  PLOT

# Thread Ratio Density Plot 
  PLOT = ggplot(data = DATA.ue, aes(y = ThreadRatio, x = Total, color = factor(DATA.ue$TruePos), alpha = 0.1))
  PLOT = PLOT + geom_jitter()
  PLOT

# Thread Ratio Density Plot 
  PLOT = ggplot(data = DATA.ue, aes(y = ThreadRatio, z = DATA.ue$exonRankInTranscript, x = Total))
  PLOT = PLOT + geom_tile()
  PLOT
  
#   PLOT = ggplot(data = DATA.ue, aes(x = TruePos, y = ThreadRatio, alpha = 0.1))
#   PLOT = PLOT + geom_jitter()
#   PLOT
