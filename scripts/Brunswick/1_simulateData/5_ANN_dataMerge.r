# ANN_dataMerge.r
#
# Merge the different data sets
# and randomize them for training / testing
#

# DATA SETS =================================================================



# TRAINING DATASET ==========================================================
# Append this dataset to entire data training set

rm(DATA_hold)

load('h1_10.Rdata'); DATA_hold = DATA
  load('h1_50.Rdata'); DATA_hold = rbind(DATA,DATA_hold)
  load('h1_100.Rdata'); DATA_hold = rbind(DATA,DATA_hold)
  load('h1_200.Rdata'); DATA_hold = rbind(DATA,DATA_hold)
  load('k562_25.Rdata'); DATA_hold = rbind(DATA,DATA_hold)
  load('k562_100.Rdata'); DATA_hold = rbind(DATA,DATA_hold)
  load('k562_200.Rdata'); DATA_hold = rbind(DATA,DATA_hold)
DATA = DATA_hold; rm(DATA_hold)


# ===================================================================
# ===================================================================
# ===================================================================

# Save Complete Data set prior to subsetting
 save(DATA,file= 'DATA_all.Rdata')

# Summary Statistics 
# Insert a file output of statistics
# about the input

# Data Randomization
# split the data into thirds:
# 2/3 for training
# 1/3 for testing the resultant ANN
  DATA$rand = sample(1:3, length(DATA$TruePos), replace=T) 
  DATA_hold = DATA  

# Save 2/3 of the data for training 
  DATA = DATA[which(DATA$rand <= 2),]
  DATA = DATA[,-35]
  save(DATA,file='DATA_train.Rdata')

# Save 1/3 of the data for testing
  DATA = DATA_hold
  DATA = DATA[which(DATA$rand == 3),]
  DATA = DATA[,-35]
  save(DATA,file='DATA_test.Rdata')
