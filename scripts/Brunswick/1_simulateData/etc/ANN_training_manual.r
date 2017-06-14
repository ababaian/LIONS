# ANN_training.R
#
# Framework to integreate epigenetic and orthogonal
# Data to the TE elements 
# 

# Libraries
library(GenomicRanges)
library(ggplot2)
library(reshape2)
#library(vioplot)

# Control Panel =====================================================

# Input 'Expressed' Transcriptome
  # Based on FPKM from reference transcriptome
  # output from fluxSimExpToGTF.r script
  expGTF = 'h1_200m.gtf'

# Input LIONS file and pc.lcsv file
  # Output from LIONS pipeline
  inLIONS = 'h1_200.lions'
  inLCSV = 'h1_200.pc.lcsv'

# Output training DATA File (.Rdata)
  DATA_OUTPUT = 'tmp.Rdata'

# Pro file from FluxSim
  FluxPRO = 'h1_200m.pro'

# REPEAT MASKER =====================================================

# Load Repeat Masker Data set
  load('rmStats.Rdata')

# RMDB parse for entry into GRanges
  # If complementary strand; swap out repStart and repLeft
  # to make everythign 'consensus' oriented
  DepleteRMDB = which(RMDB$repClass %in% c('LINE', 'SINE', 'DNA', 'LTR',
                                           'DNA?', 'LTR?', 'Other', 'Unknown'))

  RMDB = RMDB[DepleteRMDB,]

# Hold repStarts in a vector
  NegStrand = which(RMDB$strand == '-')
  RS_hold = RMDB[NegStrand, 'repStart']

# Swap columns
  RMDB[NegStrand, 'repStart'] = RMDB[NegStrand, 'repLeft']
  RMDB[NegStrand, 'repLeft'] = RS_hold

# CleanUp
  rm(RS_hold, NegStrand,DepleteRMDB)

# Build GRanges Object from RMDB
  RM <- with(RMDB,
             GRanges(genoName, IRanges(genoStart, genoEnd), c('*'), # Coords
                     repName = paste(RMDB$repName, RMDB$repClass, RMDB$repFamily),
                     swScore, milliDiv, milliDel, milliIns, RStrand = strand))
  rm(RMDB)

# Transcriptome =====================================================

# Load Transcriptome GTF File (For the simulation)
  # H1 Gencode Subset from simulation (Total Transcripts)
  # Output from subsetTranscripts.sh
  #TRANSCRIPTOME = read.csv(file = 'h1Transcriptome.gtf', header = F, sep = '\t')

  # H1 Gencode Expressed subset from simulation
  TRANSCRIPTOME = read.csv(file = expGTF, header = F, sep = '\t')


  colnames(TRANSCRIPTOME) = c('chr','source','feature','start','end',
                              'score','strand','frame','attribute') # GTF headers
  
  
  TRANSCRIPTS = TRANSCRIPTOME[which(TRANSCRIPTOME$feature == 'transcript'),] # exclude exon descriptions

# Extract All Exons
  EXONS = TRANSCRIPTOME[which(TRANSCRIPTOME$feature == 'exon'),] # All exons

# Exon Attributes
  ATTR = as.data.frame(do.call(rbind, strsplit(as.character(EXONS$attribute),split = ';')))
  # This may throw an error of columns not being the same. Just make sure that all Col 1:10 are
  # correct and you can ignore it.

# Go through the attribute matrix and reduce it to the data
  ATTR = as.matrix(ATTR)[,1:10]

    # From Gencode; use these column names. There are more but discard them
    # (10 columns)
    colnames(ATTR) = c('gene_id', 'transcript_id', 'gene_type', 'gene_status', 'gene_name',
                       'transcript_type', 'transcript_status', 'transcript_name',
                       'exon_number', 'exon_id')

  N = 1
  for (COL in colnames(ATTR)) {
    ATTR[,N] =
      as.vector(
        gsub( x = as.vector(ATTR[,N]),
              pattern = paste('',COL,''),
              replacement = '' ) )
    N = N + 1
  }

  # Rename first column
  ATTR[,1] = as.vector( gsub('gene_id ', '', ATTR[,1]))

# Recompile EXONS matrix
  EXONS = as.data.frame(cbind(as.matrix(EXONS[,-9]), ATTR))
  #levels(EXONS$strand) = c('-', '*', '+') # Set to GRanges format


# FIRST EXONS ---------------------------
# Extract First Exons
  EXONS1 = EXONS[which(EXONS$exon_number == 1),]

# TSS parse to GR
  COORDS = cbind(as.numeric(as.vector(EXONS1$start)),
                 as.numeric(as.vector(EXONS1$end)),
                 as.character(EXONS1$strand) )

# If Exon1 is on the negative strand, switch Start and End Coordinates
# to make left-most position the 'Start Site'
  COORDS = matrix(
    apply(COORDS,1,function(x){
    if (x[3] == "-"){
      return( x[ c(2,1,3) ] )
    } else {
      return( x[1:3] )
    }
    }), ncol = 3, byrow = T)

  TSS <- GRanges(seqnames = EXONS1$chr,
                 ranges = IRanges( start = as.numeric(COORDS[,1]), # Start = TSS
                                   end = as.numeric(COORDS[,1]), # End = TSS
                                   names = EXONS1$strand
                 ),
                 strand = EXONS1$strand )

#Cleanup
#  rm(EXONS1,COORDS,ATTR,N,COL)

# Spliced Exons ------------------------

# Extract only spliced transcript exons
  EXONsp = EXONS[which( duplicated(EXONS$transcript_id) |
                        duplicated(EXONS$transcript_id,fromLast = T)), ]

# Extract Spliced First Exons
  EXONsp1 = EXONsp[which(EXONsp$exon_number == 1),]


# TSS parse to GR
  COORDS = cbind(as.numeric(as.vector(EXONsp1$start)),
                 as.numeric(as.vector(EXONsp1$end)),
                 as.character(EXONsp1$strand) )

# If Exon1 is on the negative strand, switch Start and End Coordinates
# to make left-most position the 'Start Site'
  COORDS = matrix(
    apply(COORDS,1,function(x){
      if (x[3] == "-"){
        return( x[ c(2,1,3) ] )
      } else {
        return( x[1:3] )
      }
    }), ncol = 3, byrow = T)

  spTSS <- GRanges(seqnames = EXONsp1$chr,
                   ranges = IRanges( start = as.numeric(COORDS[,1]), # Start = TSS
                                     end = as.numeric(COORDS[,1]), # End = TSS
                                     names = EXONsp1$strand
                   ),
                   strand = EXONsp1$strand )

# LIONS OUTPUT ===================================================

# Load .LION File
  LION.output = read.csv(file = inLIONS, header = T, sep = '\t')

# .LION file Parse
  LION <- with(LION.output,
               GRanges(paste(c('chr'),LION.output$Chromosome,sep=''),
                       IRanges(RStart, REnd), c('*'), # Coords
                       repName = LION.output$repeatName,
                       ER_Interaction,
                       RepeatRPKM,
                       RefID,
                       assXref,
                       Contribution,
                       RStrand))

# LION rmid
  LION.rmid = LION.output$RepeatID


# Load the *.pc.lcsv [Complete RM description, pre-classification] 
  LCSV.output = read.csv(file = inLCSV, header = T, sep = '\t')
  LCSV.rmid = paste('chr', LCSV.output$Chromosome, ':', LCSV.output$RStart, sep = '')

# Assembly Intersect RM ==============================================
#
# ****** How to treat duplicates?? ******* #

# Find Intersections of RepeatMasker and Genocde TSS (and spliced TSS)
  rmUtss = as.matrix(findOverlaps(RM,TSS))
  rmUstss = as.matrix(findOverlaps(RM,spTSS))

# Create Vector of Assembly spTSS U RepeatMasker
  # and it's Recipricol
  # Note: Multiple TSSs can land in the same repeat
  spTSS_rm = spTSS[rmUstss[,2]]

  RM_spTSS = RM[rmUstss[,1]]
  RM_spTSS = RM_spTSS[-which(duplicated(RM_spTSS)),] # Remove Duplicate Repeats

# Make Repeat_ID tags for the intersection set
  RM.sptss.id = as.vector(seqnames(RM_spTSS))
    RM.sptss.id = paste(RM.sptss.id, start(RM_spTSS),sep=':')

# Clean-up
#  rm(EXONsp,EXONsp1,COORDS)

# FluxPro ===========================================================

# H1 Gencode Expressed subset from simulation
  PRO = read.csv(file = FluxPRO, header = F, sep = '\t')
  
    colnames(PRO) = c('Locus', 'TranID', 'Coding', 'TransLength', 'ExpressedFraction',
                      'ExpressedNumber','LibraryFraction', 'LibraryNumber', 'SeqFraction',
                      'SeqNumber','Covered', 'ChiSquare','CoefVar')
    LIBSIZE = sum(PRO$SeqNumber)
    
  # Calculate RPKM
  PRO$RPKM = (1e9 * (PRO$SeqNumber) / (PRO$TransLength) * (1/LIBSIZE))


# Which Transcripts are TE-initaited; which are detected by LCSV and LIONS
  # TE Initiated Transcripts
  PRO$TE = '-'
  PRO$TE[which(as.vector(PRO$TranID) %in% as.vector(EXONsp1$transcript_id[rmUstss[,2]]) )] = "+"

  # LCSV Transcripts
  SIM.sptss.id = as.vector(seqnames(RM_spTSS[SIM.tp,]))
    SIM.sptss.id = paste(SIM.sptss.id, start(RM_spTSS[SIM.tp]),sep=':')

    # I need the conversion between SIM.sptss.id and RM.sptss.id
      index1 = c(1:length(RM.sptss.id))
      index2 = which(duplicated(RM[rmUstss[,1]])) + 0.5
      mergdex = order(c(index1,index2))
    longRM.sptss.id = c(RM.sptss.id,length(index2))[mergdex]

  PRO$LCSV = "-"
  PRO$LCSV[which(as.vector(PRO$TranID) %in% as.vector(EXONsp1$transcript_id[rmUstss[,2]][which(longRM.sptss.id %in% SIM.sptss.id)]))] = "+"


# 1 Locus  chrom:start-end[W|C]  identifier of the transcriptional locus, given by the chromosome (chrom), start respectively end position, and the strand (Watson or Crick).
# 2  Transcript_ID	String	transcript identifier from the reference annotation.
# 3	Coding	[CDS|NC]	specifies whether the transcript has an annotated coding sequence (CDS) or not (NC)
# 4	Length	Integer	the mature length of the transcript after splicing out introns, disregarding the poly-A tail, as annotated in the reference annotation
# 5	Expressed Fraction	Float	fraction of RNA molecules that represent transcripts that are qualitatively equal to this RNA form
# 6	Expressed Number	Integer	absolute number of expressed RNA molecules
# 7	Library Fraction	Float	fraction of cDNA molecules in the final library that have been produced from this transcript
# Library Number	Integer	absolute number of cDNA fragments generated from this transcript
# 9	Sequenced Fraction	Float	fraction of total reads that have been sequenced from this transcript
# 10	Sequenced Number	Integer	absolute number of reads sequenced from this transcript
# 11	Covered Fraction	Float	fraction of the transcript that is covered by reads
# 12	Chi Square	Integer	chi-square goodness of fit measurement of coverage uniformity
# 13	Coefficient of Variation	Float	coefficient of variation for transcript coverage


#png(filename = paste("ROC_",INTERACTION,".png"), width = 800, height = 800)


plotTE = factor(PRO$TE, levels = c("+","-"))
plotLCSV = factor(PRO$LCSV, levels = c("+","-"))
  LCSV_pos = plotLCSV[which(PRO$TE == "+")]

  PLOT = ggplot(data = PRO, aes(x = RPKM, fill = plotTE) )
  PLOT = PLOT + geom_histogram(bins = 100)  + scale_x_log10()
  PLOT

  PLOT = ggplot(data = PRO[which(PRO$TE == "+"),], aes(x = RPKM, fill = LCSV_pos))
  PLOT = PLOT + geom_histogram()  + scale_x_log10()
  PLOT = PLOT + scale_fill_manual(values = c("dark green",'grey'))
  PLOT



# Neural Network Training  Data Set ==================================

# Beta Testing
# # All RM ID as test
#   allRMID = as.vector(seqnames(RM))
#   allRMID = paste(allRMID, start(RM),sep=':')

  # Which Simulation GTF are in .LCSV [ True Positives ]
  SIM.tp = which(RM.sptss.id %in% LCSV.rmid)
  LION.tp = which(RM.sptss.id %in% LION.output$RepeatID)
  LION.tp2 = which(LION.output$RepeatID %in% RM.sptss.id)
  LION_SIM.tp = which(LION.output$RepeatID %in% LCSV.rmid)

  # Which Simulation GTF are missing from LCSV [ False Negative]
  SIM.fn = which(!(RM.sptss.id %in% LCSV.rmid))
  LION.fn = which(!(RM.sptss.id %in% LION.output$RepeatID))

  # Note: Compare "Expression Level" between these TP and FN
  # transcripts. Intuition is that poor mapping and low expression
  # likely are associated with overall RNAseq FN

  # Which LCSV entries are Simulation GTF positive [True Positive Training Set]
  LCSV.tp = which(LCSV.rmid %in% RM.sptss.id)


  # Annotate the .LCSV file with 'True Positives's
  LCSV.output$TruePos = 0
  LCSV.output$TruePos[LCSV.tp] = 1

  LCSV.tp.output = LCSV.output[which(LCSV.output$TruePos == 1),]

  # Note: At this point there's a little pickle
  # I have a list of true positive cases
  # yet not every 'interaction' on here may neccesarily be a "True Positive"
  # interaction. ???

# Summary Statistics ========================================================

# Number of TSS in Expressed GTF
  N_tss = length(TSS) # All
  N_sp_tss = length(spTSS) # Spliced

# Number of RM TSS 
  N_rm_tss = length(rmUtss[,1]) # All ## Add UNIQUE and recalcualte
  N_rm_stss = length(rmUstss[,1]) # Spliced ## Add UNIQUE and recalcualte

# Number of RM TSS which show up in LCSV
  N_lcsv_tp = length(SIM.tp)
  N_lcsv_fn = length(SIM.fn) 

# Number of RM TSS detected by LIONS output
  N_lion = length(LION)
  N_lion_tp = length(LION.tp)
  N_lion_fn = length(LION.fn)

c(N_tss, N_sp_tss, N_rm_tss, N_rm_stss, N_lcsv_tp, N_lcsv_fn, N_lion, N_lion_tp, N_lion_fn)


c(length(which(LION.output$ER_Interaction == "Up")),
  length(which(LION.output$ER_Interaction == "UpEdge")),
  length(which(LION.output$ER_Interaction == "EInside")),
  length(which(LION.output$ER_Interaction == "RInside")))

c(length(which(LCSV.output$ER_Interaction == "Up")),
  length(which(LCSV.output$ER_Interaction == "UpEdge")),
  length(which(LCSV.output$ER_Interaction == "EInside")),
  length(which(LCSV.output$ER_Interaction == "RInside")))

c(length(which(LCSV.output$ER_Interaction == "Up" & LCSV.output$TruePos ) ),
  length(which(LCSV.output$ER_Interaction == "UpEdge" & LCSV.output$TruePos ) ),
  length(which(LCSV.output$ER_Interaction == "EInside" & LCSV.output$TruePos ) ),
  length(which(LCSV.output$ER_Interaction == "RInside" & LCSV.output$TruePos ) ))

c( length( which(LION.output$RepeatID[which(LION.output$ER_Interaction == "Up")] %in% RM.sptss.id)),
   length( which(LION.output$RepeatID[which(LION.output$ER_Interaction == "UpEdge")] %in% RM.sptss.id)),
   length( which(LION.output$RepeatID[which(LION.output$ER_Interaction == "EInside")] %in% RM.sptss.id)),
   length( which(LION.output$RepeatID[which(LION.output$ER_Interaction == "RInside")] %in% RM.sptss.id))
)

# TRAINING DATASET ==========================================================
# Append this dataset to entire data training set
#

#  if (exists(x='DATA')){
#    # Append LCSV.output to all the dataset
#    common.names <- intersect(colnames(DATA), colnames(LCSV.output))
#    DATA <- rbind(DATA[, common.names], LCSV.output[, common.names])
#  } else {
#    # Assign CHIMERA as the database
#    DATA = LCSV.output
#  }    

# ===================================================================
# ===================================================================
# ===================================================================

# Save Complete Data set prior to subsetting
#  save(DATA,file= DATA_OUTPUT)


# Summary Statistics 
  # Insert a file output of statistics
  # about the input

# # Data Randomization
# # split the data into thirds:
# # 2/3 for training
# # 1/3 for testing the resultant ANN
#   DATA$rand = sample(1:3, length(DATA$TruePos), replace=T) 
#   DATA_hold = DATA  
# 
# # Save 2/3 of the data for training 
#   DATA = DATA[which(DATA$rand <= 2),]
#   DATA = DATA[,-35]
#   save(DATA,file='DATA_train.Rdata')
# 
# # Save 1/3 of the data for testing
#   DATA = DATA_hold
#   DATA = DATA[which(DATA$rand == 3),]
#   DATA = DATA[,-35]
#   save(DATA,file='DATA_test.Rdata')
# 
