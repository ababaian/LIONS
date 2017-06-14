# rmUgtf_chromHMM.r
#
# aka: chimValid.r
# Framework to integreate epigenetic and orthogonal
# Data to the TE elements 
# 
#
# Functions
# # GTF processing
# 
# gtfParse = function(GTF) {
#   # Parse gtf Exon data to make all exons
#   # negative strand transcripts read as if they were
#   # on the positive strand
#   # (i.e. Exon 1 = first exon, start = TSS)
#   NegExons = which((GTF$strand == '-') & (GTF$class == 'exon'))
# }
#   

# Libraries
library(GenomicRanges)
library(ggplot2)
library(reshape2)
library(vioplot)

# REPEAT MASKER =====================================================

load('rmStats.Rdata')

# RMDB parse
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

# CHIMERIC OUTPUT ===================================================

# Load LION File
LION.output = read.csv(file = 'lions/k562.lion',header = T,sep = '\t')

#LION file Parse
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

# Deplete RINSIDE cases for now
#Chimera = Chimera[-which(Chimera$ER_Interaction == 'RInside')]

# CAGE ==============================================================
# Import CAGE Data from the UCSC 
# Downloaded November 26th 2013
# from UCSC table browser
CAGE.ucsc = read.table('CAGE_k562.ucsc',header=F,sep='\t')
colnames(CAGE.ucsc) = c('bin','chr','start','end','name',
                       'score','strand','FPKM','IDR','score2')

# Build CAGE Data into a GRanges Object
CAGE <- with(CAGE.ucsc,
             GRanges(chr, IRanges(start, end), strand, FPKM, IDR))

# mids are genomic points of the middle of each HMM cluster
mids = ranges(CAGE)

# Define mid-points for each CAGE HMM cluster
start(mids) = (start(ranges(CAGE))) + (width(ranges(CAGE))/2)
end(mids) = start(mids)
width(mids) = rep(1,times=(length(mids)))

# MidCAGE is middle defined CAGE
midCAGE = CAGE
ranges(midCAGE) = mids
strand(midCAGE) = Rle('*',lengths=length(midCAGE))

# sigCAGE are statistically significant CAGE clusters
# and smigCAGE are their respective midpoints
sigCAGE = CAGE[which(CAGE$IDR < 0.05)]
smigCAGE = midCAGE[which(midCAGE$IDR < 0.05)]
nsCAGE = midCAGE[-which(CAGE$IDR < 0.05)] # non-sig CAGE

# Cleanup
rm(mids, CAGE.ucsc)

# CAGE Intersect RM ==============================================
#
# ****** How to treat duplicates?? ******* #

rmUcage = as.matrix(findOverlaps(RM,CAGE))
rmUcagm = as.matrix(findOverlaps(RM,midCAGE))
rmUcasm = as.matrix(findOverlaps(RM,smigCAGE))

# Create Vector of CAGE U RepeatMasker
# and it's Recipricol

CAGE_arm = midCAGE[rmUcagm[,2]] # All rm CAGE clusters
CAGE_null = midCAGE[-rmUcage[,2]] # All non-rm CAGE clusters

CAGE_rm = sigCAGE[rmUcasm[,2]] # Significant CAGE clusters
CAGE_nsrm = CAGE_arm[-which(CAGE_arm$IDR < 0.05)] # Non-signifciant CAGE clusters

RM_cage = RM[rmUcasm[,1]] # Sig
RM_cagm = RM[rmUcagm[,1]] # ALL
RM_cagns = RM[rmUcagm[-which(CAGE_arm$IDR < 0.05),1]] # NS

# Removing Duplicate Repeats
uRM_cage = RM_cage[-which(duplicated(RM_cage)),]
uRM_cagm = RM_cagm[-which(duplicated(RM_cagm)),]
uRM_cagns = RM_cagns[-which(duplicated(RM_cagns)),]

# Repeat Family Vector
CAGE_arm$RM_cla = matrix(unlist(strsplit(RM_cagm$repName,split=' ')),byrow=T,ncol=3)[,2]
CAGE_null$RM_cla = c('null')

CAGE_rm$RM_cla = matrix(unlist(strsplit(RM_cage$repName,split=' ')),byrow=T,ncol=3)[,2]
CAGE_nsrm$RM_cla = matrix(unlist(strsplit(RM_cagns$repName,split=' ')),byrow=T,ncol=3)[,2]

# Total Repeat Family Vector (in bases)
RM_Class_bp = claStats$Bases[c(3,15,6,1,8)]
RM_Class_rel = (RM_Class_bp / sum(RM_Class_bp)) * 100

# Count Number of each RM Class in CAGE data
CAGE_Class = c(length( CAGE_rm$FPKM[which(CAGE_rm$RM_cla == 'LINE')]),
               length( CAGE_rm$FPKM[which(CAGE_rm$RM_cla == 'SINE')]),
               length( CAGE_rm$FPKM[which(CAGE_rm$RM_cla == 'LTR')]),
               length( CAGE_rm$FPKM[which(CAGE_rm$RM_cla == 'DNA')]),
               length( CAGE_rm$FPKM[which(CAGE_rm$RM_cla == 'Other')]) )
CAGE_Class_rel = (CAGE_Class / sum(CAGE_Class)) * 100

CAGE_Class_ns = c(length( CAGE_nsrm$FPKM[which(CAGE_nsrm$RM_cla == 'LINE')]),
                  length( CAGE_nsrm$FPKM[which(CAGE_nsrm$RM_cla == 'SINE')]),
                  length( CAGE_nsrm$FPKM[which(CAGE_nsrm$RM_cla == 'LTR')]),
                  length( CAGE_nsrm$FPKM[which(CAGE_nsrm$RM_cla == 'DNA')]),
                  length( CAGE_nsrm$FPKM[which(CAGE_nsrm$RM_cla == 'Other')]) )
CAGE_Class_ns_rel = (CAGE_Class_ns / sum(CAGE_Class_ns)) * 100

# CAGEuRM Intersection to LIONS ===============================================

cageUlion = as.matrix(findOverlaps(CAGE,LION))
cagmUlion = as.matrix(findOverlaps(midCAGE,LION))
casmUlion = as.matrix(findOverlaps(smigCAGE,LION))

# Parse CAGE intersection information into LION matrix
  cagePositive = numeric(length = length(LION$repName))
  cagePositive[unique(cagmUlion[,2])] = 1
  LION$cageStatus = cagePositive
  rm(cagePositive)

# Parse sigCAGE intersection information into LION matrix
  sigcagePositive = numeric(length = length(LION$repName))
  sigcagePositive[unique(casmUlion[,2])] = 1
  LION$sigCageStatus = sigcagePositive
  rm(sigcagePositive)

#Create Vector of CAGE U LIONS
# and it's Recipricol
LION_cagm = LION[cagmUlion[,2]] # All LIONS with CAGE clusters
CAGM_lion = midCAGE[cagmUlion[,1]] # All CAGE with LION clusters

RM_cage = RM[rmUcasm[,1]] # Sig
RM_cagm = RM[rmUcagm[,1]] # ALL
RM_cagns = RM[rmUcagm[-which(CAGE_arm$IDR < 0.05),1]] # NS

# cagmUlion = as.matrix(findOverlaps(uRM_cage, LION))
# casmUlion = as.matrix(findOverlaps(smigCAGE,LION))
# 
# # Parse CAGE intersection information into LION matrix
# cagePositive = numeric(length = length(LION$repName))
# cagePositive[unique(cagmUlion[,2])] = 1
# LION$cageStatus = cagePositive
# rm(cagePositive)
# 
# # Parse sigCAGE intersection information into LION matrix
# sigcagePositive = numeric(length = length(LION$repName))
# sigcagePositive[unique(casmUlion[,2])] = 1
# LION$sigCageStatus = sigcagePositive
# rm(sigcagePositive)

#Create Vector of CAGE U LIONS
# and it's Recipricol
LION_cagm = LION[cagmUlion[,2]] # All LIONS with CAGE clusters
CAGM_lion = midCAGE[cagmUlion[,1]] # All CAGE with LION clusters

RM_cage = RM[rmUcasm[,1]] # Sig
RM_cagm = RM[rmUcagm[,1]] # ALL
RM_cagns = RM[rmUcagm[-which(CAGE_arm$IDR < 0.05),1]] # NS

# Venn Diagram Matrix =========================================================

# CAGE_RM : all CAGE clusters within repeats
  # CAGE_RMns : Subset, non-signifincant clusters
  # CAGE_RMsig: Subset, significant clusters

# LION : all LION elements
  # LION_nc : Subset, non-CAGE LION elements
  # LION_csig : Subset, sigCAGE LION elements
  # LION_ns : Subset, non-significant CAGE LION elements

# Initialize Venn Diagram Matrix
  VENN = data.frame(row.names = c('CAGE_ns', 'CAGE_sig', 'LION', 'LION_Csig', 'LION_Cns'))
    VENN["LION","N"] = length( which(LION$cageStatus == 0))
    VENN["LION_Csig","N"] = length( which(LION$cageStatus == 1 & LION$sigCageStatus == 1))
    VENN["LION_Cns","N"] = length( which(LION$cageStatus == 1 & LION$sigCageStatus == 0))
    VENN["CAGE_sig","N"] = length( CAGE_rm ) - VENN["LION_Csig",'N']
    VENN["CAGE_ns","N"] = length( CAGE_arm ) - length( CAGE_rm )  - VENN["LION_Cns",'N']

# GRAPHING ====================================================================
# =============================================================================


# Set-up a data.frame to plot data
dfCAGE = data.frame(log(CAGE$FPKM),-log(CAGE$IDR))
colnames(dfCAGE) = c('FPKM','IDR')

dfCAGE_null = data.frame(log(CAGE_null$FPKM),-log(CAGE_null$IDR),CAGE_null$RM_cla)
colnames(dfCAGE_null) = c('FPKM','IDR','RM')

dfCAGE_arm = data.frame(log(CAGE_arm$FPKM),-log(CAGE_arm$IDR),CAGE_arm$RM_cla)
colnames(dfCAGE_arm) = c('FPKM','IDR','RM')

dfCAGE_all = rbind(dfCAGE_arm,dfCAGE_null)

# Plot FPKM vs IDR 
ggplot(dfCAGE, aes(FPKM, IDR)) + geom_point(alpha = 1/100)
ggplot(dfCAGE_arm, aes(FPKM, IDR) ) + xlim(c(-3,10)) + ylim(c(0,25)) + geom_point(aes(colour = factor(RM),alpha = 1/50))
ggplot(dfCAGE_null, aes(FPKM, IDR) ) + xlim(c(-3,10)) + ylim(c(0,25)) + geom_point(aes(colour = factor(RM),alpha = 1/50))


# # Plot Deviation of Repeat Classes from Expectence
# Class_Occurence = data.frame(
#   Class = factor(c("LINE","SINE","LTR","DNA","Other"),
#                  levels = c("LINE","SINE","LTR","DNA","Other")),
#   CAGE = CAGE_Class,
#   relCAGE = CAGE_Class_rel,
#   relCAGEns = CAGE_Class_ns_rel,
#   RM = RM_Class_bp,
#   relRM = RM_Class_rel,
#   delta = CAGE_Class_rel - RM_Class_rel,
#   delta_ns = CAGE_Class_ns_rel - RM_Class_rel)


# PLOT = ggplot( melt(Class_Occurence[,c('Class','delta','delta_ns')]),
#                aes(factor(Class), value, fill = variable))
# 
# PLOT = PLOT + geom_bar(stat = 'identity', position = 'dodge',
#                        colour = 'black')
# PLOT = PLOT + ggtitle("Relatiev Repeat Occurance in Significant & Non-significant CAGE clusters")
# PLOT = PLOT + xlab('') + ylab('Repeat Abundance relative to Genomic BP')
# PLOT

entry = c()
for (X in c('null','LINE','SINE','LTR','DNA','Other')){
  entry2 =
  c(length(which(dfCAGE_all$RM==X & dfCAGE_all$IDR < -log(0.05) )),
    length(which(dfCAGE_all$RM==X & dfCAGE_all$IDR > -log(0.05) )) )
  
  entry=rbind(entry,entry2)
}
entry
length(which(dfCAGE_all$RM=='LINE' & dfCAGE_all$IDR < 0.05))

# # Plot FPKM values of CAGE by Repeat Class (LINE, SINE, LTR, DNA, Other)
# vioplot( log(CAGE_arm$FPKM[which(CAGE_rm$RM_cla == 'LINE')],2),
#          log(CAGE_arm$FPKM[which(CAGE_rm$RM_cla == 'SINE')],2),
#          log(CAGE_arm$FPKM[which(CAGE_rm$RM_cla == 'LTR')],2),
#          log(CAGE_arm$FPKM[which(CAGE_rm$RM_cla == 'DNA')],2),
#          log(CAGE_arm$FPKM[which(CAGE_rm$RM_cla == 'Other')],2),
#          names=c("LINE", "SINE", "LTR", "DNA", "Other"))
#
#
# # Plot FPKM values of CAGE by Repeat Class (LINE, SINE, LTR, DNA, Other)
# vioplot( dfCAGE_all$FPKM[which(dfCAGE_all$RM == 'null' & dfCAGE_all$IDR <= 0.05)],
#          dfCAGE_all$FPKM[which(dfCAGE_all$RM == 'LINE' & dfCAGE_all$IDR <= 0.05)],
#          dfCAGE_all$FPKM[which(dfCAGE_all$RM == 'SINE' & dfCAGE_all$IDR <= 0.05)],
#          dfCAGE_all$FPKM[which(dfCAGE_all$RM == 'LTR' & dfCAGE_all$IDR <= 0.05)],
#          dfCAGE_all$FPKM[which(dfCAGE_all$RM == 'DNA' & dfCAGE_all$IDR <= 0.05)],
#          dfCAGE_all$FPKM[which(dfCAGE_all$RM == 'Other' & dfCAGE_all$IDR <= 0.05)],
#          names=c("Non-TE","LINE", "SINE", "LTR", "DNA", "Other"))
#
# # Plot FPKM values of CAGE by Repeat Class (LINE, SINE, LTR, DNA, Other)
# vioplot( dfCAGE_all$FPKM[which(dfCAGE_all$RM == 'null')],
#          dfCAGE_all$FPKM[which(dfCAGE_all$RM == 'LINE')],
#          dfCAGE_all$FPKM[which(dfCAGE_all$RM == 'SINE')],
#          dfCAGE_all$FPKM[which(dfCAGE_all$RM == 'LTR')],
#          dfCAGE_all$FPKM[which(dfCAGE_all$RM == 'DNA')],
#          dfCAGE_all$FPKM[which(dfCAGE_all$RM == 'Other')],
#          names=c("Non-TE","LINE", "SINE", "LTR", "DNA", "Other"))
#
## Plot IDR values as their negative log
# vioplot( -log(CAGE_null$IDR[which(CAGE_null$IDR < 0.05)],10),
#          -log(CAGE_rm$IDR[which(CAGE_rm$RM_cla == 'LINE')],10),
#          -log(CAGE_rm$IDR[which(CAGE_rm$RM_cla == 'SINE')],10),
#          -log(CAGE_rm$IDR[which(CAGE_rm$RM_cla == 'LTR')],10),
#          -log(CAGE_rm$IDR[which(CAGE_rm$RM_cla == 'DNA')],10),
#          -log(CAGE_rm$IDR[which(CAGE_rm$RM_cla == 'Other')],10) )
#  
# ###
# ggplot(dfCAGE_all, aes(IDR) ) +
#   stat_bin(data = subset(dfCAGE_all,RM == 'null'), fill = 'white', binwidth = 1/10) +
#   stat_bin(data = subset(dfCAGE_all,RM == 'LINE'), fill = 'red', alpha = 0.1 , binwidth = 1/10) +
#   stat_bin(data = subset(dfCAGE_all,RM == 'SINE'), fill = 'green', alpha = 0.2 , binwidth = 1/10) + 
#   stat_bin(data = subset(dfCAGE_all,RM == 'LTR'), fill = 'orange', alpha = 0.3 , binwidth = 1/10) +
#   stat_bin(data = subset(dfCAGE_all,RM == 'DNA'), fill = 'gray', alpha = 0.5 , binwidth = 1/10) +scale_y_log10()

###
sigCAGE_all=dfCAGE_all[which(dfCAGE_all$IDR > -log(0.05)),]

ggplot(sigCAGE_all, aes(FPKM) ) +
  stat_bin(data = subset(sigCAGE_all,RM == 'null'),fill = ' dark blue', binwidth = 1/5 , right = T) +
  stat_bin(data = subset(sigCAGE_all,RM != 'null'),fill = 'orange', binwidth = 1/5 , right = T) + scale_y_log10() +
#   stat_bin(data = subset(sigCAGE_all,RM == 'LINE'), geom = 'line', colour = 'red', binwidth = 0.3) +
#   stat_bin(data = subset(sigCAGE_all,RM == 'SINE'), geom = 'line', colour = 'green', binwidth = 0.3) + 
#   stat_bin(data = subset(sigCAGE_all,RM == 'LTR'), geom = 'line', colour = 'orange', binwidth = 0.3) +
#   stat_bin(data = subset(sigCAGE_all,RM == 'DNA'), geom = 'line', colour = 'purple', binwidth = 0.3) + 
#   stat_bin(data = subset(sigCAGE_all,RM == 'Other'), geom = 'line', colour = 'gray', binwidth = 0.3) + scale_y_log10() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# Histogram of IDRs
ggplot(dfCAGE_all, aes(IDR) ) +
  stat_bin(data = subset(dfCAGE_all,RM == 'null'), binwidth = 1/5, 
           drop = FALSE, right = TRUE, alpha = 1, fill = "dark blue")  +
  stat_bin(data = subset(dfCAGE_all,RM != 'null'), binwidth = 1/5, 
           drop = FALSE, right = TRUE, alpha = 1, fill = "orange")  +  
  xlim(c(0,24)) +  scale_y_log10() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# Histogram of Significant FPKM
sigCAGE_all=dfCAGE_all[which(dfCAGE_all$IDR > -log(0.05)),]
ggplot(sigCAGE_all, aes(FPKM) ) +
  stat_bin(data = subset(sigCAGE_all,RM == 'null'),fill = ' dark blue', binwidth = 1/5 , right = T) +
  stat_bin(data = subset(sigCAGE_all,RM != 'null'),fill = 'orange', binwidth = 1/5 , right = T) + scale_y_log10() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# # Reproducible CAGE Expression by TE class
# ggplot(sigCAGE_all, aes(RM, FPKM) ) +
#   geom_violin() +
#   geom_boxplot(width=.1) +
#   #scale_y_log10() +
#   theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
#                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# 
# 
# ggplot(sigCAGE_all, aes(RM) ) +
#   geom_bar() + 
#   scale_y_log10(limits = c(1, 1000000)) +
#   theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
#                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# 
# # All CAGE expression by TE Class
# ggplot(dfCAGE_all, aes(RM, FPKM) ) +
#   geom_violin() +
#   geom_boxplot(width=.1) +
#   scale_y_log10() +
#   theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
#                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggplot(sigCAGE_all[which(sigCAGE_all$RM %in% c('DNA','LINE','LTR','Other','SINE','null')),], aes(FPKM, fill = RM)) +
  geom_histogram(position="fill", binwidth = 1/5)

ggplot() + 
  geom_bar(data = sigCAGE_all[which(sigCAGE_all$RM %in% c('DNA','LINE','LTR','Other','SINE','null')),], aes(FPKM, colour = RM),
           aes(x = factor(RM),fill = factor(FPKM)),
           position = "fill")

# Cumulative Ditribution Function
ggplot(sigCAGE_all[which(sigCAGE_all$RM %in% c('DNA','LINE','LTR','Other','SINE','null')),], aes(FPKM, colour = RM) ) +
  stat_ecdf() +
  labs(title="Empirical Cumulative \n Density Function",
       y = "F(Repeats)", x="log(FPKM)")+
  scale_colour_manual(values=c('purple','red','orange','grey','green','blue')) +
  theme_classic()

# Students T-test
t.test(dfCAGE_all$IDR[which(dfCAGE_all$RM=='null')],dfCAGE_all$IDR[which(dfCAGE_all$RM!='null')],var.equal = F)
t.test(sigCAGE_all$FPKM[which(sigCAGE_all$RM=='null')],sigCAGE_all$FPKM[which(sigCAGE_all$RM!='null')],var.equal = F)
  