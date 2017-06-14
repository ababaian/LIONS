# rmUgtf_assembly.r
#
# Framework to integreate epigenetic and orthogonal
# Data to the TE elements 
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

# ASSEMBLY ==========================================================

# Load Assembly GTF File
# from Cufflinks

ASSEMBLY = read.csv(file = 'assembly/k562.assembly.gtf', header = F, sep ='\t')
  colnames(ASSEMBLY) = c('chr','source','feature','start','end',
                         'score','strand','frame','attribute') # GTF headers

# Transcripts 
TRANSCRIPTS = ASSEMBLY[which(ASSEMBLY$feature == 'transcript'),] # exclude exon descriptions

# Extract All Exons
EXONS = ASSEMBLY[which(ASSEMBLY$feature == 'exon'),] # All exons

# Exon Attributes
ATT = as.data.frame(do.call(rbind, strsplit(as.character(EXONS$attribute),split = ';')))
  colnames(ATT) = c('gene_id','transcript_id','exon_number',
                    'FPKM','frac','conf_lo','conf_hi','cov')
  ATT = as.matrix(ATT)

N = 1
for (COL in colnames(ATT)) {
  ATT[,N] =
        as.vector(
        gsub( x = as.vector(ATT[,N]),
        pattern = paste('',COL,''),
        replacement = '' ) )
  N = N + 1
}
ATT[,1] = as.vector( gsub('gene_id ', '', ATT[,1]))

EXONS = as.data.frame(cbind(as.matrix(EXONS[,-9]), ATT))
  levels(EXONS$strand) = c('-', '*', '+') # Set to GRanges format

# # Assembly of Exon 1 GRangse Objects
# COORDS = cbind(EXONS1$start, EXONS1$end)
# ASSex1 <- GRanges(seqnames = EXONS1$chr,
#                ranges = IRanges( start = apply(COORDS, 1, min),
#                                  end = apply(COORDS, 1, max),
#                                  names = EXONS1$strand
#                                  ) ,
#               strand = EXONS1$strand,
#               FPKM = as.numeric(EXONS1$FPKM) )

# Extract and parse All Exon 1s in assembly

# Extract All First Exons --------------------------------
  EXONS1 = EXONS[which(EXONS$exon_number == 1),]

# TSS parse to GR
  COORDS = cbind(EXON1$start, EXONS$end)

  TSS <- GRanges(seqnames = EXONS1$chr,
                 ranges = IRanges( start = COORDS[,1], # Start = TSS
                                   end = COORDS[,1], # End = TSS
                                   names = EXONS1$strand
                 ) ,
                 strand = EXONS1$strand,
                 FPKM = as.numeric(EXONS1$FPKM) )

  # only stranded TSS
    strTSS = TSS[-which(as.vector(strand(TSS1)) == "*"),]
    strTSS$EStrand = as.vector(strand(strTSS)) # Parse into a vector
    # strand(strTSS) = c("*")

# Extract only spliced transcript exons
  EXONsp = EXONS[which( duplicated(EXONS$transcript_id) |
                        duplicated(EXONS$transcript_id,fromLast = T)), ]

# Extract Spliced First Exons
  EXONsp1 = EXONsp[which(EXONsp$exon_number == 1),]

#Assembly file Parse
COORDS = cbind(EXONS$start, EXONS$end)

# TSS file Parse
  COORDS = cbind(EXON1$start, EXONS$end)
  TSS <- GRanges(seqnames = EXONS1$chr,
                    ranges = IRanges( start = COORDS[,1], # Start = TSS
                                      end = COORDS[,1], # End = TSS
                                      names = EXONS1$strand ) ,
                    strand = EXONS1$strand,
                    FPKM = as.numeric(EXONS1$FPKM) )
  
  # only stranded TSS
  strTSS = TSS[-which(as.vector(strand(TSS1)) == "*"),]
    strTSS$EStrand = as.vector(strand(strTSS)) # Parse into a vector
    # strand(strTSS) = c("*")

# Spliced TSS file Parse
  spTSS <- GRanges(seqnames = EXONsp1$chr,
                 ranges = IRanges( start = COORDS[,1], # Start = TSS
                                   end = COORDS[,1], # End = TSS
                                   names = EXONS1$strand
                 ) ,
                 strand = EXONS1$strand,
                 FPKM = as.numeric(EXONS1$FPKM) )
  
  # only stranded TSS
  strTSS = TSS[-which(as.vector(strand(TSS1)) == "*"),]
  strTSS$EStrand = as.vector(strand(strTSS)) # Parse into a vector
  # strand(strTSS) = c("*")

# Clean-up
rm(ATT,COORDS,EXONS1,TRANSCRIPTS,ASSEMBLY)

# CHIMERIC OUTPUT ===================================================

# Load LION File
LION.output = read.csv(file = 'lions/k562.lion', header = T, sep = '\t')

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

# Assembly Intersect RM ==============================================
#
# ****** How to treat duplicates?? ******* #

  rmUtss = as.matrix(findOverlaps(RM,TSS))
  rmUstss = as.matrix(findOverlaps(RM,strTSS))

# Create Vector of Assembly strTSS U RepeatMasker
# and it's Recipricol
# Note: Multiple TSSs can land in the same repeat

  strTSS_rm = strTSS[rmUstss[,2]]

  RM_strTSS = RM[rmUstss[,1]]
    RM_strTSS = RM_strTSS[-which(duplicated(RM_strTSS)),] # Remove Duplicate Repeats


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

