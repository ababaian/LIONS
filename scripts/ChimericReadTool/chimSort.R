# chimSort.R
# March 23 2015 - Modified for LIONS
# ---------------------------------------------------------
# Usage:
#   chimSort.R <ChimericResults> <FilteredOutputName> <Number Exonic Reads>
# 
# Read ChimericReadTool.sh output
#   Parse data into Rdata
#   Calculate derived values
#   Sort/Classify Chimera Cases
#   Output Rdata of sorted chimeric cases
# ---------------------------------------------------------

# CONTROL PANEL ===============================================================

  STDIN = commandArgs(trailingOnly = TRUE)

  INPUT = STDIN[1] # Chimera_Results

  OUTPUT = STDIN[2] # Chimeric Output

  EXONIC_READS = as.numeric(STDIN[3]) # Number of exonic reads in library

# Selection Criteria
# 3 classes of Selection Criteria are defined
  # Upstream TE Interaction
    # Interaction = Upstream
    # Chimeric Reads >= 5
    # ThreadsRatio >= 10
    # Contribution >= 0.1
    # UpstreamCoverage >= 0.5

  # Upper Edge TE Interaction
    # Interaction = UpEdge
    # ExonNumber = 1
    # Chimeric Reads >= 5
    # ThreadsRatio >= 10

  # Exon Contained in TE
    # Interaction = EInside
    # Chimeric Reads >= 5
    # ThreadsRatio >= 10

  # Repeat Contained in TE
    # Interaction = RInside
    # Chimeric Reads >= 5
    # ThreadsRatio >= 10
    # DownThread >= 10

# Chimeric Read Support (total)
  # Number of supporting chimeric reads required (Total)
  # Variable Number: 3 or 1/20 RPM
   scREADS=max(as.numeric(STDIN[4]), round( EXONIC_READS / 20000000 ) ) # >=

# Threads Ratio
  # DownThread / UpThread, set to value if UpThread =0
   scTHREAD=as.numeric(STDIN[5]) # >=
   scDownThread=as.numeric(STDIN[6]) # >=

# RPKM Cut-off to consider an exon 'expressed'
   scRPKM=as.numeric(STDIN[7]) # >=

# Contribution Ratio: Ratio of peak coverage between Exon/TE
   scCONTR=as.numeric(STDIN[8]) # >=

# Upstream Repeat Coverage Ratio to remove readthroughs
  # Peak TE coverage / Peak Upstream Coverage
   scUPCOV=as.numeric(STDIN[9]) # >=

# Upstream Exon Expression Ratio
   scUPEXON=as.numeric(STDIN[10]) # >=

# Splice Partner Classification
# Repeat Rank > 0 AND RepeatExonic

# IMPORT ======================================================================

# Read output from ChimericReadTool.sh
ChimeraTable = read.csv(file=INPUT,
                        header=T,
                        sep='\t')

# CALCULATED VALUES ===========================================================

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



# CHIMERA FILTRATION ==========================================================

# Global Filters (Apply to all cases)
  # Initialization
    print('Applying Global Filters to Chimeric Interaction Table')
    print(paste('     Input Entries   : ', nrow(ChimeraTable)))

# CASE 1: Upstream Interactions
  RETAIN_UP = which( ChimeraTable$ER_Interaction == 'Up' &
                       ChimeraTable$Total >= scREADS &
                       ChimeraTable$ThreadRatio >= scTHREAD &
                       ChimeraTable$Contribution >= scCONTR &
                       ChimeraTable$UpCov >= scUPCOV)

# CASE 2: Up Edge Interactions
  RETAIN_UPEDGE = which( ChimeraTable$ER_Interaction == 'UpEdge' &
                           ChimeraTable$Total >= scREADS &
                           ChimeraTable$exonRankInTranscript == 1 &
                           ChimeraTable$ExonInGene == 2 &
                           ChimeraTable$Contribution >= scCONTR/5 &
                           # Order of magnitude reduction since more direct
                           # evidence of promoter is avaiable
                           ChimeraTable$ThreadRatio >= scTHREAD )

# CASE 3: Exon Inside Interactions
  RETAIN_EINSIDE = which( ChimeraTable$ER_Interaction == 'EInside' &
                            ChimeraTable$Total >= scREADS &
                            ChimeraTable$ExonInGene == 2 &
                            ChimeraTable$ThreadRatio >= scTHREAD )

# CASE 4: Repeat Inside Interaction
  RETAIN_RINSIDE = which( ChimeraTable$ER_Interaction == 'RInside' &
                          ChimeraTable$Total >= scREADS &
			  ChimeraTable$ExonInGene == 2 &
			  ChimeraTable$ThreadRatio >= scTHREAD &
			  ChimeraTable$DownThread >= scDownThread )


  ChimeraOut = ChimeraTable[c(RETAIN_UPEDGE, RETAIN_EINSIDE, RETAIN_UP, RETAIN_RINSIDE),]
                      
# PARSE and OUTPUT ============================================================

    print(paste('     Output Entries   : ', nrow(ChimeraOut)))
  # Sort Rows
    SORT=4 # Coordinates
    ChimeraOut = ChimeraOut[order(ChimeraOut[,SORT]),]

  # Start of Repeat (For comparisons)
    ChimeraOut$RepeatID = paste('chr',ChimeraOut$Chromosome,':',ChimeraOut$RStart,sep='')

  # Add input to last column
  #  (for comparing multiple lists)
  
  # LIBRARY = unlist(strsplit(INPUT,split = '/'))[1]
    LIBRARY = INPUT
    ChimeraOut = cbind(ChimeraOut, LIBRARY)
    
  # Write output CSV
    write.table(ChimeraOut,
                file = OUTPUT,
                quote = F,
                sep = '\t',
                row.names = F,
                col.names = T)
  

# End of Script :D
