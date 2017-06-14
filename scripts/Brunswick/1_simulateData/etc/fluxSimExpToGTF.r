# fluxSimExpToGTF.r
#
# Input .pro file from Flux Simulator
# and filber by expression level (FPKM)
# Return a subset GTF file of 'expressed transcripts only'


# FPKM = ( 10^9 * Fragment_Number ) / ( Exon_Len   * Reads_in_Library) 
#

# PRO Table description
  # 1  Locus	chrom:start-end[W|C]	identifier of the transcriptional locus, given by the chromosome (chrom), start respectively end position, and the strand (Watson or Crick).
  # 2	Transcript_ID	String	transcript identifier from the reference annotation.
  # 3	Coding	[CDS|NC]	specifies whether the transcript has an annotated coding sequence (CDS) or not (NC)
  # 4	Length	Integer	the mature length of the transcript after splicing out introns, disregarding the poly-A tail, as annotated in the reference annotation
  # 5	Expressed Fraction	Float	fraction of RNA molecules that represent transcripts that are qualitatively equal to this RNA form
  # 6	Expressed Number	Integer	absolute number of expressed RNA molecules
  # 7	Library Fraction	Float	fraction of cDNA molecules in the final library that have been produced from this transcript
  # 8 Library Number	Integer	absolute number of cDNA fragments generated from this transcript
  # 9	Sequenced Fraction	Float	fraction of total reads that have been sequenced from this transcript
  # 10	Sequenced Number	Integer	absolute number of reads sequenced from this transcript
  # 11	Covered Fraction	Float	fraction of the transcript that is covered by reads
  # 12	Chi Square	Integer	chi-square goodness of fit measurement of coverage uniformity
  # 13	Coefficient of Variation	Float	coefficient of variation for transcript coverage


# Control Panel =====================================================

# Input Flux Sim Pro file
  inPRO = 'fluxSimFiles/h1_50m.pro'

# Total Transcriptome GTF
  inGTF = 'h1Transcriptome.gtf'

# Output Expressed Transccriptome GTF
  outGTF = 'h1_50_expressed.gtf'

# IMPORT PRO TABLE ==================================================

# Import .pro file from FluxSim run. This contains all the data on the simulation
SIMULATION = read.table(file = inPRO, header = F, sep = '\t')
  colnames(SIMULATION) = c('Locus', 'TranscriptID', 'Coding,NC', 'MatureLength',
                           'ExpressedFrac','ExpressedNum','LibFrac','LibNum',
                           'SeqFrac','SeqNum','Coverage','Chi2','Variation')

TotalReads = sum(SIMULATION$SeqNum)

# Re-arranged FPKM formula for theoretical simualtion output
  # to prevent numeric overflow. Numerically the same.
  FPKM = (10^9/TotalReads) * (SIMULATION$SeqNum / SIMULATION$MatureLength)

# Expressed Gene Set (FPKM > 1)
  Expressed = which(FPKM >= 1)

# Expressed TranscriptID List
  ExpTranID = as.character(SIMULATION$TranscriptID[Expressed])

# IMPORT TRANSCRIPTOME ==============================================
# Load Transcriptome GTF File (For the simulation)
# GENCODE based
  TRANSCRIPTOME = read.csv(file = inGTF, header = F, sep = '\t')
  colnames(TRANSCRIPTOME) = c('chr','source','feature','start','end',
                              'score','strand','frame','attribute') # GTF headers

# All Attributes
  ATTR = as.data.frame(do.call(rbind, strsplit(as.character(TRANSCRIPTOME$attribute),split = ';')))
    # From Gencode; use these column names. There are more but discard them
    # (10 columns)
    colnames(ATTR) = c('gene_id', 'transcript_id', 'gene_type', 'gene_status', 'gene_name',
                       'transcript_type', 'transcript_status', 'transcript_name',
                       'exon_number', 'exon_id')

# Go through the attribute matrix and reduce it to the data
  ATTR = as.matrix(ATTR)[,1:10]
  
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

# Parse

# Transcriptome TRANID Vector
  gtfTranID = ATTR[,2]

# SUBSET TRANSCRIPTOME  ==============================================

Hits = c()

# Populate Recurrence Matrix
  for (Key in ExpTranID) { # for every column (Key Library)
    Query = ATTR[,2]
    Hits = c(Hits, which(Key == Query))
  }

# Expressed Transcritome Subset
  expressedTranscriptome = TRANSCRIPTOME[Hits,]

# Write the GTF File
  write.table(expressedTranscriptome, file = outGTF,
              col.names = F, row.names = F, quote = F, sep = '\t')


# Write Table of Expressed Transcript ID Output
# write.table(ExpTranID,file = 'h1_100_ExpressedTranscripts.txt',col.names = F, row.names = F, quote = F)

# Plot of FPKM distribution
#   library(ggplot2)
#   ggplot(data.frame(FPKM), aes(x = log(FPKM) ) ) + geom_density(kernel = "gaussian") + scale_x_continuous(limits = c(-10,10)) #+ scale_x_log10()
