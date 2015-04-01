#!/bin/bash
# initializeScripts.sh
# -----------------------------------------------
# source initializeScripts.sh
# 
# Checks that all the scripts for LIONS are in place and read/executable

# Script list----------------------------

# LIONS initialization
	# initializeScripts.sh
	# initializeBin.sh
	# initializeRes.sh

# eastLion.sh
	# runTH2.sh

# westLion.sh
	# ./scripts/RNAseqPipeline
	#  RNAseqMaster.sh
	#  RNAseqCoverageCalculator.sh
	#  RNAgetRes.sh
	#  RNAseqPipeline.sh
	#  RPKM.sh
	#  TE_stats.sh
	#  WIG_RPKM.sh
	#  buildResourceGTF.sh
	#  buildResources.sh
	#  EnsembleResourceGenerator.sh
	#  GencodeGenerate.sh
	#  RefSeqGenerate.sh
	#  RepeatMaskerGenerate.sh

	# ./scripts/ChimericReadTool
	#  ChimericReadTool.sh
	#  exon_scan.sh
	#  intervalTree.py
	#  chimericReadSearch.py


# Functions -------------------------------------

# Test if script exist and can be read

FCHECK_rs='if [ -s $FILE -a -r $FILE ]; then echo "... $FILE found."; else echo "     $FILE not found (empty or non-readable)."; echo " ===== ERROR 2: MISSING REQUISITE FILE ===== "; exit 2; fi'

FILE=''

# Check Script list -------------------------
echo "     Check that LIONS scripts exist and are read/executable"

echo ''

# Initializations
	FILE='Initialize/initializeLIONS.sh'; eval $FCHECK_rs
	FILE='Initialize/initializeScripts.sh'; eval $FCHECK_rs
	FILE='Initialize/initializeBin.sh'; eval $FCHECK_rs
	FILE='Initialize/initializeRes.sh'; eval $FCHECK_rs

# East Lion
	FILE='eastLion.sh'; eval $FCHECK_rs

# West Lion
	FILE='westLion.sh'; eval $FCHECK_rs

# RNAseq Pipeline
	FILE='RNAseqPipeline/RNAseqMaster.sh'; eval $FCHECK_rs
	FILE='RNAseqPipeline/RNAseqCoverageCalculator.sh'; eval $FCHECK_rs
	FILE='RNAseqPipeline/RNAgetRes.sh'; eval $FCHECK_rs
	FILE='RNAseqPipeline/RNAseqPipeline.sh'; eval $FCHECK_rs
	FILE='RNAseqPipeline/RPKM.sh'; eval $FCHECK_rs
	FILE='RNAseqPipeline/TE_stats.sh'; eval $FCHECK_rs
	FILE='RNAseqPipeline/WIG_RPKM.sh'; eval $FCHECK_rs

	FILE='RNAseqPipeline/resourceGeneration/buildResources.sh'; eval $FCHECK_rs
	FILE='RNAseqPipeline/resourceGeneration/buildResourceGTF.sh'; eval $FCHECK_rs
	FILE='RNAseqPipeline/resourceGeneration/RepeatMaskerGenerate.sh'; eval $FCHECK_rs
	FILE='RNAseqPipeline/resourceGeneration/RefSeqGenerate.sh'; eval $FCHECK_rs

# Chimeric Read Tool
	FILE='ChimericReadTool/ChimericReadTool.sh'; eval $FCHECK_rs
	FILE='ChimericReadTool/exon_scan.sh'; eval $FCHECK_rs
	FILE='ChimericReadTool/intervalTree.py'; eval $FCHECK_rs
	FILE='ChimericReadTool/chimericReadSearch.py'; eval $FCHECK_rs
	FILE='ChimericReadTool/chimIntersect.sh'; eval $FCHECK_rs
	FILE='ChimericReadTool/chimIntLookup.R'; eval $FCHECK_rs
	FILE='ChimericReadTool/chimSort.R'; eval $FCHECK_rs

# Chimeric Analysis Scripts
	FILE='ChimericAnalysis/chimGroup.R'; eval $FCHECK_rs

echo ''
# End of script : )
