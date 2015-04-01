#!/bin/bash
# initializeBin.sh
# -----------------------------------------------
# source setBinPath.sh
# 
# Define the PATHS for software required by LIONS
# also creates shortcuts to software in ~/LIONS/bin

# Set Software list----------------------------
# For all software below, enter in where LIONS should access it
# Samtools (Less than version 1.x)
	export SAMTOOLS='samtools_0.1.18'
	export BAM2FASTX='bam2fastx'
	export TOPHAT2='tophat2'
	export BOWTIE2='bowtie2'
	export BOWTIE_BUILD='bowtie-build'
	export CUFFLINKS='cufflinks'
	export PYTHON3='python3'
	export JAVA='java'
	export RSCRIPT='Rscript'
	

# Functions -------------------------------------

# Checks if running the command set in $WARE exists
# If it doesn't exits with error 4: bin missing

wareExists='if [ -z $(command -v $WARE) ]; then echo ... $WARE not found. Exiting with status 4; exit 4; else echo ... $WARE found.; fi'

mkLink='ln -fs $(command -v $WARE) $BASE/bin/$WARE'
# shortcut link doesn't work. Can't figure out why not

py3Exists='$PYTHON3 -c "import $WARE"; if [ $? = 0 ]; then echo " ... $WARE python3 module found"; else echo " $WARE python3 module not found! ERROR 6"; echo " (note: If this is pysam, install it from ~/LIONS/software/pysam.zip"; exit 6; fi'

# Check Software list -------------------------
# 
echo "     Check that system software requirements are available and working."

echo ''
	# Samtools
	WARE=$SAMTOOLS; eval $wareExists #; eval $mkLink
		ln -fs $(command -v $WARE) $BASE/bin/samtools

	# Bam2Fastx
	WARE=$BAM2FASTX; eval $wareExists
		ln -fs $(command -v $WARE) $BASE/bin/bam2fastx

	# Tophat2
	WARE=$TOPHAT2; eval $wareExists
		ln -fs $(command -v $WARE) $BASE/bin/tophat2

	# Bowtie2
	WARE=$BOWTIE2; eval $wareExists
		ln -fs $(command -v $WARE) $BASE/bin/bowtie2

	# Bowtie-build
	WARE=$BOWTIE_BUILD; eval $wareExists
		ln -fs $(command -v $WARE) $BASE/bin/bowtie-build

	# Cufflinks
	WARE=$CUFFLINKS; eval $wareExists
		ln -fs $(command -v $WARE) $BASE/bin/cufflinks

	# Python3
	WARE=$PYTHON3; eval $wareExists
		ln -fs $(command -v $WARE) $BASE/bin/python3

	# Java
	WARE=$JAVA; eval $wareExists
		ln -fs $(command -v $WARE) $BASE/bin/java

	# Rscript
	WARE=$RSCRIPT; eval $wareExists
		ln -fs $(command -v $WARE) $BASE/bin/Rscript

echo ''

# Check that Python3 Modules are installed and function

echo "     Check Python3 Modules are installed"

	# Run 'python3 -c 'import $WARE' and check if the module $WARE exists in python3 else exit 6


	WARE='csv'; eval $py3Exists
	WARE='setuptools'; eval $py3Exists
	WARE='pysam'; eval $py3Exists
	WARE='sys'; eval $py3Exists
	WARE='pickle'; eval $py3Exists
	WARE='os'; eval $py3Exists
	WARE='pprint'; eval $py3Exists
	WARE='threading'; eval $py3Exists
	WARE='collections'; eval $py3Exists
	WARE='multiprocessing'; eval $py3Exists
	WARE='datetime'; eval $py3Exists
	WARE='re'; eval $py3Exists


# End of script : )
