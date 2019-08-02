#!/bin/bash
# parameter.sh

set -e

# ===================================================================
# LIONS analysis pipeline - commandline arguments
# ===================================================================
#
# Analyze an input .bam RNAseq file for transcripts which initiate in
# transposable elements and create an annotation file. Compare TE
# files between biological groups.
#
# Details can be found in README
#
echo ''
echo ''
echo '==============================================================='
echo '========================= L I O N S ==========================='
echo '==============================================================='
echo '''                             _   _
                           _/ \|/ \_
                          /\\/   \//\
                          \|/<\ />\|/   *RAWR*
                          /\   _   /\  /
                          \|/\ Y /\|/
                           \/|v-v|\/
                            \/\_/\/
'''
echo ''

# import default parameters

# Read parameter file (imports run parameters)
    echo " Importing default file:"
    echo "      ./LIONS/controls/parameter.ctrl"
    export PARAMETER="./controls/parameter.ctrl"
    echo ''
# Run parameter script
    source $PARAMETER # works in bash only

usage='''Commandline syntax:

parameter.sh <options>

where:
  -h | --help                    Show this usage help
  -P | --parameter <parameter>   Define user parameter.ctrl file, N.B. This option must be defined first.
  -b | --base <path>             Path to working directory
  -p | --project <project_name>  Project identifier
  -i | --inputlist <file>        CSV File containing list of input files and sample data <libName> <libPath> <group> csv file
  -I | --index <INDEX>           Genome resource set to use (see README)
  --callsettings                 Pre-set TE-initiation classification settings
  --geneset <ucsc annotation>    UCSC-format annotation file for reference geneset
  --repeatmasker <ucsc>          UCSC-format annotation file for RepeatMasker
  --systemctrl                   System-specfic control file
  --bowtie                       Options to be passed to Bowtie2 during alignment
  --alignbypass                  Bypass alignment, "1" to bypass, "0" to carry out alignment 
  --inread                       Tophat Inner read distance
  --threads                      Tophat - Number of threads
  --denovo
  --gtfguide                     Guide file upon which to assemble
  --minfrag                      Minimum fragments per transfrag, default is 10
  --multifrag                    Multi-mapping transfrag fragments, default is 0.75;
  --mintrim                      Minimum coveraged to attempt 3` trimming, default is 10
  --trimdrop                     Minimum Trim dropoff in fraction
  --merger                       Merge Radius in bp, default is 50
  --quality                      RNAseq Analysis Quality - q = quality cutoff; F bam flags to discard, i.e. "q10.F772"
  --cggrouprecurrence            Recurrent Cutoff, Must be present in at least this many Group 2 libraries
  --cgspecificity                Specificity Cutoff, Must be absent in at least this many Group 1 libraires

Example:


'''

TEMP=`getopt -o b:p:i:hI:P: -l parameter:,base:,project:,inputlist:,threads:,callsettings:,index:,geneset:,repeatmasker:,systemctrl:,bowtie:,alignbypass:,inread:,denovo:,gtfguide:,minfrag:,multifrag:,mintrim:,trimdrop:,merger:,quality:,cggrouprecurrence:,cgspecificity: \
	-- "$@"`
#	-n 'parameter.sh' -- "$@"`
echo "$TEMP"
eval set -- "$TEMP"

if [ "$#" -eq "2" ]; then
    echo "Only one argument, $2";
    export PARAMETER="$2"
    source $2;
else
while true; do
  case "$1" in 
    -h | --help )
	echo "$usage"; exit ;;
    -P | --parameter )
# Users parameter.ctrl file
        echo "Multiple arguments, $2";
        export PARAMETER="$2"
        source $2;
        shift 2 ;;
    -b | --base )
# The LIONS home folder
	export BASE="$2"; # Base folder for ~/LIONS/
        shift 2 ;;
    -p | --project )
# PROJECT -------------------------------------------------
	export PROJECT="$2"; # Project Name
	shift 2 ;;
    -i | --inputlist )
	export INPUT_LIST="$2"; # <libName> <libPath> <group> csv file
        shift 2 ;;
	# <libPath> is one of
	#     A) Paired end bam file. '/home/libPath.bam'
	#     B) Sorted paired fastQ files. Comma seperated.
	#        '/home/lib.fq1,/home/lib.fq2'
    --callsettings )
	export CALLSETTINGS="$2"; shift 2 ;;
	# Pre-set TE-initiation classification settings
        # - 'oncoexapatation'  : Detect high abundance TE isoforms
        # - 'transcriptomeANN' : Artifical Neural Network based classifier
        # - 'screenTE'         : High sensitivity, low specificity detection
        # - 'driverTE'         : TE-initiations as main drivers of gene-expression
        # = 'custom'           : Use custom settings defined below

# RESOURCES 
    -I | --index )
	export INDEX="$2"; # Genome resource set to use
	shift 2 ;;
	#  Should be a ~/LIONS/resources/<INDEX>/ folder
	#  See README file 'RESOURCE INITIALIZATION' to set up
	#  the folder prior to running LIONS the first time
    --geneset )
	export GENESET="$2"; # Reference Gene Set. UCSC format
        shift 2 ;;
	# in ~/LIONS/resources/<INDEX>/annotation/<GENESET>	
    --repeatmasker )
	export REPEATMASKER="$2"; # RepeatMasker. UCSC format
	shift 2 ;;
	# in ~/LIONS/resources/<INDEX>/repeat/<REPEATMASKER>

# SYSTEM --------------------------------------------------
	# For each system you'd like to run LIONS on you can edit
	# one <system>.ctrl file
    --systemctrl )
	export SYSTEMCTRL="$2";
	source $SYSTEMCTRL;
	shift 2 ;;


# EAST LION Parameters-----------------------------------------------
# Bowtie 2
	# Additional alignment parameters passed on to Tophat2
    --bowtie )
	export BOWTIE="$2"; shift 2 ;;

# Tophat 2
	# Alignment Bypass
	# 0: Calculate a new alignment for the input
	# 1: Do not re-calculate the alignment, simply create
	#    symbolic link to the bam file in <input.list> within
	#    the LIONS folder architecture
    --alignbypass )
	export ALIGNBYPASS="$2"; shift 2 ;;

    --inread )
	# Inner Read distance [-r]
	# 200 - 75 - 75 = 50
	export INREAD="$2"; 
	export ctrlTH2=" $BOWTIE -p $THREADS -r $INREAD --report-secondary-alignments";
	shift 2 ;;

    --threads )
	export THREADS="$2";
	export ctrlTH2=" $BOWTIE -p $THREADS -r $INREAD --report-secondary-alignments";
	shift 2 ;;

# Cufflinks
	# Ab initio (guided) or De Novo assembly?
    --denovo )
	export deNovo="$2"; # 1 = De Novo; 0 = Ab Initio
        shift 2 ;;

    --minfrags )
	#Minimum fragments per transfrag
	# default is 10
	export minFrag="--min-frags-per-transfrag $2"; shift 2 ;;

    --multiflag )
	# Multi-mapping transfrag fragments
	# Default is 0.75;
	export multiFrag="--max-multiread-fraction $2"; shift 2 ;;

    --mintrim )
	# Minimum coveraged to attempt 3` trimming
	# Default is 10
	export minTrim="--trim-3-avgcov-thresh $2"; shift 2 ;;

    --trimdrop )
	# Minimum Trim dropoff in fraction
	export trimDrop="--trim-3-dropoff-frac=$2"; shift 2 ;;

    --merger )
	# Merge Radius in bp
	# Default is 50 bp
	export mergeR="--overlap-radius $2"; shift 2 ;;

	# Multiple Mapping read correction
	# include [-u]
	

# RNAseqPipeline

	# RNAseq Analysis Quality
		# q = quality cutoff; F bam flags to discard
		# i.e. 'q10.F772'
		# i.e. 'q1.F1796' 
    --quality )
	export QUALITY="$2"; shift 2 ;; 

    -- ) shift ; break ;;
    * ) break ;;

  esac
done


	# Check that assembly guide GTF file exists if it's being used
	if [ $deNovo == '0' ]
	then
		if [ -z $guide ]
		then # Guide file upon which to assemble
		  export guide="$RESOURCES/annotation/guide.gtf";
		fi

		# Check of Guide file doesn't exists; exit
		if [ ! -s $guide ]
		then
		echo " ERROR 5 - Missing resource file"
		echo " Ab initio assembly selected (using a GTF guide -g)"
		echo " but the guide file $guide does not exist or cannot be read"
		echo " check parameter.ctrl and/or initialize a guide file"
		exit 5
		fi
	fi


	if [ $deNovo == '1' ]
	then # De Novo Assembly
		export ctrlCL2="-u -p $THREADS $minFrag $multiFrag $minTrim $trimDrop $mergeR"
	else # Ab Initio Assembly
		export ctrlCL2="-u -p $THREADS -g $guide"
	fi

# Chimeric Read Tool
	# self-contained

# ChimSort

	# Sort Bypass: If a sorted .lions file is already present
	# do you want to re-calcualate the 'initiations'?
	# 0 = Re-calcualte 'initiation' using parameters set below
	# 1 = Don't re-calculate 'initiation' if it exists
	export SORTBYPASS='0'

	# These parameters define what is defined as a 
	# "TE-initiated transcript" and what is exonization/termination are
	# note: See chimSort.R script to see how these are used for sorting
        # in more detail

if [ $CALLSETTINGS == 'custom' ]
then
	# Define call settings manually

	# Number of Chimeric Reads Required (total)
	export crtReads='3' # Variable: <Value below> or 1/20 RPM
	
	# Thread Ratio
	export crtThread='10' # >=
	export crtDownThread='10' # |number| required 

	# RPKM cut-off to consider an exon 'expressed'
	export crtRPKM='1' # >=

	# Contribution of TE to total transcript expression (Exon / TE)
	export crtContribution='0.1' # >=

	# Expression immediatly upstream TE
	export crtUpstreamCover='2' # >=

	# Expression of exons upstream of involved exon
	export crtUpstreamExonCover='1.5' # >=

	# Splice Partner Classification
	# Repeat Rank > 0 AND RepeatExonic
	
	# Final Custom Export String
	export CRT=$( echo $crtReads $crtThread $crtDownThread $crtRPKM $crtContribution $crtUpstreamCover $crtUpstreamExonCover )

elif [ $CALLSETTINGS == 'oncoexaptation' ]
then
	# Oncoexaptation Export String
	export CRT=$( echo '3 10 10 1 0.1 2 1.5' )

elif [ $CALLSETTINGS == 'screenTE' ]
then
	# Screen TE-initiations Export String
	export CRT=$( echo '2 5 5 1 0.05 2 1' )

elif [ $CALLSETTINGS == 'driverTE' ]
then
	# driverTE export string. Same as Oncoexaptation but a higher
	# contribution score (0.75 vs 0.1) required as  the 'driver'
	# of expression.
	export CRT=$( echo '5 10 10 1 0.75 2 1.5' )

elif [ $CALLSETTINGS == 'transcritomeANN' ]
then
	# ANN classifier
	# standard 'onco-exaptation' classification will be done
	# and the ANN classifier is 'on top' of that
	export CRT=$( echo '3 10 10 1 0.1 2 1.5' )

	# ANN model to use (stored in ~/LIONS/resources/ANN
	export ANNMODEL="$RESOURCES/ANN/ANN_Model_160919.Rdata"
	
fi


# WEST LION Parameters-----------------------------------------------
# chimGroup.R
	# These parameters define recurrent and specific
	# 'TE initiated' events
	# 
	# The third column of <input.list> defines libraries as one of
	# 1: Control / Normal
	# 2: Experimental / Cancer
	# 3: Other / unused

	# Recurrent Cutoff
	# Must be present in at least this many Group 2 libraries
	export cgGroupRecurrence='1'

	# Specificity Cutoff
	# Must be absent in at least this many Group 1 libraires
	export cgSpecificity='1'

	# chimGroup command string
	export CG=$(echo $cgGroupRecurrence $cgSpecificity)

fi
##########################################################
#                     Run analysis                       #
##########################################################

# Run Initialization Script
echo ' running initializeLIONS.sh'

	bash $SCRIPTS/Initialize/initializeLIONS.sh

echo ' initialization completed successfully.'
echo ''


# EAST LION =========================================================
echo ''
echo '                     E A S T       L I O N                     '
echo ''
echo ' ./LIONS/scripts/eastLion.sh ' 
echo '==============================================================='
echo '  Align reads to genome and perform TE-initiation analysis'
echo ''
cd $pDIR #./LIONS/projects/<projectName>

# Initialize Summit Log file
touch summitLog_$RUNID

# Loop through each library in input file
iterN=$(wc -l $INPUT_LIST | cut -f1 -d' ' -)

for nLib in $(seq $iterN)
do
	# Extract row of entries from input list
	rowN=$(sed -n "$nLib"p $INPUT_LIST)

        # test to ignore empty lines in $INPUT_LIST
        if [[ ! -z "$rowN" ]]
        then 
		# Library Name
		libName=$(echo $rowN | cut -f1 -d' ')

		echo " Iteration $nLib: $libName ------------------------------------------"
		echo "      run: $QSUB eastLion.sh $libName"

		if [ ! -e $pDIR/$libName/$libName.lions ]
		then
		# Lions output for this library doesn't exist
		# so let's make it.

			if [ $CLUSTER == '1' ]
			then # Cluster QSUB
				$QSUB $SCRIPTS/eastLion.sh $libName

			else # Local (no) QSUB
				$SCRIPTS/eastLion.sh $libName
			fi

		elif [ $SORTBYPASS = '0' ]
		then
		# Lions output already exists but
		# East Lion bypass is set to false, re-calculate lions file
		# so let's make it.

			if [ $CLUSTER == '1' ]
			then # Cluster QSUB
				$QSUB $SCRIPTS/eastLion.sh $libName

			else # Local (no) QSUB
				$SCRIPTS/eastLion.sh $libName
			fi

		else
		# East Lion file already exists and bypass is true (1)
		# Skip the east lion

			echo "   East Lions has previously been completed. "
			lionSuccess='1'
			echo $libName $lionSuccess $(date) >> $pDIR/summitLog_$RUNID
		fi

		echo " ... run complete -------------------------------------------"
		echo ''
		echo ''
	fi
done

# Check that all libraries have completed
	#iterN is the number of libraries
summitN=$(wc -l $pDIR/summitLog_$RUNID | cut -f1 -d' ' )

while [  $summitN -lt $iterN ]  # Not all EAST LION iterations have completed
do 

	# Verbose
	echo " $summitN / $iterN East Lion scripts completed. Waiting..."
	date

	# Wait 10 minutes
	sleep 600s # Actual

	# Recalculate summitN
	summitN=$(wc -l $pDIR/summitLog_$RUNID | cut -f1 -d' ')

done

# All runs completed

# Check if they are each succesful
for Log in $(cut -f2 $pDIR/summitLog_$RUNID)
do
	if [ $Log = '0' ];
	then
		echo " ERROR 15: One of the East Lion Libraries didn't finish."
		exit 15
	fi
done

# Clear summit log
rm $pDIR/summitLog_$RUNID

echo ''
echo ' All EAST LION scripts have completed. '


# WEST LION =========================================================
echo ''
echo '                     W E S T       L I O N                     '
echo ''
echo ' ./LIONS/scripts/westLion.sh ' 
echo '==============================================================='
echo '  Group and analyze lions files of TE-initiation events'
echo ''
cd $pDIR #./LIONS/projects/<projectName>

# Run West Lions Script
	echo ' run: westLion.sh'
	bash $SCRIPTS/westLion.sh

echo ''
echo ' WEST LION scripts have completed. '

