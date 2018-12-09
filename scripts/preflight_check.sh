#!/bin/bash
# Usage: .lions.sh <parameter.ctrl (opt.)>
set -e
# ===================================================================
# LIONS analysis pipelinei -pre-run file check
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
echo '================= L I O N S Pre-Flight Checks ================='
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

# INITIALIZATION ===================================================
# Start-up script which checks all requisites are operational for
# LIONS to run. Also initializes the project space  
# *** WRITE LAST ***

# Read parameter file (imports run parameters)
if [ -z $1 ]
then
    echo " No parameter input file specified. Importing default file:"
    echo "      ./LIONS/controls/parameter.ctrl"
    export PARAMETER="controls/parameter.ctrl"
    echo ''
else
    echo " Import parameter file."
    echo "    Project Parameters: ./$1"
    export PARAMETER=$1
    echo ''
fi
    # Run parameter script
    source $PARAMETER # works in bash only

# Run Initialization Script
echo ' running initializeLIONS.sh'

    bash $SCRIPTS/Initialize/initializeLIONS.sh

echo ' initialization completed successfully.'
echo ''
