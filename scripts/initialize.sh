#!/bin/bash
# initialize.sh
#
# 1) Check neccesary bin/scripts exist and are accesible
# 2) Create project workspace
# 3) (optional/ not implemented) run toy data to ensure pipeline works

echo ""
echo " INITIALIZING LIONS ..."
echo " ===================================================="
echo "     Project Name: $PROJECT "
echo "     Genome: $INDEX "
echo "     System: $SYSTEM"
echo "     LIONS dir: $BASE"

#cd $BASE #Go to base folder
 
# FILE CHECK ----------------------------------------------

# Function: Checks $FILE exists + permissions, returns that it does
# or returns a does not exist error 2 and exits

FCHECK_rs='if [ -s $FILE -a -r $FILE ]; then echo "     $FILE found."; else echo "     $FILE not found (empty or non-readable)."; echo " ===== ERROR 2: MISSING REQUISITE FILE ===== "; exit 2; fi'

FCHECK_x='if [ -s $FILE -a -x $FILE ]; then echo "     $FILE found."; else echo "     $FILE not found (empty or non-executable)."; echo " ===== ERROR 2: MISSING REQUISITE FILE ===== "; exit 2; fi'

FILE='' # File to check

# Resource Check

eval $FCHECK_rs

# PROJECT WORKSPACE ---------------------------------------

echo " ---------- PROJECT WORKSPACE ---------- "

# Make project folder 
#	mkdir -p projects/sandbox

# Data Demo ------------------------------------------
# NOT IMPLEMENTED ***
##
