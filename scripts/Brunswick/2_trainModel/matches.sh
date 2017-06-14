#!/bin/bash
#
# Go to working directory for chimNeuro
# and run an instance of the Rscript
#

INTERACTION=$1 # ER Interaction which will be trained

DIR=$2 # Brunswick directory 

TRAINING_DATA=$3

#cd /projects/magerlab/ababaian # Apollo 
#cd /scratch/magerlab/ababaian/chimNeuro # Genesis

Rscript chimNeuro.R $1 $2 $3 &
Rscript chimNeuro.R $1 $2 $3


