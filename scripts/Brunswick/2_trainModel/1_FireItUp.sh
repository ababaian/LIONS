#!/bin/bash
#
# Fire it up!
# train some neural networks!
#
# This script will run on a genesis node
# and call the matches.sh multiple times to
# script to run on genesis

# CONTROL PANEL ============================================

# Directory of Brunswick scripts: "2_trainModel" folder

	BRUNSWICK='/projects/magerlab/ababaian/Brunswick/2_trainModel'

# Training Data (relative path from 2_trainModel folder)

	TRAINING_DATA='DATA_train.Rdata' #


# QSUB  Commands
# Leave blank if you're not using qsub
	QSUB=''
	#QSUB="qsub -S /bin/bash -V -q centos5.q -m e -M ababaian@bccrc.ca -pe ncpus 2 -l mem_free=2G -l mem_token=2G"

# Number of nodes to request
# set to 1 if not using cluster
	NNODES='1'	
	#NNODES='10'	
	
# ==========================================================


# Directory of Brunswick script base

	cd $BRUNSWICK # Apollo

	#cd /scratch/magerlab/ababaian/chimNeuro # Genesis Cluster

# Make directories to run ChimNeuro into
	mkdir -p UpEdge
	mkdir -p Up
	mkdir -p EInside
	mkdir -p RInside

# Set the for loop for how many bots you want
# in the neural network training
for x in $(seq 1 $NNODES )
do
	echo "Starting Job $x"
	
	# QSUB matches.sh <Interaction> <Directory>
	$QSUB matches.sh UpEdge $BRUNSWICK $TRAINING_DATA
	$QSUB matches.sh Up $BRUNSWICK $TRAINING_DATA
	$QSUB matches.sh EInside $BRUNSWICK $TRAINING_DATA
	$QSUB matches.sh RInside $BRUNSWICK $TRAINING_DATA


done

