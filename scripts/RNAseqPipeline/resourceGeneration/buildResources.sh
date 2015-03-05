#!/bin/bash

for s in "mm9" "hg18" "hg19"
do
	name=`echo $s"rs_0313"`
	echo $name
	mkdir -p $name
	cd $name
	sh /home/jlever/bin/RefSeqGenerate.sh $s $name
	cd -
done
