cat ANN.list | while read line; do Rscript ANN_training.r $line; done
