cat input.list | while read line; do Rscript fluxSimExpToGTF.r $line; done
