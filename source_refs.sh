#! /bin/bash -e


# ===================================================================
# LIONS analysis pipeline - commandline arguments
# ===================================================================
#
# To download known (currently hg19 or hg38) references and associated 
# resources for the appropriate online repositories into directoryies 
# containing the required sub-directory structure
#
# Details can be found in README
#
echo ''
echo ''
echo '==============================================================='
echo '============== L I O N S Resource Download ===================='
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


if [ $# -ne 1 ]
then
  cat << EOF

incorrect number of parameters.

usage:
  source_refs.sh <reference>

  where <reference> is either hg38 or hg19

EOF
   exit 1
fi

wgt=`command -v wget`
crl=`command -v curl`

ref=$1

cmd=""

if [ wgt ]
then 
  cmd="$wgt"
elif [ crl ]
then
  cmd="$crl -O"
else
  cat << EOF

unable to find wget or curl

EOF
exit 1
fi


if [ $ref == "hg38" ]
then 
  # genomeName = hg38 ========================
  mkdir -p hg38/{genome,annotation,repeat}
  cd hg38/genome
  # Genome
  echo "Downloading Reference: hg38.fa.gz"
  $cmd http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
  # Gene Annotation (RefSeq)
  cd ../annotation
  echo "Downloading RefGene resource: refseq_hg38.ucsc.gz"
  $cmd https://s3-us-west-2.amazonaws.com/lionproject/resources/hg38/refseq_hg38.ucsc.gz
  gunzip -v refseq_hg38.ucsc.gz
  # Repeat Masker
  cd ../repeat
  echo "Downloading Repeat Masker resource: rm_hg38.ucsc.gz"
  $cmd https://s3-us-west-2.amazonaws.com/lionproject/resources/hg38/rm_hg38.ucsc.gz
  gunzip -v rm_hg38.ucsc.gz
  echo "hg38 successfully downloaded"

elif [ $ref == "hg19" ]
then 
  # genomeName = hg19 ========================
  mkdir -p hg19/{genome,annotation,repeat}
  cd hg19/genome
  # Genome
  echo "Downloading Reference: hg19.2bit"
  $cmd http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
  # Gene Annotation (RefSeq)
  cd ../annotation
  echo "Downloading RefGene resource: refseq_hg19.ucsc"
  $cmd https://s3-us-west-2.amazonaws.com/lionproject/resources/hg19/refSeq_hg19.ucsc.zip
  unzip refSeq_hg19.ucsc.zip
  rm -v refSeq_hg19.ucsc.zip
  # Repeat Masker
  cd ../repeat
  echo "Downloading Repeat Masker resource: rm_hg19.ucsc"
  $cmd https://s3-us-west-2.amazonaws.com/lionproject/resources/hg19/rm_hg19.ucsc.zip
  unzip rm_hg19.ucsc.zip
  rm -v rm_hg19.ucsc.zip
  echo "hg19 successfully downloaded"

else
  cat << EOF

$1 not a known reference.

Please note, this script currently only recognises human references hg19 and hg38.

EOF

fi

