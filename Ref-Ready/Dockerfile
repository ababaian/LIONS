# use Lions image as base
FROM lions

# Please note this Dockerfile assumes that the user has already created the Lions docker container locally


# Download hg19 references

RUN cd /LIONS-docker/resources && \
  /LIONS-docker/source_refs.sh hg19

# Download example files

RUN mkdir /LIONS-data && \
  cd /LIONS-data && \
  wget https://www.encodeproject.org/files/ENCFF001RFH/@@download/ENCFF001RFH.fastq.gz -O GM12878.rep1.R1.fastq.gz && \
  wget https://www.encodeproject.org/files/ENCFF001RFG/@@download/ENCFF001RFG.fastq.gz -O GM12878.rep1.R2.fastq.gz && \
  wget https://www.encodeproject.org/files/ENCFF001RED/@@download/ENCFF001RED.fastq.gz -O K562.rep1.R1.fastq.gz && \
  wget https://www.encodeproject.org/files/ENCFF001RDZ/@@download/ENCFF001RDZ.fastq.gz -O K562.rep1.R2.fastq.gz && \
  echo $'gm12878	/LIONS-data/GM12878.rep1.R1.fastq.gz,/LIONS-data/GM12878.rep1.R2.fastq.gz	1\nk562	/LIONS-data/K562.rep1.R1.fastq.gz,/LIONS-data/K562.rep1.R2.fastq.gz	2\n ' > /LIONS-docker/controls/input.list

WORKDIR /LIONS-docker

