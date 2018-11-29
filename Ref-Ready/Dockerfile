# use Lions image as base
#FROM lions
FROM 642dfca95d52 

# Please note this Dockerfile assumes that the user has already created the Lions docker container locally


# Download hg19 references

RUN cd /LIONS-docker/resources && \
  /LIONS-docker/source_refs.sh hg19

# Download example files

RUN mkdir /LIONS-data && \
  cd /LIONS-data

RUN echo $'kmh2	/LIONS-data/kmh1.bam	1\n\
kmh3	/LIONS-data/kmh1.bam	2\n '\
> /LIONS-docker/controls/input.list

WORKDIR /LIONS-docker
