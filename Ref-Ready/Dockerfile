# use Lions image as base
#FROM lions
FROM 642dfca95d52 

# Please note this Dockerfile assumes that the user has already created the Lions docker container locally


# Download hg19 references

RUN cd /LIONS-docker/resources

RUN /LIONS-docker/source_refs.sh hg19

# Download example files

#RUN cd /LIONS-data



WORKDIR /LIONS-docker
