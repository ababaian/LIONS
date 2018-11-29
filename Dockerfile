#Download base image Centos 7
FROM centos:7

# Update Centos
#RUN yum update
RUN yum -y install yum-utils make wget xz-devel unzip epel-release gcc-gfortran libXt-devel libcurl-devel vim less  && \
    yum clean all
RUN yum-builddep -y python && \
    yum clean all

# Install pnglib for wigToBigWig to run
RUN yum -y install libpng12.x86_64 libpng12-devel.x86_64 && \
    yum clean all

# Install Python 3 and pysam
RUN wget https://www.python.org/ftp/python/3.5.0/Python-3.5.0.tgz && \
    tar xf Python-3.5.0.tgz && \
    cd Python-3.5.0 && \
    ./configure && \
    make && \
    make install && \
    cd .. && \
    rm -vr Python-3.5.0 Python-3.5.0.tgz 

RUN pip3 install --upgrade pip && \
    pip3 install pysam

# Install Bowtie2
RUN wget https://excellmedia.dl.sourceforge.net/project/bowtie-bio/bowtie2/2.3.4.1/bowtie2-2.3.4.1-linux-x86_64.zip && \
    unzip bowtie2-2.3.4.1-linux-x86_64.zip && \
    cp bowtie2-2.3.4.1-linux-x86_64/bowtie2* /usr/local/bin/. && \
    rm -rv bowtie2-2.3.4.1-linux-x86_64 bowtie2-2.3.4.1-linux-x86_64.zip

# Install Tophat2
run wget http://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz && \
    tar -zvxf tophat-2.1.1.Linux_x86_64.tar.gz && \
    cd tophat-2.1.1.Linux_x86_64 &&\
    find . -perm -111 -type f -exec cp '{}' /usr/local/bin/. \; && \
    cd .. && \
    rm -rv tophat-2.1.1.Linux_x86_64 tophat-2.1.1.Linux_x86_64.tar.gz

# Install OpenJDK 9
RUN yum -y install java-1.8.0-openjdk-headless java-1.8.0-openjdk-devel && \
    yum clean all

# Install Samtools
#RUN wget  https://github.com/samtools/htslib/releases/download/1.8/htslib-1.8.tar.bz2 https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2 && \
#    tar -jxf samtools-1.8.tar.bz2 && \
#    tar -jxf htslib-1.8.tar.bz2 && \
#    cd htslib-1.8 && \
#   ./configure && \
#    make && \
#    make install && \
#    cd ../samtools-1.8 && \
#    ./configure && \
#    make && \
#    make install && \
#    cd ..

RUN wget https://excellmedia.dl.sourceforge.net/project/samtools/samtools/0.1.18/samtools-0.1.18.tar.bz2 && \
    tar -jxf samtools-0.1.18.tar.bz2 && \
    cd samtools-0.1.18 && \
    make && \
    find . -perm -111 -type f -exec cp '{}' /usr/local/bin/. \; && \
    cd .. && \
    rm -vr samtools-0.1.18 samtools-0.1.18.tar.bz2
    #ln /usr/local/bin/samtools /usr/local/bin/samtools_0.1.18

# Install R

RUN yum-builddep -y R && \
    yum clean all

RUN wget https://cloud.r-project.org/src/base/R-3/R-3.5.0.tar.gz && \
    tar -zxf R-3.5.0.tar.gz && \
    cd R-3.5.0 && \
    ./configure && \
    make && \
    make install && \
    cd .. && \
    rm -rv R-3.5.0.tar.gz R-3.5.0 

# Install Bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz && \
    tar -zxvf bedtools-2.25.0.tar.gz && \
    cd bedtools2 && \
    make && \
    cp bin/* /usr/local/bin/. && \
    cd .. && \
    rm -vr bedtools-2.25.0.tar.gz bedtools2

# Install Cufflinks
RUN wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz && \
    tar -zxf cufflinks-2.2.1.Linux_x86_64.tar.gz && \
    cd cufflinks-2.2.1.Linux_x86_64 && \
    find . -perm -111 -type f -exec cp '{}' /usr/local/bin/. \; && \
    cd .. && \
    rm -vr cufflinks-2.2.1.Linux_x86_64 cufflinks-2.2.1.Linux_x86_64.tar.gz 


# Get LIONS
RUN wget https://github.com/ababaian/LIONS/archive/master.zip && \
    unzip master.zip && \
    rm -v master.zip && \
    mv /LIONS-master /LIONS-docker 
    
WORKDIR /LIONS-docker
