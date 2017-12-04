# LIONS
## Detecting TE-initiated transcripts from paired-end RNAseq

#### [LIONS Manuscript on bioRxiv](https://www.biorxiv.org/content/early/2017/06/13/149864)

 LIONS is a bioinformatic analysis pipeline which brings together a  few pieces of software and some home-brewed scripts to annotate a paired-end RNAseq library against a reference TE annotation set (such as Repeat Masker)

 `East Lion` scripts a bam file input, re-aligns it to a genome,  builds an ab initio assembly using Tophat2. This assembly is then  proccessed and local read searches are done at the 5' ends to find  additional transcript start sites and quality control the 5' ends of the assembly. The output is a file-type <library>.lions which annotates the intersection between the assembly, a reference gene set and repeat set.

 `West Lion` scripts compiles different `.lions` files, groups them into biological catagories (i.e. Cancer vs. Normal or Treatment vs. Control) and compares and analyzes the data to create graphs and meaningful interpretation of the data.

### Installation

1. Download the [LIONS repo](https://github.com/ababaian/LIONS/archive/master.zip)

2. Install the dependencies for `LIONS`
	- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
	- [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/)
	- [Python3](https://www.python.org/downloads/)
	- [pip 9.0+](https://pypi.python.org/pypi/pip)
	- [pysam](http://pysam.readthedocs.io/en/latest/api.html)
	- [R](https://www.r-project.org/)
	- [samtools v0.1.18](https://sourceforge.net/projects/samtools/files/samtools/0.1.18/)
	- [bedtools v2.18+](http://bedtools.readthedocs.io/en/latest/)
	- [Java OpenJDK](http://openjdk.java.net/)

3. Initialize the 'Parameter Files' for your system for `LIONS`
	1. `$LIONS_PATH/controls/<system>.sysctrl`: System-specific variables
	2. `$LIONS_PATH/controls/<project>.ctrl`: Project-specific variables
	3. `$LIONS_PATH/controls/<input>.ctrl`: List of RNA-seq file inputs for project

4. Add Reference / Annotation files for `LIONS`

4. Populate the resource files:
	**NOTE:** UCSC files are downloaded from: [UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgTables)). There is an `example` folder with example of what files should look like. 

	1. In `$LIONS_PATH/resources/<genomeName>/genome/` add a <genomeName>.fa genome sequence file
	2. In `$LIONS_PATH/resources/<genomeName>/repeat/` UCSC annotation for RepeatMasker for <genomeName>
	3. (Optional) In `$LIONS_PATH/resources/<genomeName>/annotation/` UCSC annotation for protein-coding genes

5. Run the master `lions.sh` in bash: 
	```
	bash $LIONS_PATH/lions.sh <$LIONS_PATH/controls/parameter.ctrl>
	```

#### [LIONS Folder Map](./docs/DIR_MAP.md)

#### [LIONS Error Codes](./docs/ERROR_CODES.md)

If you have any questions please email me: [Artem Babaian](mailto:ababaian@bccrc.ca).
I'll do my best to respond and help get this working!

