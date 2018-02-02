# LIONS output definitions

## Ouptut File-types
LIONS produces several outputs from different stages of the analysis apart from
the standard outputs one would expect (.bam / .gtf).

`<library>.lcsv` / `.pc.lcsv`
	These are LIONS CSV files; that is the raw calculations for all major
	numeric operations.
	This includes ALL TE-exon interactions types (Initiation, Exonization and
	Termination). As such there are usually hundreds of thousands of TEs which
	have read fragments joining them to some assembled exon.

	The `.pc.` pre-suffix means the data has been intersected to the input
	set of protein coding genes.

	Use this file for re-calculating "TE-Initiations" with new parameters.

`<library>.lion`
	This is the filtered set of TE-exon interactions which have been classified
	as "TE-Initiations" or TE transcription start sites. This is per-library
	input.

`<project>.lions` / 
	A merged file of several `.lion` files combining biological groups defined
	in the `input.list`. A good example of this is merging 10 cancer libraries
	and 10 normal libraries and outputing only those TE-initiations which are
	in at least 20% of Cancer and no Normal libraries. These parameters can be
	changed in the `paramter.ctrl` input.

`<project>.rslions`	
	The `rs` is for Recurrent and Specific TE-initiations only. That is if you
	compare the set of libraries 1 (Normal) vs set 2 (Cancer), this contains only
	those TE-initiations which occur multiple times in Cancer (recurrant) and do
	not occur in Normal (specific). As defined by `$cgGroupRecurrence` and
	`$cgSpecificity` in the `paramter.ctrl` file.

`<project>.inv.rslions`	
	The `.inv.` pre-suffix is simply the **inverse** of the `.rslion` file. So
	instead of "Cancer vs. Normal", "Normal vs. Cancer". A neccesary control if
	one makes any conclusions based on enrichment/depletion.


## Output Columns
Most columns should be self-explanatory, some are not.

transcriptID: Unique identifier for the transcript (isoform). Usually taken
	from the assembly/reference transcriptome

exonRankInTranscript: For each TE-exon interaction combination (row) which exon
	in the 'transcriptID' is this row referring to

repeatName: The <repeat_name>:<repeat_class>:<repeat_family> taken from input
	set

coordinates: Useful coordinates for visualizing the interaction. It starts/ends
	in the exon and repeat so when opening in a visualization tool you can see
	the reads spanning this area.

ER_Interaction: The type of relative intersection in the genome between the exon
	and the repeat. Definitions are relative to the exon. Can be "Up", "UpEdge",
	"EInside", "RInside", "Down", "DownEdge".

IsExonic: ??

ExonsOverlappingWithRepeat: A list of <transcriptID:exonRank> which overlap
	the repeat.

ER / DR / DE / DD / Total: A count of the number of TE-Exon sequence fragments
	which join this rows TE and Exon. ER means that one end overlaps the Exon
	and one end overlaps the Repeat exclusively, DD means that both ends of the
	fragment overlap both (dual) exon and repeat ...

Chromosome / EStart / EEnd / EStrand: Start, end and strand of the exon

RStart / REnd / RStrand: Start, end and strand of the repeat

RepeatRank: Relative exon/intron position of the repeat to the contig

UpExonStart / UpExonEnd: Coordinates used for calculating expression of
	genome immediatly adjacent an exon boundary. Useful for quantifying read-
	through or spurious transcriptional events.

UpThread: The number of read 'threads' going upstream of the exon.
	See Manuscript for a figure explaining this.

DownThread: The number of read 'threads' going downstream of the exon.
	See Manuscript for a figure explaining this.

ExonRPKM: RPKM calculation for this exon

ExonMax: The maximum coverage count reached within the exon boundaries. Often
	more reliable measure of expression then RPKM for small exons.

UpExonRPKM / UpExonMax: The expression of the exon immediatly upstream of the
	one this row is referring to. (i.e. Exon 1 expression if the row refers
	to Exon 2). Useful for quantifying the relative increase in expression
	when a TE is acting as an alternative promoter into a downstream exon.

RepeatRPKM / RepeatMaxCoverage:	Expression level within repeat boundaries.

UpstreamRepeatRPKM / UpstreamRepeatMaxCoverage: The expression adjacent to the
	repeat, a test for background expression levels.

RefID: When intersecting to a reference gene set, the gene symbol of any genes
	which intersect the area between the Exon-Repeat coordinates.

RefStrand: Strand of the reference genes defined above

assXref: The strand-relationship between the reference gene and the contig exon
	This accounts for anti-sense long non-coding RNA (as), or transcripts
	which run anti-sense to the reference gene. (s) is sense and (c) means
	complex, often some combination of multiple genes. (u) means it could not
	be determined.

Contribution: An estimate of the promoter contribution of this Repeat TSS to
	the expression of gene in total. Calculated with ExonMax and UpExonMax.

UpCov: Ratio of the coverage adjacent to an exon and the exon expression

UpEoxnRatio: Ratio of hte expression of the exon and it's upstream exon

ThreadRatio: DownThread / UpThread. Set to [10] if dividing by zero.

RepeatID: A unique Identifer for each Repeat in the genome (left-most
	coordinate). Can repeat and thus be used for determining one repeat
	inititating a trancsript in different assemblies.

LIBRARY: Library from which this repeat-exon interaction was calculated from.
