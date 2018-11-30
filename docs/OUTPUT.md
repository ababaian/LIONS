# LIONS output definitions

## Output File-types
LIONS produces several outputs from different stages of the analysis apart from the standard outputs one would expect (.bam / .gtf).

`<library>.lcsv` / `.pc.lcsv` These are LIONS CSV files; that is the raw calculations for all major numeric operations. This includes ALL TE-exon interactions types (Initiation, Exonization and Termination). As such there are usually hundreds of thousands of TEs which have read fragments joining them to some assembled exon.

 The `.pc.` pre-suffix means the data has been intersected to the input set of protein coding genes.
 Use this file for re-calculating "TE-Initiations" with new parameters.

`<library>.lion` This is the filtered set of TE-exon interactions which have been classified as "TE-Initiations" or TE transcription start sites. This is per-library input.

`<project>.lions` /  A merged file of several `.lion` files combining biological groups defined in the `input.list`. A good example of this is merging 10 cancer libraries and 10 normal libraries and outputting only those TE-initiations which are in at least 20% of Cancer and no Normal libraries. These parameters can be changed in the `parameter.ctrl` input.

`<project>.rslions`	 The `rs` is for Recurrent and Specific TE-initiations only. That is if you compare the set of libraries 1 (Normal) vs set 2 (Cancer), this contains only those TE-initiations which occur multiple times in Cancer (recurrant) and do not occur in Normal (specific). As defined by `$cgGroupRecurrence` and `$cgSpecificity` in the `parameter.ctrl` file.

`<project>.inv.rslions`	The `.inv.` pre-suffix is simply the **inverse** of the `.rslion` file. So instead of "Cancer vs. Normal", "Normal vs. Cancer". A necessary control if one makes any conclusions based on enrichment/depletion.


## Output Columns
For a simplified output see: `~/LIONS/scripts/lions2bed.sh` conversion script.

### `.lion` & `.lions`
Descriptions of each column in the `.lion(s)` output format. This format contains many data fields which can be used for downstream analysis. 

1. *transcriptID*: Unique identifier for the transcript (transfrag) of origin. Taken from the assembly or reference transcriptome

2. *exonRankInTranscript*: The exon # of the *transcriptID* for which this row (unique TE-Exon interaction) is referring

3. *repeatName:* The <repeat_name>:<repeat_class>:<repeat_family> taken from input set of TEs (RepeatMasker UCSC)

4. *coordinates:* Coordinates for visualizing the TE-Exon interaction. These coordinates connect the transcript-exon and TE for which this row refers

5. *ER_Interaction:* The type of relative intersection in the genome between the exon and the TE, defined relative to the exon. Can be "Up", "UpEdge", "EInside", "RInside", "Down", "DownEdge"

6. *IsExonic:* <Yes/No> Does the TE overlap with a known exon?

7. *ExonsOverlappingWithRepeat:* A list of <transcriptID:exonRank> which overlap the TE

8. *ER / DR / DE / DD / Total:* A count of the number of TE-Exon read fragments joining this TE and Exon. E(xon), R(epeat) and D(ual) relate to how each of the two paired-end reads of a single fragment intersect the transcript-exon, TE or both, respectively

13. *Chromosome / EStart / EEnd / EStrand:* Start, end and strand of the transcript- exon

16. *RStart / REnd / RStrand:* Start, end and strand of the repeat/TE

20. *RepeatRank:* The relative exon/intron position of the repeat with respect to transcript-exon. -1 means the TE is upstream of the first exon. If <IsExonic> is `No` and the RepeatRank is 5, then the TE is in the fifth intron of the transcript and does not intersect exon 5

21. *UpExonStart / UpExonEnd:* Coordinates of the most immediate upstream exon of the transcript. If exonRankInTranscript = 8, then this is the coordinates of for the 7th exon in that transcript. If exonRank = 1, then it returns exon 1 coordinates

23. *UpThread / DownThread:* The number of read 'threads' going upstream of the exon. See User Manual {ref} for a figure explaining this.

25. *ExonInGene:* Deprecated

26. *ExonRPKM:* RPKM calculation for the transcript-exon only (not the entire transcript)

27. *ExonMax:* The maximum coverage count reached within the exon boundaries. Often more reliable measure of expression then RPKM for small exons.

28. *UpExonRPKM / UpExonMax:* The expression of the exon immediately upstream of the one this row is referring to. (i.e. Exon 1 expression if the row refers to Exon 2). Useful for quantifying the relative increase in expression when a TE is acting as an alternative promoter into a downstream exon.

30. *RepeatRPKM / RepeatMaxCoverage:* Expression within TE boundaries.

32. *UpstreamRepeatRPKM / UpstreamRepeatMaxCoverage:* The level of expression adjacent to the TE in the genome, a measure of background expression or spurious transcription at this locus

34. *RefID / RefStrand :* The gene symbol and strand for any known genes which intersect the area between the TE and the transcript-exon coordinates. Taken from the input reference gene set. Under the standard pipeline, this is restricted to protein-coding genes.

36. *assXref:* The strand-relationship transcript-exon and the reference gene. Can be sense (s), anti-sense (as), intergenic (i), complex (c) or unknown (u). Complex often refers to cases where multiple genes are present.

37. *Contribution:* An estimate of the promoter contribution of this TE-initiation to the expression of gene in total. Calculated as ExonMax/UpExonMax.

38. *UpCov:* Ratio of the coverage adjacent to the repeat and the repeat itself. Calculated as RepeatMax / UpstreamRepeatMax.

39. *UpExonRatio:* Ratio of the expression of the exon and its upstream exon.

40. *ThreadRatio:* DownThread / UpThread. Set to [10] if dividing by zero.

41. *RepeatID:* A unique Identifer for each Repeat in the genome (left-most coordinate). Can repeat and thus be used for determining one repeat inititating a trancsript in different assemblies.

42. *LIBRARY:* Library from which this repeat-exon interaction was calculated from.

### `.rslions` & `.inv.rs.lions`
The `recurrant` and `specific` TE-initiations from the biological grouping of libraries
as defined in `input.list` (Set 2 vs. Set 1). This is a collapsed table of `.lions` 
in which common TE-initiation events are collapsed into a single entry.

The `inv.rs.lions` file is the equivalent inverse analysis (Set 1 vs. Set 2).

11. *Normal_occ:* Number of times this TE-initiation was found in the "normal" set of libraries (Set #1)

12. *Cancer_occ:* Number of times this TE-initiation was found in the "cancer" set of libraries (Set #2)

13. *Library:* A semi-colon seperated list of the LIBRARY identifiers in which this TE-initiation was found in.